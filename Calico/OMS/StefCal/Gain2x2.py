# -*- coding: utf-8 -*-
import numpy
import math
import cmath
import operator
import scipy.ndimage.filters
import Kittens.utils

from MatrixOps import *

_verbosity = Kittens.utils.verbosity(name="gain2x2");
dprint = _verbosity.dprint;
dprintf = _verbosity.dprintf;


verbose_baselines = (); # set([('CS001HBA0','CS003HBA0'),('CS003HBA0','CS001HBA0')]);
#verbose_baselines_corr = set(); # set([('0','5'),('5','0')]);
verbose_baselines_corr = (); set([('CS001HBA0','CS003HBA0'),('CS003HBA0','CS001HBA0')]);
verbose_element = 0,0;

verbose_stations = (); # set(('CS001HBA0','CS001HBA1'));
#verbose_stations = set([p[0] for p in verbose_baselines]+[p[1] for p in verbose_baselines]);
verbose_stations_corr = (); # set([p[0] for p in verbose_baselines_corr]+[p[1] for p in verbose_baselines_corr]);

class Gain2x2 (object):
  """Support class to handle a set of subtiled gains in the form of diagonal G matrices.
    For a data shape of N1,N2,..., and a subtiling of M1,M2,..., we have the following shapes in play:
      gainparm_shape = N1/M1,N2/M2,... Let's call this K1,K2,...
    Then, data is reshaped into
      subtiled_shape =  K1,M1,K2,M2,...
    to which we can apply parms with an index of [:,numpy.newaxis,:,numpy.neqaxis] etc.,
    and then collapse every second axis.
    Thus we define the following methods:
      tile_data(x):        reshapes datashape into subtiled_shape
      untile_data(x):      reshapes subtiled_shape into datashape
      tile_gain(x):        reshapes gain into subtiled shape (by inserting a numpy.newaxis for every second axis)
      reduce_subtiles(x):  reduces subtiled_shape to gain_shape by summing every second axis

    If the common case of a 1,1,... subtiling, all these can be an identity
  """;
  polarized = True;
  nparm = 4;

  def __init__ (self,original_datashape,datashape,subtiling,solve_ifrs,epsilon,conv_quota,
                twosided=False,verbose=0,
                smoothing=None,bounds=None,init_value=1,**kw):
    """original_datashape gives the unpadded datashape. This is the one that ought to be used to estimate
    convergence quotas. datashape is the real, padded, datashape. subtiling gives the subtiling."""
    _verbosity.set_verbose(verbose);
    _verbosity.enable_timestamps(True,modulo=6000);
    self._solve_ifrs = solve_ifrs;
    self._epsilon = epsilon;
    self._twosided = twosided;
    self.datashape = datashape;
    self.subtiling = subtiling;
    self.smoothing = smoothing;
    self.gainshape = tuple([ nd/nt for nd,nt in zip(datashape,subtiling) ]);
    self._bounds = (min(bounds),max(bounds)) if bounds else None;
    # if subtiling is 1,1,... then override methods with identity relations
    if max(subtiling) == 1:
      self.tile_data = self.untile_data = self.tile_gain = self.reduce_subtiles = identity_function;
      self.expand_gain_to_datashape = array_to_vells;
    else:
      self.subtiled_shape = [];
      self.gain_expansion_slice = [];
      for ng,nt in zip(self.gainshape,subtiling):
        self.subtiled_shape += [ng,nt];
        self.gain_expansion_slice += [slice(None),numpy.newaxis];
      self.subtiled_axes = range(1,len(datashape)*2,2)[-1::-1];
    # init empty parms
    self._antennas = set();
    for ifr in self._solve_ifrs:
      self._antennas.update(ifr);
    self._zero  = numpy.zeros(self.gainshape,dtype=complex);
    self._unity = numpy.ones(self.gainshape,dtype=complex);
    self._nullflag =  numpy.zeros(self.gainshape,dtype=bool);
    self.gainflags = {};
    # init_value=1: init each parm with the _unity array
    # if init_value is a dict, use it to initialize each array with a different value
    # presumably this happens when going to the next tile -- we use the previous tile's solution
    # as a starting point
    if isinstance(init_value,dict):
      # init values will write into arrays, so make a copy
      self.gain = dict([(p,(self._unity.copy(),self._zero.copy(),self._zero.copy(),self._unity.copy())) for p in self._antennas]);
      for p,values in init_value.iteritems():
        gmat = self.gain.get(p);
        if gmat is not None:
          for g,value in zip(gmat,values):
            if value.ndim == 1:
              g[numpy.newaxis,...] = value;
            else:
              # just being defensive if shapes are different 
              slc = tuple([ slice(0,min(a,b)) for a,b in zip(g.shape,value.shape) ]);
              g[slc] = value[slc];
    # else assume scalar init value, and use it to initialize default array
    # subsequent steps create new arrays, so sufficient to use the same initial value object for all antennas
    else:
      default = numpy.empty(self.gainshape,dtype=complex);
      default[...] = init_value;
      self.gain = dict([ (p,(default,self._zero,self._zero,default)) for p in self._antennas ]);
    # setup various counts and convergence targets
    # total slots being solved for (which can include padding)
    self.total_slots  = reduce(operator.mul,self.gainshape);
    unpadded_gainshape = tuple([ int(math.ceil(nd/float(nt))) for nd,nt in zip(original_datashape,subtiling) ]);
    real_slots = reduce(operator.mul,unpadded_gainshape);
    real_target = round(real_slots*conv_quota);
    self.convergence_target = real_target + (self.total_slots - real_slots);
    self._reset();
    dprint(1,"convergence target %d of %d real slots (%d of %d padded slots)"%(real_target,real_slots,
        self.convergence_target,self.total_slots));

  # define methods
  def tile_data (self,x):
    return x.reshape(self.subtiled_shape) if not (numpy.isscalar(x) or x.size==1) else x;
  def untile_data (self,x):
    return x.reshape(self.datashape) if not (numpy.isscalar(x) or x.size==1) else x;
  def tile_gain (self,x):
    return x[self.gain_expansion_slice] if not (numpy.isscalar(x) or x.size==1) else x;
  def reduce_subtiles (self,x):
    if numpy.isscalar(x) or x.size==1:
      return x;
    for ax in self.subtiled_axes:
      x = x.sum(ax);
    return x;
  def expand_gain_to_datashape (self,x,datashape,expanded_dataslice):
    if expanded_dataslice:
      a = numpy.empty(self.subtiled_shape,dtype=complex);
      a[...] = self.tile_gain(x);
      b = meq.complex_vells(datashape);
      b[...] = self.untile_data(a)[expanded_dataslice];
      return b;
    else:
      a = meq.complex_vells(self.subtiled_shape);
      a[...] = self.tile_gain(x);
      return self.untile_data(a);

  def _reset (self):
    self._apply_cache = {};
    self._apply_inverse_cache = {};
    self._residual_cache = {};
    self._residual_inverse_cache = {};
    self._gmat = {};
    self._ginv = {};
    self._gconj = {};
    self._ginvconj = {};

  def get_values (self):
    return self.gain.copy();

  def set_values (self,gain):
    self._reset();
    self.gain = gain;

  def _get_matrix (self,p,q,data):
    if (p,q) in data:
      return map(self.tile_data,data[p,q]);
    elif (q,p) in data:
      return matrix_conj(map(self.tile_data,data[q,p]));
    else:
      return None;

  def _get_conj_matrix (self,p,q,data):
    if (p,q) in data:
      return matrix_conj(map(self.tile_data,data[p,q]));
    elif (q,p) in data:
      return map(self.tile_data,data[q,p]);
    else:
      return None;

  def iterate (self,lhs,rhs,flagmask,verbose=0,niter=0,weight=None,averaging=None,omega=None,feed_forward=False):
    self._reset();
    # iterates G*lhs*G^H -> rhs
    # G updates from step 0 and 1 go here
    gain0dict = {};
    gain1dict = {};
    nflag = 0;
    flags_per_antenna = {};
    for step,(gain0,gain1) in enumerate([(self.gain,gain0dict),(gain0dict,gain1dict)]):
      active_gain = gain0.copy() if omega is not None else gain0;
      # loop over all antennas
      for p in self._antennas:
        g0p = gain0[p];
        # build up sums
        sum_dv = NULL_MATRIX();
        sum_vhv = NULL_MATRIX();
        for q in self._antennas:
          if (p,q) in self._solve_ifrs or (q,p) in self._solve_ifrs:
            # get D/D^H and M/M^H, depending on p<q or p>q
            if not self._twosided or not step:
              m = self._get_conj_matrix(p,q,lhs);
              d = self._get_matrix(p,q,rhs);
            else:
              m = self._get_matrix(q,p,lhs);
              d = self._get_conj_matrix(q,p,rhs);
            if weight is not None:
              ww = weight.get((p,q)) or weight.get((q,p),0);
            else:
              ww = 1;
            if m is None or d is None or not ww:
              continue;
#            print p,q,weight;
            # multiply and accumulate
            v = vw = matrix_multiply(map(self.tile_gain,active_gain[q]),m);
            if ww != 1:
              vw  = matrix_scale(v,ww);
            dv  = matrix_multiply(d,vw);
            vhv = matrix_multiply(matrix_conj(v),vw);
#           if p == '0':
#             print "%s%s:11"%(p,q),"VHV",vhv[3][verbose_element];
            if (p,q) in verbose_baselines:
              print "%s%s"%(p,q),"D",[ g[verbose_element] for g in d[0],d[3] ],"M",[ g[verbose_element] for g in m[0],m[3] ];
              print "%s%s"%(p,q),"Gq",[ g[verbose_element] for g in gain0[q][0],gain0[q][3] ],"V",[ g[verbose_element] for g in v[0],v[3] ];
              print "%s%s"%(p,q),"DV",[ g[verbose_element] for g in dv[0],dv[3] ],"VHV",[ g[verbose_element] for g in vhv[0],vhv[3] ];
            dv1 = map(self.reduce_subtiles,dv);
            vhv1 = map(self.reduce_subtiles,vhv);
            #if self.smoothing:
              #vhv1 = [ x if is_null(x) else scipy.ndimage.filters.gaussian_filter(x.real,self.smoothing,mode='constant')
                          #for x in vhv1 ];
              #dv1 = [ x if is_null(x) else
                      #scipy.ndimage.filters.gaussian_filter(x.real,self.smoothing,mode='constant')
                      #+1j*scipy.ndimage.filters.gaussian_filter(x.imag,self.smoothing,mode='constant')
                          #for x in dv1 ];
              ## averaging may have smeared dv into invalid slots, so reset these back to zero
              #if flagmask:
                #fmask = flagmask.get((p,q));
                #if fmask is None:
                  #fmask = flagmask.get((q,p));
                #if fmask is not None:
                  #fmask = self.tile_data(fmask);
                  #for x in dv1:
                    #x[fmask] = 0;
                  #for x in vhv1:
                    #x[fmask] = 0;
            matrix_add1(sum_dv,dv1);
            matrix_add1(sum_vhv,vhv1);
#        if p == '0':
#          print "%s:1"%p,"sum VHV",sum_vhv[3][verbose_element];
        # accumulation done, now invert and multiply
        #print p,"SUM VHV:",sum_vhv;
        # if sum is null, then we had no rhs for this station
        if all([is_null(x) for x in sum_vhv]):
          gain1[p] = g0p;
          continue;
        # smooth
        if self.smoothing:
          sum_vhv = [ x if is_null(x) else scipy.ndimage.filters.gaussian_filter(x.real,self.smoothing,mode='constant')
                      for x in sum_vhv ];
          sdv = [];
          for x in sum_dv:
            if is_null(x):
              sdv.append(x);
            else:
              sdv.append( scipy.ndimage.filters.gaussian_filter(x.real,self.smoothing,mode='constant')
                          +1j*scipy.ndimage.filters.gaussian_filter(x.imag,self.smoothing,mode='constant') );
          sum_dv = sdv;
        #
        g1p = gain1[p] = matrix_multiply(sum_dv,matrix_invert(sum_vhv));
        if p in verbose_stations:
#            print "S%d"%step,p,"sum DV",[ g[verbose_element] for g in sum_dv ],"sum VHV",[ g[verbose_element] for g in sum_vhv ];
          print "S%d"%step,p,"G'"," ".join([ "%.5g@%.3g"%(abs(g[verbose_element]),cmath.phase(g[verbose_element])) for g in g1p ]);
          print "S%d"%step,p,"G'"," ".join([ "%.7g%+.4g"%(g[verbose_element].real,g[verbose_element].imag) for g in g1p ]);
        # take mean with previous value
        if averaging == 1 or ( averaging == 2 and step ):
          if omega == 0.5:
            matrix_scale1(matrix_add1(g1p,g0p),0.5);
          else:
            matrix_add1(matrix_scale1(g1p,omega),matrix_scale(g0p,(1-omega)));
          if p in verbose_stations:
            print "SA",p," ".join([ "%.5g@%.3g"%(abs(g[verbose_element]),cmath.phase(g[verbose_element])) for g in gain1[p] ]);
            print "SA",p," ".join([ "%.7g%+.4g"%(g[verbose_element].real,g[verbose_element].imag) for g in gain1[p] ]);
        # mask out infs/nans
        flag0 = self.gainflags.get(p,False);
        masks = [ (~numpy.isfinite(g1))|flag0 for g1 in g1p ];
        for g1,g0,m in zip(g1p,g0p,masks):
          g1[m] = g0[m];
        # flag out-of-bounds gains
        if self._bounds and niter>2:
          absg = map(abs,g1p);
          # oob[i] will be True if gain element #i is out of bounds, and not masked above
          oob = [ ((x<self._bounds[0])|(x>self._bounds[1]))&~m for x,m in zip(absg,masks) ];
          # reduce to single mask for all 4 gains
          flag = reduce(operator.or_,oob);
          nfl = flag.sum();
          if nfl:
            flags_per_antenna[p] = nfl;
            if flag0 is False:
              self.gainflags[p] = flag;
            else:
              flag0 |= flag;
            nflag += nfl;
            # reset gains to unity
            for g in g1p:
              g[flag] = 1;
            # reset model and data to 0
            flag = self.tile_gain(flag);
            for q in self._antennas:
              pq = (p,q) if (p,q) in lhs else (q,p);
              for hs in lhs,rhs:
                xx = hs.get(pq);
                if xx is not None:
                  for x in xx:
                    x1 = self.tile_data(x);
                    if not is_null(x1):
                      x1[flag] = 0;
        # feed forward G', if enabled
        if feed_forward:
          active_gain[p] = g1p;

    ## constrain phases (since we have an inherent phase ambiguity) -- set the sum of the xx phases to 0
    #phi0 = sum([ numpy.angle(gg[0]) for gg in gain1dict.itervalues() ])/len(gain1dict);
    #print "correcting G for sum of phases:",phi0[0];
    #ephi0 = numpy.exp(-1j*phi0);
    #for gg in gain1dict.itervalues():
      #for g in gg:
        #g *= ephi0;
    if flags_per_antenna:
      dprint(3,"new gain-flags: "," ".join(["%s:%d"%x for x in sorted(flags_per_antenna.iteritems())]));

    square = lambda x:(x*numpy.conj(x)).real;
    deltanorm_sq = sum([ sum([ square(g1-g0) for g0,g1 in zip(self.gain[p],gain1) ]) for p,gain1 in gain1dict.iteritems() ]);
    gainnorm_sq  = sum([ sum([ square(g1) for g1 in gain1 ]) for p,gain1 in gain1dict.iteritems() ]);
    self.gainnorm = numpy.sqrt(gainnorm_sq).max();
#    # print per-antenna max
#    print "--- iter %d max updates per-antenna --"%niter;
#    for i,p in enumerate(sorted(self.gain.keys())):
#      print i,p,[ "%.2g"%math.sqrt((square(g1-g0)/gainnorm_sq).max()) for g0,g1 in zip(self.gain[p],gain1dict[p]) ],\
#                [ "%.2g %.2g"%(math.sqrt(square(g1).min()),math.sqrt(square(g1).max())) for g1 in gain1dict[p] ],\
#                [ fl.sum() for fl in self.gainflags[p] ];
    # find how many have converged
    self.delta_sq = deltanorm_sq/gainnorm_sq;
    self.num_converged = (self.delta_sq <= self._epsilon**2).sum();
    self.delta_max = math.sqrt(self.delta_sq.max());
    self.gain = gain1dict;
    # print norm (maye generate norm as a diagnostic?)
    return (self.num_converged >= self.convergence_target),self.delta_max,self.delta_sq,nflag;

  def _G (self,p):
    g = self._gmat.get(p);
    if g is None:
      g = self._gmat[p] = map(self.tile_gain,self.gain[p]);
    return g;

  def _Gconj (self,p):
    g = self._gconj.get(p);
    if g is None:
      g = self._gconj[p] = matrix_conj(self._G(p));
    return g;

  def _Ginv (self,p,reg=0):
    gi = self._ginv.get((p,reg));
    if gi is None:
      a,b,c,d = self._G(p);
      gi = self._ginv[p,reg] = matrix_invert((a+reg,b,c,d+reg));
      if p in verbose_stations_corr:
        print p,"G",[ g[verbose_element] for g in a,b,c,d ];
        print p,"Ginv",[ g[verbose_element] for g in gi ];
#      print "Inverting",p;
    return gi;

  def _Ginvconj (self,p,reg=0):
    gi = self._ginvconj.get((p,reg));
    if gi is None:
      gi = self._ginvconj[p,reg] = matrix_conj(self._Ginv(p,reg));
    return gi;

  def get_last_timeslot (self):
    """Returns dict of p->g, where p is a parm ID, and g is the gain solution for the last timeslot.
    This can be used to initialize new gain objects (for init_value)"""
    return dict([((p,i,j),self.gain[p,i,j][-1,...]) for p in self._antennas for i,j in IJ2x2 ]);

  def residual (self,lhs,rhs,pq):
    """Returns residual R = Gp*lhs*conj(Gq) - rhs, tiled into subtile shape.
    Computes it on-demand, if not already cached""";
    r = self._residual_cache.get(pq);
    if r is None:
      c = self.apply(lhs,pq,cache=True);
      r = self._residual_cache[pq] = matrix_sub(c,rhs[pq]);
      if pq in verbose_baselines_corr:
        print pq,"D",[ 0 if is_null(g) else g[verbose_element] for g in lhs[pq] ];
        print pq,"M",[ 0 if is_null(g) else g[verbose_element] for g in rhs[pq] ];
        print pq,"C",[ 0 if is_null(g) else g[verbose_element] for g in c ];
        print pq,"R",[ 0 if is_null(g) else g[verbose_element] for g in r ];
    else:
      if pq in verbose_baselines_corr:
        print pq,"R cached",[ 0 if is_null(g) else g[verbose_element] for g in r ];
    return r;

  def residual_inverse (self,lhs,rhs,pq,regularize=0):
    """Returns residual R = Gp^{-1}*lhs*Gq^{-1H} - rhs, tiled into subtile shape.
    Computes it on-demand, if not already cached""";
    r = self._residual_inverse_cache.get(pq);
    if r is None:
      c = self.apply_inverse(lhs,pq,cache=True,regularize=regularize);
      r = self._residual_inverse_cache[pq] = matrix_sub(c,rhs[pq]);
      if pq in verbose_baselines_corr:
        print pq,"D",[ 0 if is_null(g) else g[verbose_element] for g in lhs[pq] ];
        print pq,"M",[ 0 if is_null(g) else g[verbose_element] for g in rhs[pq] ];
        print pq,"C",[ 0 if is_null(g) else g[verbose_element] for g in c ];
        print pq,"Ri",[ 0 if is_null(g) else g[verbose_element] for g in r ];
    return r;

  def apply (self,lhs,pq,cache=False):
    """Returns G*lhs*Gq^H."""
    p,q = pq;
    appl = self._apply_cache.get(pq);
    if appl is None or not cache:
      lhs = self._get_matrix(p,q,lhs);
      appl = map(self.untile_data,matrix_multiply(self._G(p),matrix_multiply(lhs,self._Gconj(q))));
      if cache:
        self._apply_cache[pq] = appl;
    return appl;

  def apply_inverse (self,rhs,pq,cache=False,regularize=0):
    """Returns corrected data Gp^{-1}*D*Gq^{H-1}."""
    p,q = pq;
    appl = self._apply_inverse_cache.get(pq);
    if appl is None or not cache:
      appl = map(self.untile_data,matrix_multiply(self._Ginv(p,regularize),
        matrix_multiply(map(self.tile_data,rhs[pq]),self._Ginvconj(q,regularize))));
      if cache:
        self._apply_inverse_cache[pq] = appl;
      if pq in verbose_baselines_corr:
        print pq,"RHS",[ 0 if is_null(g) else g[verbose_element] for g in rhs[pq] ];
        print pq,"APPLINV",[ 0 if is_null(g) else g[verbose_element] for g in appl ];
    return appl;

  def get_2x2_gains (self,datashape,expanded_dataslice):
    gainsdict = {};
    for p,gmat in self.gain.iteritems():
      gainsdict[p] = [ self.expand_gain_to_datashape(x,datashape,expanded_dataslice) for x in gmat ];
    return gainsdict;


