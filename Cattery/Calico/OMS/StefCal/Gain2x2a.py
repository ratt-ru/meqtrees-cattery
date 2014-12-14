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

from DataTiler import DataTiler    

class Gain2x2a (DataTiler):
  """Support class to handle a set of subtiled gains in the form of 2x2 G matrices.
  Version 2x2a: add proper flag accounting during smoothing.
  """;
  polarized = True;
  nparm = 4;

  def __init__ (self,original_datashape,datashape,subtiling,solve_ifrs,opts,
                init_value=1,verbose=0,
                force_subtiling=False,**kw):
    """original_datashape gives the unpadded datashape. This is the one that ought to be used to estimate
    convergence quotas. datashape is the real, padded, datashape. subtiling gives the subtiling."""
    # init base datatiler class
    DataTiler.__init__(self,datashape,subtiling,original_datashape=original_datashape,force_subtiling=force_subtiling);
    _verbosity.set_verbose(verbose);
    _verbosity.enable_timestamps(True,modulo=6000);
    self._solve_ifrs = solve_ifrs;
    self.opts = opts;
    # init empty parms
    self._antennas = set();
    for ifr in self._solve_ifrs:
      self._antennas.update(ifr);
    self._zero  = numpy.zeros(self.subshape,dtype=complex);
    self._unity = numpy.ones(self.subshape,dtype=complex);
    self._nullflag =  numpy.zeros(self.subshape,dtype=bool);
    self.gainflags = {};
    # init_value=1: init each parm with a unity Jones term
    # if init_value is a dict, use it to initialize each array with a different value
    # as a starting point
    if isinstance(init_value,dict):
      # init values will write into arrays, so make a copy
      self.gain = dict([(p,[self._unity.copy(),self._zero.copy(),self._zero.copy(),self._unity.copy()]) for p in self._antennas]);
      for p,values in init_value.iteritems():
        gmat = self.gain.get(p);
        if gmat is not None:
          for i,(g,value) in enumerate(zip(gmat,values)):
            if is_null(value):
              gmat[i] = 0;
            elif numpy.isscalar(value) or value.ndim == 1:
              g[numpy.newaxis,...] = value;
            else:
              # just being defensive if shapes are different 
              slc = tuple([ slice(0,min(a,b)) for a,b in zip(g.shape,value.shape) ]);
              g[slc] = value[slc];
    # else assume scalar init value, and use it to initialize default array
    # subsequent steps create new arrays, so sufficient to use the same initial value object for all antennas
    else:
##      default = numpy.empty(self.subshape,dtype=complex);
##      default[...] = init_value;
      self.gain = dict([ (p,(init_value,0,0,init_value)) for p in self._antennas ]);
##      self.gain = dict([ (p,(default,self._zero,self._zero,default)) for p in self._antennas ]);
    # setup convergence targets
    self.convergence_target = round(self.real_slots*opts.convergence_quota);
    self._reset();
    dprint(1,"convergence target %d of %d real slots"%(self.convergence_target,self.real_slots));

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

  def iterate (self,lhs,rhs,bitflags,bounds=None,verbose=0,niter=0,weight=None):
    self._reset();
    # iterates G*lhs*G^H -> rhs
    # G updates from step 0 and 1 go here
    gain0dict = {};
    gain1dict = {};
    nflag = 0;
    flags_per_antenna = {};
    gaindiff2 = {};
    for step,(gain0,gain1) in enumerate([(self.gain,gain0dict),(gain0dict,gain1dict)]):
      active_gain = gain0.copy() if self.opts.omega is not None else gain0;
      # loop over all antennas
      for p in self._antennas:
        pmask = self.gainflags.get(p,False);
        g0p = gain0[p];
        # build up sums
        sum_dv = NULL_MATRIX();
        sum_vhv = NULL_MATRIX();
        for q in self._antennas:
          if (p,q) in self._solve_ifrs or (q,p) in self._solve_ifrs:
            # get D/D^H and M/M^H, depending on p<q or p>q
            m = self._get_conj_matrix(p,q,lhs);
            d = self._get_matrix(p,q,rhs);
            if m is None or d is None:
              continue;
            if weight is not None:
              ww = weight.get((p,q),weight.get((q,p),None));
              if is_null(ww):
                continue;
              m = matrix_scale(m,ww);
              d = matrix_scale(d,ww);
            # get applicable flags -- these are same shape as gains
            pqmask = pmask|self.gainflags.get(q,False);
            # get bitflag mask -- same shape as the data
            bfmask = bitflags.get((p,q),0)!=0;
            # zero flagged elements in data, model
            if numpy.any(bfmask):
              bfmask = self.tile_data(bfmask);
              # copy matrices, unless we applied weight above, in which case we already have the copy
              if weight is None:
                d = matrix_copy(d);
                m = matrix_copy(m);
              for idm,dm in enumerate((d,m)):
                for x in dm:
                  if not is_null(x):
#                    print ('m' if idm else 'd'),p,q,x.shape,bfmask.shape
                    x[bfmask] = 0;
            # get the current gain
            g = map(self.tile_subshape,active_gain[q]);
#            print p,q,niter,step,": M",[ is_null(x) for x in m ],"D",[ is_null(x) for x in d ],"G",[ is_null(x) for x in g ];
            # multiply and accumulate
            v = matrix_multiply(g,m);
            dv  = matrix_multiply(d,v);
            vhv = matrix_multiply(matrix_conj(v),v);
#            print "V",[ is_null(x) for x in v ],"DV",[ is_null(x) for x in dv ],"VHV",[ is_null(x) for x in vhv ];
#            print m[1],d[1],v[1];
            if (p,q) in verbose_baselines:
              print "%s%s"%(p,q),"D",[ g[verbose_element] for g in d[0],d[3] ],"M",[ g[verbose_element] for g in m[0],m[3] ];
              print "%s%s"%(p,q),"Gq",[ g[verbose_element] for g in gain0[q][0],gain0[q][3] ],"V",[ g[verbose_element] for g in v[0],v[3] ];
              print "%s%s"%(p,q),"DV",[ g[verbose_element] for g in dv[0],dv[3] ],"VHV",[ g[verbose_element] for g in vhv[0],vhv[3] ];
            # reduce tiles back to gain shape
            dv1 = map(self.reduce_tiles,dv);
            vhv1 = map(self.reduce_tiles,vhv);
            # mask out flagged elements
            if numpy.any(pqmask):
              for mat in dv1,vhv1:
                for x in mat:
                  if not is_null(x):
                    x[pqmask] = 0;
            # add
            matrix_add1(sum_dv,dv1);
            matrix_add1(sum_vhv,vhv1);
        if all([is_null(x) for x in sum_vhv]):
          gain1[p] = g0p;
          continue;
#        print p,q,niter,step,": SDV",[ is_null(x) for x in sum_dv ],"SVHV",[ is_null(x) for x in sum_vhv ];
        if self.opts.real_only:
          dv  = [ x.real if not is_null(x) else x for x in sum_dv ];
          vhv = [ x.real if not is_null(x) else x for x in sum_vhv ];
        # smooth with gaussian, if enabled
        if self.opts.smoothing:
          sum_vhv = [ x if is_null(x) else 
                        scipy.ndimage.filters.gaussian_filter(x.real,self.opts.smoothing,mode='constant')
                        +1j*scipy.ndimage.filters.gaussian_filter(x.imag,self.opts.smoothing,mode='constant') 
                      for x in sum_vhv ];
          sum_dv = [ x if is_null(x) else
                        scipy.ndimage.filters.gaussian_filter(x.real,self.opts.smoothing,mode='constant')
                        +1j*scipy.ndimage.filters.gaussian_filter(x.imag,self.opts.smoothing,mode='constant') 
                     for x in sum_dv ];
        # invert and do update
        inv_vhv = matrix_invert(sum_vhv);
        g1p = gain1[p] = matrix_multiply(sum_dv,inv_vhv);
#        print p,q,niter,step,": IVHV",[ is_null(x) for x in inv_vhv ],"G1P",[ is_null(x) for x in g1p ];
        
        if p in verbose_stations:
#            print "S%d"%step,p,"sum DV",[ g[verbose_element] for g in sum_dv ],"sum VHV",[ g[verbose_element] for g in sum_vhv ];
          print "S%d"%step,p,"G'"," ".join([ "%.5g@%.3g"%(abs(g[verbose_element]),cmath.phase(g[verbose_element])) for g in g1p ]);
          print "S%d"%step,p,"G'"," ".join([ "%.7g%+.4g"%(g[verbose_element].real,g[verbose_element].imag) for g in g1p ]);
        # take mean with previous value
        if self.opts.average == 1 or (self.opts.average == 2 and step):
          if self.opts.omega == 0.5:
            matrix_scale1(matrix_add1(g1p,g0p),0.5);
          else:
            matrix_add1(matrix_scale1(g1p,self.opts.omega),matrix_scale(g0p,(1-self.opts.omega)));
          if p in verbose_stations:
            print "SA",p," ".join([ "%.5g@%.3g"%(abs(g[verbose_element]),cmath.phase(g[verbose_element])) for g in gain1[p] ]);
            print "SA",p," ".join([ "%.7g%+.4g"%(g[verbose_element].real,g[verbose_element].imag) for g in gain1[p] ]);
        # mask out infs/nans
        flag0 = self.gainflags.get(p,False);
        masks = [ (~numpy.isfinite(g1))|flag0 for g1 in g1p ];
        for g1,g0,m in zip(g1p,g0p,masks):
          if not is_null(g1):
            g1[m] = (g0[m] if not numpy.isscalar(g0) else g0);
        # compute norm of gain diff
        if step:
          square = lambda x:(x*numpy.conj(x)).real;
          gaindiff2[p] = gd2 = sum([ square(g1-g0) for g0,g1 in zip(g0p,g1p) ]);
          # flag out-of-bounds gains
          if bounds:
            lower,upper = bounds;
            absg = map(abs,g1p);
            # oob[i] will be True if gain element #i is out of bounds, and not masked above
            oob = [ ((x<(lower or 0))|(upper and x>upper)) for x in (absg[0],absg[3]) ];
            # reduce to single mask for all 4 gains
            flag = reduce(operator.or_,oob);
            nfl = flag.sum();
            if nfl:
              flags_per_antenna[p] = nfl;
              # add mask to gain flags
              if p in self.gainflags:
                self.gainflags[p] |= flag;
              else:
                self.gainflags[p] = flag;
              nflag += nfl;
              # reset gains to unity and diff to zero
              gd2[flag] = 0;
              for g0,g1,default in zip(self.gain[p],g1p,(1,0,0,1)):
                if not numpy.isscalar(g0):
                  g0[flag] = 1;
                if not numpy.isscalar(g1):
                  g1[flag] = 1;
        # feed forward G', if enabled
        if self.opts.feed_forward:
          active_gain[p] = g1p;

    if flags_per_antenna:
      dprint(3,"new gain-flags: "," ".join(["%s:%d"%x for x in sorted(flags_per_antenna.iteritems())]));

    deltanorm_sq = sum(gaindiff2.itervalues());
    gainnorm_sq  = sum([ sum([ square(g1) for g1 in gain1 ]) for p,gain1 in gain1dict.iteritems() ]);
    self.gainnorm = numpy.sqrt(gainnorm_sq).max();
    # find how many have converged
    self.delta_sq = deltanorm_sq/gainnorm_sq;
    self.delta_sq[gainnorm_sq==0] = 0;
    self.converged_mask = self.delta_sq <= self.opts.epsilon**2;
    self.num_converged = self.converged_mask.sum() - self.padded_slots;
    self.delta_max = math.sqrt(self.delta_sq.max());
    self.gain = gain1dict;
    return (self.num_converged >= self.convergence_target),self.delta_max,self.delta_sq,nflag;
    
  def get_converged_mask (self):
    """Returns mask (of shape datashape) showing which t/f slots have been converged""";
    cmask = numpy.zeros(self.datashape,bool);
    self.tile_data(cmask)[...] = self.tile_subshape(self.converged_mask);
    return cmask;

  def get_gainflags (self,flagnull=False,lower=None,upper=None):
    """works out gain flags based on upper and lower bounds (also taking into account previously raised flags)
    Returns array of gain-based flags of shape datashape.""";
    for p,g0p in self.gain.iteritems():
      mask = False;
      # accumulate mask
      for g0 in g0p[0],g0p[3]:
        if flagnull:
          mask |= (g0==0);
        if lower:
          mask |= (abs(g0)<lower);
        if upper:
          mask |= (abs(g0)>upper);
      # got anything?
      if numpy.any(mask):
        if p in self.gainflags:
          self.gainflags[p] |= mask;
        else:
          self.gainflags[p] = mask;
    # reformat entries to datashape
    outflags = {};
    for p,gf in self.gainflags.iteritems():
      fl = numpy.zeros(self.datashape,bool);
      self.tile_data(fl)[...] = self.tile_subshape(gf);
      outflags[p] = fl;
    return outflags;

  def _G (self,p):
    g = self._gmat.get(p);
    if g is None:
      g = self._gmat[p] = map(self.tile_subshape,self.gain[p]);
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
      gi = self._ginv[p,reg] = matrix_maskinf(matrix_invert((a+reg,b,c,d+reg)));
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

  def residual (self,lhs,rhs,pq,tiler=None,cache=True):
    """Returns residual R = Gp*lhs*conj(Gq) - rhs, tiled into subtile shape.
    Computes it on-demand, if not already cached""";
    cache = cache and tiler;
    r = self._residual_cache.get(pq);
    if not cache or r is None:
      c = self.apply(lhs,pq,cache=cache,tiler=tiler);
      r = matrix_sub(c,rhs[pq]);
      if cache:
        self._residual_cache[pq] = r;
      if pq in verbose_baselines_corr:
        print pq,"D",[ 0 if is_null(g) else g[verbose_element] for g in lhs[pq] ];
        print pq,"M",[ 0 if is_null(g) else g[verbose_element] for g in rhs[pq] ];
        print pq,"C",[ 0 if is_null(g) else g[verbose_element] for g in c ];
        print pq,"R",[ 0 if is_null(g) else g[verbose_element] for g in r ];
    else:
      if pq in verbose_baselines_corr:
        print pq,"R cached",[ 0 if is_null(g) else g[verbose_element] for g in r ];
    return r;

  def residual_inverse (self,lhs,rhs,pq,regularize=0,tiler=None,cache=True):
    """Returns residual R = Gp^{-1}*lhs*Gq^{-1H} - rhs, tiled into subtile shape.
    Computes it on-demand, if not already cached""";
    cache = cache and tiler;
    r = self._residual_inverse_cache.get(pq);
    if not cache or r is None:
      c = self.apply_inverse(lhs,pq,cache=cache,regularize=regularize,tiler=tiler);
      r = matrix_sub(c,rhs[pq]);
      if cache:
        self._residual_inverse_cache[pq] = r;
      if pq in verbose_baselines_corr:
        print pq,"D",[ 0 if is_null(g) else g[verbose_element] for g in lhs[pq] ];
        print pq,"M",[ 0 if is_null(g) else g[verbose_element] for g in rhs[pq] ];
        print pq,"C",[ 0 if is_null(g) else g[verbose_element] for g in c ];
        print pq,"Ri",[ 0 if is_null(g) else g[verbose_element] for g in r ];
    return r;

  def apply (self,lhs,pq,cache=False,tiler=False):
    """Returns G*lhs*Gq^H."""
    p,q = pq;
    cache = cache and tiler;
    tiler = tiler or self;
    appl = self._apply_cache.get(pq) if cache else None;
    if appl is None:
      lhs = self._get_matrix(p,q,lhs);
      appl = map(tiler.untile_data,matrix_multiply(self._G(p),matrix_multiply(lhs,self._Gconj(q))));
      if cache:
        self._apply_cache[pq] = appl;
    return appl;

  def apply_inverse (self,rhs,pq,cache=False,regularize=0,tiler=None):
    """Returns corrected data Gp^{-1}*D*Gq^{H-1}."""
    p,q = pq;
    cache = cache and tiler;
    tiler = tiler or self;
    appl = self._apply_inverse_cache.get(pq) if cache else None;
    if appl is None:
      appl = map(tiler.untile_data,matrix_multiply(self._Ginv(p,regularize),
        matrix_multiply(map(tiler.tile_data,rhs[pq]),self._Ginvconj(q,regularize))));
      if cache:
        self._apply_inverse_cache[pq] = appl;
      if pq in verbose_baselines_corr:
        print pq,"RHS",[ 0 if is_null(g) else g[verbose_element] for g in rhs[pq] ];
        print pq,"APPLINV",[ 0 if is_null(g) else g[verbose_element] for g in appl ];
    return appl;

  def get_2x2_gains (self,datashape,expanded_dataslice,tiler=None):
    tiler = tiler or self;
    gainsdict = {};
    for p,gmat in self.gain.iteritems():
      gainsdict[p] = [ tiler.expand_subshape(x,datashape,expanded_dataslice) for x in gmat ];
    return gainsdict;

  def check_finiteness (self,data,label,complete=False):
    # check for NANs in the data
    for i,d in enumerate(data):
      if type(d) is numpy.ndarray:
        fin = numpy.isfinite(d);
        if not fin.all():
          dprintf(0,"%s %s: %d slots are INF/NAN\n",label,i,d.size-fin.sum());
          if not complete:
            break;


