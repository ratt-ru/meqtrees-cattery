# -*- coding: utf-8 -*-
import numpy
import math
import scipy.ndimage.filters

from MatrixOps import *
from Timba.Meq import meq

verbose_baselines = ()#set([('0','C'),('C','0')]);
verbose_baselines_corr = ()#set([('0','C'),('C','0')]);
verbose_element = 5,0;

verbose_stations = set([p[0] for p in verbose_baselines]+[p[1] for p in verbose_baselines]);
verbose_stations_corr = set([p[0] for p in verbose_baselines_corr]+[p[1] for p in verbose_baselines_corr]);

def pq_direct_conjugate (p,q,data):
  """Helper function:
  Converts p,q into a (p,q) equation key, and conversion functions for direct and conjugate values.
  If p,q is in data, then this is just an identity relation
  If q,p is in data, then swaps components around, and uses conjugate
  """;
  if (p,q) in data:
    return (p,q),identity_function,numpy.conj
  elif (q,p) in data:
    return (q,p),numpy.conj,identity_function;
  else:
    return None,None,None;

class SubtiledDiagGain (object):
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
  def __init__ (self,datashape,subtiling,solve_ifrs,epsilon,conv_quota,init_value=1,bounds=None,smoothing=None):
    self._solve_ifrs = solve_ifrs;
    self._epsilon = epsilon;
    self.datashape = datashape;
    self.subtiling = subtiling;
    self.gainshape = tuple([ nd/nt for nd,nt in zip(datashape,subtiling) ]);
    self.smoothing = smoothing;
    # if subtiling is 1,1,... then override methods with identity relations
    if max(subtiling) == 1:
      self.tile_data = self.untile_data = self.tile_gain = self.reduce_subtiles =  identity_function;
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
    self._parms = parms = set();
    for ifr in self._solve_ifrs:
      self._antennas.update(ifr);
      self._parms.update([(p,i) for p in ifr for i in range(2)]);
    self._unity = numpy.ones(self.gainshape,dtype=complex);
    # init_value=1: init each parm with the _unity array
    if init_value == 1:
      self.gain = dict([ (pp,self._unity) for pp in parms ]);
    # if init_value is a dict, use it to initialize each array with a different value
    # presumably this happens when going to the next tile -- we use the previous tile's solution
    # as a starting point
    elif isinstance(init_value,dict):
      self.gain = dict([ (pp,self._unity) for pp in parms ]);
      for p,value in init_value.iteritems():
        g = self.gain.get(p);
        if g is not None:
          if value.ndim == 1:
            g[numpy.newaxis,...] = value;
          else:
            g[...] = value;
    # else assume scalar init value, and use it to initialize default array
    else:
      default = numpy.empty(self.gainshape,dtype=complex);
      default[...] = init_value;
      self.gain = dict([ (pp,default) for pp in parms ]);
    # setup various counts and convergence targets
    self.total_parms = reduce(lambda a,b:a*b,self.gainshape);
    self.convergence_target = round(self.total_parms*conv_quota);
    self._reset();

  # define methods
  def tile_data (self,x):
    return x.reshape(self.subtiled_shape) if not numpy.isscalar(x) else x;
  def untile_data (self,x):
    return x.reshape(self.datashape) if not numpy.isscalar(x) else x;
  def tile_gain (self,x):
    return x[self.gain_expansion_slice] if not numpy.isscalar(x) else x;
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

  def reduce_subtiles (self,x):
    if not numpy.isscalar(x):
      for ax in self.subtiled_axes:
        x = x.sum(ax);
    return x;

  def _reset (self):
    self._residual_cache = {};
    self._residual_inverse_cache = {};
    self._apply_cache = {};
    self._apply_inverse_cache = {};
    self._gpgq = {};
    self._gpgq_inv = {};
    self._gp_inv = {};

  def iterate (self,lhs,rhs,niter=0,verbose=0):
    """Does one iteration of Gp*lhs*Gq^H -> rhs""";
    self._reset();
    gain0 = {};  # dict of gain parm updates from steps 0 and 1
    gain1 = {};
    #
#    if verbose:
#      print [ (rhs[key].min(),lhs[key].min()) for key in self._solve_ifrs ];
    for step,g0,g1 in (0,self.gain,gain0),(1,gain0,gain1):
      # converge from one side (solve for Gp)
      for p,i in self.gain.keys():
        # build up sums
        sum_reim = numpy.zeros(self.gainshape,dtype=complex);
        sum_sq = 0;
        for q,j in self.gain.keys():
          pq,direct,conjugate = pq_direct_conjugate(p,q,rhs);
          if pq in self._solve_ifrs:
            m,d = conjugate(lhs[pq][i*2+j]),direct(rhs[pq][i*2+j]);
            mh = self.tile_data(m)*self.tile_gain(g0[q,j]);
            dmh = self.tile_data(d)*mh;
            mh2 = abs(mh)**2;
  #          if p == '0' and i==1:
  #            print "S0 %s%s:%s%s"%(p,q,i,j),"VHV",mh2[verbose_element];
            if (p,q) in verbose_baselines and i==j:
              print "S%d %s%s:%s%s"%(step,p,q,i,j),"D",d[verbose_element],"MH",m[verbose_element], \
                  "Gq",g0[q,j][verbose_element],"V",mh[verbose_element],"DV",dmh[verbose_element],"VHV",mh2[verbose_element];
            sum_reim += self.reduce_subtiles(dmh);
            sum_sq += self.reduce_subtiles(mh2);
        # generate update
        if self.smoothing:
          sum_sq = scipy.ndimage.filters.gaussian_filter(sum_sq.real,self.smoothing,mode='constant');
          sum_reim.real = scipy.ndimage.filters.gaussian_filter(sum_reim.real,self.smoothing,mode='constant');
          sum_reim.imag = scipy.ndimage.filters.gaussian_filter(sum_reim.imag,self.smoothing,mode='constant');
        if verbose:
          print step,p,i,sum_sq.min(),sum_sq.max();
  #      if p == '0' and i==1:
  #        print "S0 %s:%s"%(p,i),"sum VHV",sum_sq[verbose_element];
        if p in verbose_stations:
          print "S%d %s:%s"%(step,p,i),"sum DV",sum_reim[verbose_element],"sum VHV",sum_sq[verbose_element];
          print "S%d %s:%s"%(step,p,i),"G'",(sum_reim/sum_sq),(sum_reim/sum_sq)[verbose_element];
        mask = sum_sq==0;
  #      gain1[p,i] = (sum_reim/sum_sq + self.gain[p,i])/2;
        g1[p,i] = sum_reim/sum_sq;
        g1[p,i][mask] = g0[p,i][mask];
        if p in verbose_stations:
          print "S%d %s:%s"%(step,p,i),"G''",g1[p,i],g1[p,i][verbose_element];
    # after two steps, take the average of the previous two solutions, i.e. G(i) and G(i-1)
    for pi,g0 in gain0.iteritems():
      gain1[pi] = (gain1[pi] + g0)/2;
    # compare to self.gain (i.e. G(i-2))
#    print "step 1 G:0",gain2['0',0][verbose_element],gain2['0',1][verbose_element];
    # compute ||G_new-G_old|| and ||G_new|| for each time/freq slot
    square = lambda x:x*numpy.conj(x);
    deltanorm_sq = sum([square(gain1[pp] - self.gain[pp]) for pp in self.gain.iterkeys()]);
    gainnorm_sq = sum([square(gain1[pp]) for pp in self.gain.iterkeys()]);
    # find how many have converged
    self.delta_sq = deltanorm_sq/gainnorm_sq;
    self.num_converged = (self.delta_sq <= self._epsilon**2).sum();
    self.delta_max = numpy.sqrt(self.delta_sq.max());
    self.gain = gain1;
    # print norm (maybe generate norm as a diagnostic?)
    print "%d slots converged, max delta is %f"%(self.num_converged,self.delta_max);
    return (self.num_converged >= self.convergence_target),self.delta_max,self.delta_sq,0;

  def gpgq (self,pq,i,j):
    """Returns Gp*conj(Gq), tiled into subtile shape.
    Computes it on-demand, if not already cached""";
    g = self._gpgq.get((pq,i,j));
    if g is None:
      g = self._gpgq[pq,i,j] = self.tile_gain( self.gain.get((pq[0],i),self._unity)*
                                                numpy.conj(self.gain.get((pq[1],j),self._unity)) );
    return g;

  def gpgq_inv (self,pq,i,j,regularize=0):
    """Returns 1/((Gp+reg)*conj(Gq+reg)), tiled into subtile shape.
    Computes it on-demand, if not already cached""";
    g = self._gpgq_inv.get((pq,i,j));
    if g is None:
      gpgqi = [];
      p,q = pq;
      for pi in (p,i),(q,j):
        gpi = self._gp_inv.get(pi);
        if gpi is None:
          gpi = self._gp_inv[pi] =  1/(self.gain.get(pi,self._unity)+reg);
          if pi[0] in verbose_stations_corr:
            print "G:%s:%s"%pi,self.gain[pi][verbose_element];
            print "Ginv:%s:%s"%pi,gpi[verbose_element];
        gpgqi.append(gpi);
      g = self._gpgq_inv[pq,i,j] = self.tile_gain(gpgqi[0]*numpy.conj(gpgqi[1]));
    return g;

  def gain_keys (self):
    return self._parms;

  def get_last_timeslot (self):
    """Returns dict of p->g, where p is a parm ID, and g is the gain solution for the last timeslot.
    This can be used to initialize new gain objects (for init_value)"""
    return dict([(key,value[-1,...]) for key,value in self.gain.iteritems() ]);

  def reset_residuals (self):
    self._residual = {};

  def residual (self,lhs,rhs,pq):
    """Returns residual R = apply(lhs) - rhs, tiled into subtile shape.
    Computes it on-demand, if not already cached""";
    res = self._residual_cache.get(pq);
    if res is None:
      corr = self.apply(lhs,pq,index=True,cache=True);
      self._residual_cache[pq] = res = matrix_sub(corr,rhs[pq]);
      if pq in verbose_baselines_corr:
        print pq,"D",[ 0 if is_null(g) else g[verbose_element] for g in lhs[pq] ];
        print pq,"M",[ 0 if is_null(g) else g[verbose_element] for g in rhs[pq] ];
        print pq,"C",[ 0 if is_null(g) else g[verbose_element] for g in corr ];
        print pq,"R",[ 0 if is_null(g) else g[verbose_element] for g in res ];
    return res;

  def residual_inverse (self,lhs,rhs,pq,regularize=0):
    """Returns residual R = apply_inverse(lhs) - rhs, tiled into subtile shape.
    Computes it on-demand, if not already cached""";
    res = self._residual_inverse_cache.get(pq);
    if res is None:
      corr = self.apply_inverse(lhs,pq,cache=True,regularize=regularize);
      self._residual_inverse_cache[pq] = res = matrix_sub(corr,rhs[pq]);
      if pq in verbose_baselines_corr:
        print pq,"D",[ 0 if is_null(g) else g[verbose_element] for g in lhs[pq] ];
        print pq,"M",[ 0 if is_null(g) else g[verbose_element] for g in rhs[pq] ];
        print pq,"C",[ 0 if is_null(g) else g[verbose_element] for g in corr ];
        print pq,"R",[ 0 if is_null(g) else g[verbose_element] for g in res ];
    return res;

  def apply (self,lhs,pq,index=True,cache=False):
    """Returns lhs with gains applied: Gp*lhs*Gq^H."""
    if index:
      lhs = lhs[pq];
    appl = self._apply_cache.get(pq) if cache else None;
    if appl is None:
      appl = [ 0 if is_null(d) else self.untile_data(self.tile_data(d)*self.gpgq(pq,i,j))
                                   for d,(i,j) in zip(lhs,IJ2x2) ];
      if cache:
        self._apply_cache[pq] = appl;
      if pq in verbose_baselines_corr:
        print pq,"LHS",[ 0 if is_null(g) else g[verbose_element] for g in lhs ];
        print pq,"APPL(LHS)",[ 0 if is_null(g) else g[verbose_element] for g in appl ];
    return appl;

  def apply_inverse (self,rhs,pq,cache=False,regularize=0):
    """Returns rhs with inverse gains applied: Gp^{-1}*rhs*Gq^{H-1}."""
    appl = self._apply_inverse_cache.get(pq) if cache else None;
    if appl is None:
      mod = rhs[pq];
      appl = [ 0 if is_null(m) else self.untile_data(self.tile_data(m)*
                  self.gpgq_inv(pq,i,j,regularize=regularize)) for m,(i,j) in zip(mod,IJ2x2) ];
      if cache:
        self._apply_inverse_cache[pq] = appl;
    return appl;

  def get_2x2_gains (self,datashape,expanded_dataslice):
    gainsdict = {};
    for p in self._antennas:
      gx = self.gain.get((p,0),self._unity);
      gy = self.gain.get((p,1),self._unity);
      gainsdict[p] = [ self.expand_gain_to_datashape(gx,datashape,expanded_dataslice),meq.sca_vells(0),meq.sca_vells(0),
                       self.expand_gain_to_datashape(gy,datashape,expanded_dataslice) ];
    return gainsdict;


