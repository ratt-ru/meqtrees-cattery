# -*- coding: utf-8 -*-
import numpy
import math
import scipy.ndimage.filters

from MatrixOps import *

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
  def __init__ (self,datashape,subtiling,solve_ifrs,epsilon,conv_quota,init_value=1,smoothing=None):
    self._solve_ifrs = solve_ifrs;
    self._epsilon = epsilon;
    self.datashape = datashape;
    self.subtiling = subtiling;
    self.gainshape = tuple([ nd/nt for nd,nt in zip(datashape,subtiling) ]);
    self.smoothing = smoothing;
    # if subtiling is 1,1,... then override methods with identity relations
    if max(subtiling) == 1:
      self.tile_data = self.untile_data = self.tile_gain = self.reduce_subtiles = identity_function;
    else:
      self.subtiled_shape = [];
      self.gain_expansion_slice = [];
      for ng,nt in zip(self.gainshape,subtiling):
        self.subtiled_shape += [ng,nt];
        self.gain_expansion_slice += [slice(None),numpy.newaxis];
      self.subtiled_axes = range(1,len(datashape)*2,2)[-1::-1];
    # init empty parms
    self._parms = parms = set();
    for ifr in self._solve_ifrs:
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
    self._residual_cache = {};
    self._apply_cache = {};
    self._apply_inverse_cache = {};
    self._gpgq = {};

  # define methods
  def tile_data (self,x):
    return x.reshape(self.subtiled_shape) if not numpy.isscalar(x) else x;
  def untile_data (self,x):
    return x.reshape(self.datashape) if not numpy.isscalar(x) else x;
  def tile_gain (self,x):
    return x[self.gain_expansion_slice] if not numpy.isscalar(x) else x;
  def reduce_subtiles (self,x):
    if not numpy.isscalar(x):
      for ax in self.subtiled_axes:
        x = x.sum(ax);
    return x;

  def iterate (self,lhs,rhs,first_iter=False,verbose=0):
    """Does one iteration of Gp*lhs*Gq^H -> rhs""";
    self._residual_cache = {};
    self._apply_cache = {};
    self._apply_inverse_cache = {};
    self._gpgq = {};
    self._gpgq_inv = {};
    self._gp_inv = {};
    gain1 = {};  # dict of gain parm updates from this step
    #
#    if verbose:
#      print [ (rhs[key].min(),lhs[key].min()) for key in self._solve_ifrs ];
    # converge from one side (solve for Gp)
    for p,i in self.gain.keys():
      # build up sums
      sum_reim = numpy.zeros(self.gainshape,dtype=complex);
      sum_sq = 0;
      for q,j in self.gain.keys():
        pq,direct,conjugate = pq_direct_conjugate(p,q,rhs);
        if pq in self._solve_ifrs:
          m,d = conjugate(lhs[pq][i*2+j]),direct(rhs[pq][i*2+j]);
          mh = self.tile_data(m)*self.tile_gain(self.gain[q,j]);
          dmh = self.tile_data(d)*mh;
          mh2 = abs(mh)**2;
#          if p == '0' and i==1:
#            print "S0 %s%s:%s%s"%(p,q,i,j),"VHV",mh2[verbose_element];
          if (p,q) in verbose_baselines and i==j:
            print "%s%s:%s%s"%(p,q,i,j),"D",d[verbose_element],"MH",m[verbose_element], \
                "Gq",self.gain[q,j][verbose_element],"V",mh[verbose_element],"DV",dmh[verbose_element],"VHV",mh2[verbose_element];
          sum_reim += self.reduce_subtiles(dmh);
          sum_sq += self.reduce_subtiles(mh2);
      # generate update
      if self.smoothing:
        sum_sq = scipy.ndimage.filters.gaussian_filter(sum_sq.real,self.smoothing,mode='constant');
        sum_reim.real = scipy.ndimage.filters.gaussian_filter(sum_reim.real,self.smoothing,mode='constant');
        sum_reim.imag = scipy.ndimage.filters.gaussian_filter(sum_reim.imag,self.smoothing,mode='constant');
      if verbose:
        print p,i,sum_sq.min(),sum_sq.max();
#      if p == '0' and i==1:
#        print "S0 %s:%s"%(p,i),"sum VHV",sum_sq[verbose_element];
      if p in verbose_stations:
        print "S0 %s:%s"%(p,i),"sum DV",sum_reim[verbose_element],"sum VHV",sum_sq[verbose_element];
        print "S0 %s:%s"%(p,i),"G'",(sum_reim/sum_sq),(sum_reim/sum_sq)[verbose_element];
      mask = sum_sq==0;
#      gain1[p,i] = (sum_reim/sum_sq + self.gain[p,i])/2;
      gain1[p,i] = sum_reim/sum_sq;
      gain1[p,i][mask] = self.gain[p,i][mask];
      if p in verbose_stations:
        print "S0 %s:%s"%(p,i),"G''",gain1[p,i],gain1[p,i][verbose_element];
#    print "step 0 G:0",gain1['0',0][verbose_element],gain1['0',1][verbose_element];
    if verbose:
      print "done left iter";
    gain2 = {};  # dict of gain parm updates from this step
    # then solve for Gq
    for q,j in self.gain.keys():
      # build up sums
      sum_reim = numpy.zeros(self.gainshape,dtype=complex);
      sum_sq = 0;
      for p,i in self.gain.keys():
        pq,direct,conjugate = pq_direct_conjugate(p,q,rhs);
        if pq in self._solve_ifrs:
          m,d = direct(lhs[pq][i*2+j]),conjugate(rhs[pq][i*2+j]);
          mh = self.tile_data(m)*self.tile_gain(gain1[p,i]);
          dmh = self.tile_data(d)*mh;
          mh2 = abs(mh)**2;
#          if q == '0' and j==1:
#            print "S1 %s%s:%s%s"%(p,q,i,j),"VHV",mh2[verbose_element];
          if (p,q) in verbose_baselines and i==j:
#          if q == 'A' and i==j and j==1:
            print "%s%s:%s%s"%(p,q,i,j),"DH",d[verbose_element],"M",m[verbose_element], \
                "Gp",gain1[p,i][verbose_element],"V",mh[verbose_element],"DHV",dmh[verbose_element],"VHV",mh2[verbose_element];
          sum_reim += self.reduce_subtiles(dmh);
          sum_sq += self.reduce_subtiles(mh2);
          # smoothing?
#      if q == '0' and j==1:
#        print "S1 %s:%s"%(q,j),"sum VHV",sum_sq[verbose_element];
      # generate update
      if self.smoothing:
        sum_sq = scipy.ndimage.filters.gaussian_filter(sum_sq.real,self.smoothing,mode='constant');
        sum_reim.real = scipy.ndimage.filters.gaussian_filter(sum_reim.real,self.smoothing,mode='constant');
        sum_reim.imag = scipy.ndimage.filters.gaussian_filter(sum_reim.imag,self.smoothing,mode='constant');
      if verbose:
        print q,j,sum_sq.min(),sum_sq.max();
      if q in verbose_stations:
        print "S1 %s:%s"%(q,j),"sum DV",sum_reim[verbose_element],"sum VHV",sum_sq[verbose_element];
        print "S1 %s:%s"%(q,j),"G'",(sum_reim/sum_sq),(sum_reim/sum_sq)[verbose_element];
      mask = sum_sq==0;
      gain2[q,j] = (sum_reim/sum_sq + gain1[q,j])/2;
      gain2[q,j][mask] = gain1[q,j][mask];
      if q in verbose_stations:
        print "S1 %s:%s"%(q,j),"G''",gain2[q,j],gain2[q,j][verbose_element];
    # check for convergence
    if verbose:
      print "done right iter";
#    print "step 1 G:0",gain2['0',0][verbose_element],gain2['0',1][verbose_element];
    # compute ||G_new-G_old|| and ||G_new|| for each time/freq slot
    deltanorm_sq = sum([(gain2[pp] - self.gain[pp])**2 for pp in self.gain.iterkeys()]);
    gainnorm_sq = sum([gain2[pp]**2 for pp in self.gain.iterkeys()]);
    # find how many have converged
    self.num_converged = (deltanorm_sq/gainnorm_sq <= self._epsilon**2).sum();
    self.maxdiff = numpy.sqrt(deltanorm_sq.max());
    self.gain = gain2;
    return self.num_converged >= self.convergence_target;

  def gpgq (self,pq,i,j):
    """Returns Gp*conj(Gq), tiled into subtile shape.
    Computes it on-demand, if not already cached""";
    g = self._gpgq.get((pq,i,j));
    if g is None:
      g = self._gpgq[pq,i,j] = self.tile_gain( self.gain.get((pq[0],i),self._unity)*
                                                numpy.conj(self.gain.get((pq[1],j),self._unity)) );
    return g;

  def gpgq_inv (self,pq,i,j,reg=0):
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

  def apply (self,lhs,pq,index=True,cache=False):
    """Returns lhs with gains applied: Gp*lhs*Gq^H."""
    if index:
      lhs = lhs[pq];
    corr = self._apply_cache.get(pq) if cache else None;
    if corr is None:
      corr = [ 0 if is_null(d) else self.untile_data(self.tile_data(d)*self.gpgq(pq,i,j))
                                   for d,(i,j) in zip(lhs,IJ2x2) ];
      if cache:
        self._apply_cache[pq] = corr;
      if pq in verbose_baselines_corr:
        print pq,"UNCORR",[ 0 if is_null(g) else g[verbose_element] for g in lhs ];
        print pq,"CORR",[ 0 if is_null(g) else g[verbose_element] for g in corr ];
    return corr;


  def apply_inverse (self,rhs,pq,cache=False):
    """Returns rhs with inverse gains applied: Gp^{-1}*rhs*Gq^{H-1}."""
    corr = self._apply_inverse_cache.get(pq) if cache else None;
    if corr is None:
      mod = rhs[pq];
      corr = [ 0 if is_null(m) else self.untile_data(self.tile_data(m)*self.gpgq_inv(pq,i,j)) for m,(i,j) in zip(mod,IJ2x2) ];
      if cache:
        self._apply_inverse_cache[pq] = corr;
    return corr;
    
#  def get_2x2_gains (self):
#    for 
#    return dict([(pq,self.tile_gain(self.gain.get((pq[0],i),self._unity) ]);
    

