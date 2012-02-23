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

class Subtiled2x2Gain (object):
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
  def __init__ (self,datashape,subtiling,solve_ifrs,epsilon,conv_quota,
    smoothing=None,init_value=1):
    self._solve_ifrs = solve_ifrs;
    self._epsilon = epsilon;
    self.datashape = datashape;
    self.subtiling = subtiling;
    self.smoothing = smoothing;
    self.gainshape = tuple([ nd/nt for nd,nt in zip(datashape,subtiling) ]);
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
    self._antennas = set();
    for ifr in self._solve_ifrs:
      self._antennas.update(ifr);
    self._zero  = numpy.zeros(self.gainshape,dtype=complex);
    self._unity = numpy.ones(self.gainshape,dtype=complex);
    self.gain = dict([(p,(self._unity,self._zero,self._zero,self._unity)) for p in self._antennas]);
    # init_value=1: init each parm with the _unity array
    # if init_value is a dict, use it to initialize each array with a different value
    # presumably this happens when going to the next tile -- we use the previous tile's solution
    # as a starting point
    if isinstance(init_value,dict):
      for p,values in init_value.iteritems():
        gmat = self.gain.get(p);
        if gmat is not None:
          for g,value in zip(gmat,values):
            if value.ndim == 1:
              g[numpy.newaxis,...] = value;
            else:
              g[...] = value;
    # else assume scalar init value, and use it to initialize default array
    else:
      default = numpy.empty(self.gainshape,dtype=complex);
      default[...] = init_value;
      self.gain = dict([ (p,(default,self._zero,self._zero,default)) for p in self._antennas ]);
    # setup various counts and convergence targets
    self.total_parms = len(self._antennas)*4*reduce(lambda a,b:a*b,self.gainshape);
    self.convergence_target = int(self.total_parms*conv_quota);
    self._reset();

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

  def _reset (self):
    self._residual = {};
    self._corrupt = {};
    self._gmat = {};
    self._ginv = {};
    self._gconj = {};
    self._ginvconj = {};

  def _get_matrix (self,p,q,data):
    if (p,q) in data:
      return map(self.tile_data,data[p,q]);
    elif (q,p) in data:
      return map(numpy.conj,map(self.tile_data,data[q,p]));
    else:
      return None;

  def _get_conj_matrix (self,p,q,data):
    if (p,q) in data:
      return matrix_conj(map(self.tile_data,data[p,q]));
    elif (q,p) in data:
      return matrix_transpose(map(self.tile_data,data[q,p]));
    else:
      return None;

  def iterate (self,data,model,verbose=0,first_iter=False):
    self._reset();
    # even step:
    #   need to compute V_i = G_i M^H_i
    #   then sum D_i V_i and V^H_i V_i
    # odd step:
    #   need to compute V_i = G_i M_i
    #   then sum D^H_i V_i and V^H_i V_i
    gain0 = self.gain.copy();
    for step in 0,1:
      # G update goes here
      gain1 = {};
      # loop over all antennas
      for p in self._antennas:
        # build up sums
        sum_dv = sum_vhv = NULL_MATRIX;
        for q in self._antennas:
          if (p,q) in self._solve_ifrs or (q,p) in self._solve_ifrs:
            # get D/D^H and M/M^H, depending on p<q or p>q
            if not step:
              m = self._get_conj_matrix(p,q,model);
              d = self._get_matrix(p,q,data);
            else:
              m = self._get_matrix(q,p,model);
              d = self._get_conj_matrix(q,p,data);
            if m is None or d is None:
              continue;
            # multiply and accumulate
            v = matrix_multiply(map(self.tile_gain,gain0[q]),m);
            dv = matrix_multiply(d,v);
            vhv = matrix_multiply(matrix_conj(v),v);
#           if p == '0':
#             print "%s%s:11"%(p,q),"VHV",vhv[3][verbose_element];
            if (p,q) in verbose_baselines:
              print "%s%s"%(p,q),"D",[ g[verbose_element] for g in d[0],d[3] ],"M",[ g[verbose_element] for g in m[0],m[3] ];
              print "%s%s"%(p,q),"Gq",[ g[verbose_element] for g in gain0[q][0],gain0[q][3] ],"V",[ g[verbose_element] for g in v[0],v[3] ];
              print "%s%s"%(p,q),"DV",[ g[verbose_element] for g in dv[0],dv[3] ],"VHV",[ g[verbose_element] for g in vhv[0],vhv[3] ];
            sum_dv = matrix_add(sum_dv,map(self.reduce_subtiles,dv));
            sum_vhv = matrix_add(sum_vhv,map(self.reduce_subtiles,vhv));
#        if p == '0':
#          print "%s:1"%p,"sum VHV",sum_vhv[3][verbose_element];
        # accumulation done, now invert and multiply
        #print p,"SUM VHV:",sum_vhv;
        # if sum is null, then we had no data for this station
        if all([is_null(x) for x in sum_vhv]):
          gain1[p] = gain0[p];
        else:
          # smooth
          if self.smoothing:
            sum_vhv = [ x if is_null(x) else scipy.ndimage.filters.gaussian_filter(x.real,self.smoothing,mode='constant')
                        for x in sum_vhv ];
            for x in sum_dv:
              if not is_null(x):
                x.real = scipy.ndimage.filters.gaussian_filter(x.real,self.smoothing,mode='constant');
                x.imag = scipy.ndimage.filters.gaussian_filter(x.imag,self.smoothing,mode='constant');
          #
          g1 = gain1[p] = matrix_multiply(sum_dv,matrix_invert(sum_vhv));
          if p in verbose_stations:
            print "S%d"%step,p,"sum DV",[ g[verbose_element] for g in sum_dv ],"sum VHV",[ g[verbose_element] for g in sum_vhv ];
            print "S%d"%step,p,"G'",g1[3],[ g[verbose_element] for g in g1 ];
          # take mean with previous value, and mask out infs/nans
          for g,g0 in zip(g1,gain0[p]):
            g += g0;
            g /= 2;
            mask = ~numpy.isfinite(g);
            g[mask] = g0[mask];
      # end of step -- update gains
      gain0 = gain1;
      #print "step",step,"G:0",[ g[verbose_element] for g in gain0['0'] ];
    # check for convergence
    self.gaindiff = {};
    self.num_converged = 0;
    for p in self._antennas:
      self.gaindiff[p] = 0;
      for n,(i,j) in enumerate(IJ2x2):
        g1 = gain0[p][n];
        delta = abs(g1 - self.gain[p][n]);
        self.num_converged += (delta < self._epsilon).sum();
        self.gaindiff[p] = max(self.gaindiff[p],abs(delta).max());
    self.maxdiff = max(self.gaindiff.itervalues());
    self.gain = gain0;
    #print "final",step,"G:0",[ g[verbose_element] for g in self.gain['0'] ];
    return self.num_converged >= self.convergence_target;

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

  def reset_residuals (self):
    self._residual = {};

  def residual (self,data,model,pq):
    """Returns residual R = D - Gp*M*conj(Gq), tiled into subtile shape.
    Computes it on-demand, if not already cached""";
    r = self._residual.get(pq);
    if r is None:
      c = self.corrupt(model,pq,True)
      r = self._residual[pq] = matrix_sub(data[pq],c);
      if pq in verbose_baselines_corr:
        print pq,"D",[ 0 if is_null(g) else g[verbose_element] for g in data[pq] ];
        print pq,"M",[ 0 if is_null(g) else g[verbose_element] for g in model[pq] ];
        print pq,"C",[ 0 if is_null(g) else g[verbose_element] for g in c ];
        print pq,"R",[ 0 if is_null(g) else g[verbose_element] for g in r ];
    return r;

  def correct (self,data,pq,index=True,regularize=0):
    """Returns corrected data Gp^{-1}*D*Gq^{H-1}."""
    p,q = pq;
    data = data[pq] if index else data;
    corr = map(self.untile_data,matrix_multiply(self._Ginv(p,regularize),
      matrix_multiply(map(self.tile_data,data),self._Ginvconj(q,regularize))));
#    print "Correcting",p,q;
    if pq in verbose_baselines_corr:
      print pq,"UNCORR",[ 0 if is_null(g) else g[verbose_element] for g in data ];
      print pq,"CORR",[ 0 if is_null(g) else g[verbose_element] for g in corr ];
    return corr;

  def corrupt (self,model,pq,cache=False):
    """Returns corrupted model Gp*M*Gq."""
    p,q = pq;
    c = self._corrupt.get(pq);
    if c is None or not cache:
      m = self._get_matrix(p,q,model);
      c = map(self.untile_data,matrix_multiply(self._G(p),matrix_multiply(m,self._Gconj(q))));
      if cache:
        self._corrupt[pq] = c;
    return c;

