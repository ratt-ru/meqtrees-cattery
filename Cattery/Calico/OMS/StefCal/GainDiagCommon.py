# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

import numpy
import math
import operator
import itertools
import scipy.ndimage.filters

from .MatrixOps import *

from Timba.Meq import meq
import Kittens.utils

_verbosity = Kittens.utils.verbosity(name="gaindiag");
dprint = _verbosity.dprint;
dprintf = _verbosity.dprintf;


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

from .GainDiag import GainDiag    
from .DataTiler import DataTiler    

square = lambda x:(x*numpy.conj(x)).real;

class GainDiagCommon (GainDiag):
  """Support class to handle a set of subtiled gains in the form of diagonal G matrices.
  """;
  polarized = False;
  nparm = 2;

  def __init__ (self,original_datashape,datashape,subtiling,solve_ifrs,opts,
                init_value=1,verbose=0,
                force_subtiling=False,
                **kw):
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
    self._unity = numpy.ones(self.subshape,dtype=complex if not opts.real_only else float);
    # init_value=1: init each parm with the _unity array
    # subsequent steps create new arrays, so sufficient to use the same initial value object for all antennas
    if init_value == 1:
      self.gain = [ self._unity,self._unity ];
    # if init_value is a dict, use it to initialize each array with a different value
    # presumably this happens when going to the next tile -- we use the previous tile's solution
    # as a starting point
    elif isinstance(init_value,dict):
      self.gain = [ self._unity.copy(),self._unity.copy() ];
      for i in range(2):
        value = init_value.get(i);
        if value is not None:
          g = self.gain[i];
          if value.ndim == 1:
            g[numpy.newaxis,...] = value;
          else:
            # just being defensive if shapes are different 
            slc = tuple([ slice(0,min(a,b)) for a,b in zip(g.shape,value.shape) ]);
            g[slc] = value[slc];
    # else assume scalar init value, and use it to initialize default array
    else:
      default = numpy.empty(self.subshape,dtype=complex);
      default[...] = init_value;
      self.gain = [ default, default ];
    # setup gain flags
    self.gainflags = False
    # setup convergence targets
    self.convergence_target = round(self.real_slots*opts.convergence_quota);
    self._reset();
    dprint(1,"convergence target %d of %d real slots"%(self.convergence_target,self.real_slots));

  def _reset (self):
    self._residual_cache = {};
    self._residual_inverse_cache = {};
    self._apply_cache = {};
    self._apply_inverse_cache = {};
    self._gpgq = {};
    self._gpgq_inv = {};
    self._gp_inv = {};

  def get_values (self):
    return  [ g.copy() for g in self.gain ];

  def set_values (self,gain):
    self._reset();
    self.gain = gain;

  def iterate (self,lhs,rhs,bitflags,bounds=None,verbose=0,niter=0,weight=None):
    """Does one iteration of Gp*lhs*Gq^H -> rhs""";
    self._reset();
    gain0 = [0,0];  # dict of gain parm updates from steps 0 and 1
    gain1 = [0,0];
    # pre-averaged differences
    gaindiff2 = 0
    #
    # get gain flag mask -- this will be the same shape as the gains
    pqmask = self.gainflags;
    # this does two iterations at a time              
    for step,g0,g1 in (0,self.gain,gain0),(1,gain0,gain1):
      # converge from one side (solve for Gp)
      for i in range(2):
        # build up sums for this jones element
        sum_reim = 0;
        sum_sq = 0;
        for p,q,j in itertools.product(self._antennas,self._antennas,list(range(2))):
          pq,direct,conjugate = pq_direct_conjugate(p,q,rhs);
          if pq in self._solve_ifrs:
            # get weights, data, model
            ww = weight.get(pq,None) if weight else 1;
            if ww is None:
              continue;
            m,d = conjugate(lhs[pq][i*2+j]),direct(rhs[pq][i*2+j]);
            if is_null(m) or is_null(d):
              continue;
            # get bitflag mask -- same shape as the data
            bfmask = bitflags.get(pq,0)!=0;
            # apply bitflag mask to make zero model/data
            m,d = m*ww,d*ww;
            if numpy.any(bfmask):
              m[bfmask] = 0;
              d[bfmask] = 0;
            # compute update
            mh = self.tile_data(m)*self.tile_subshape(g0[j]);
            dmh = self.tile_data(d)*mh;
            mh2 = abs(mh)**2;
            # average over tiles
            dmh = self.reduce_tiles(dmh);
            mh2 = self.reduce_tiles(mh2);
            # mask out flagged gain elements
            if not is_null(pqmask):
              dmh[pqmask] = 0;
              mh2[pqmask] = 0;
            # sum
            sum_reim += dmh;
            sum_sq += mh2;
        if self.opts.real_only:
          sum_reim = sum_reim.real;
        # generate update
        if self.opts.smoothing:
          sum_sq = scipy.ndimage.filters.gaussian_filter(sum_sq.real,self.opts.smoothing,mode='mirror');
          if self.opts.real_only:
            sum_reim = scipy.ndimage.filters.gaussian_filter(sum_reim,self.opts.smoothing,mode='mirror');
          else:
            sum_reim.real = scipy.ndimage.filters.gaussian_filter(sum_reim.real,self.opts.smoothing,mode='mirror');
            sum_reim.imag = scipy.ndimage.filters.gaussian_filter(sum_reim.imag,self.opts.smoothing,mode='mirror');
        gold = g0[i];
        # null sumsq means no valid data (or gain flagged)
        if is_null(sum_sq):
          gnew = g1[i] = gold;
          mask = True;
        else:
          # null sumsq in some slot means null model, so keep the gain constant there
          # (this is ok -- null model means simply that the slot was flagged)
          gnew = g1[i] = sum_reim/sum_sq;
          mask = sum_sq==0;
          gnew[mask] = gold[mask];
        # inf/nan gains means something else is very wrong, better print a diagnostic
        mask = (~mask)&(~numpy.isfinite(gnew));
        if mask.any():
          gnew[mask] = gold[mask];
          dprint(2,"%s: %d values reset due to INF/NAN"%(p,mask.sum()));
        if p in verbose_stations:
          print("S%d %s:%s"%(step,p,i),"G''",g1[p,i],g1[p,i][verbose_element]);
        # take difference at second step
        gaindiff2 += square(gnew - gold)
        gaindiff2[pqmask] = 0;
        # apply solution averaging
        if self.opts.average == 1 or (self.opts.average == 2 and step):
          gnew += gold;
          gnew /= 2;
    # apply gain flags based on bounds
    num_flagged = 0;
    if bounds:
      lower,upper = bounds;
      mask = False;
      for i,g0 in enumerate(gain1):
        # accumulate mask
        if lower:
          mask |= (abs(g0)<lower);
        if upper:
          mask |= (abs(g0)>upper);
      # got anything?
      if numpy.any(mask):
        for i in range(2):
          gain1[i][mask] = 1;
          self.gain[i][mask] = 1; 
        num_flagged += mask.sum();
        gaindiff2[mask] = 0;
        self.gainflags |= mask;
    # now sum the ||delta-G||^2 over all gains
    deltanorm_sq = gaindiff2;
    # norm-squared of new gain solution, per each t/f slot
    gainnorm_sq = sum(map(square,gain1));
    self.gainnorm = numpy.sqrt(gainnorm_sq).max();
    # find how many have converged
    self.delta_sq = deltanorm_sq/gainnorm_sq;
    self.delta_sq[gainnorm_sq==0] = 0;
    self.converged_mask = self.delta_sq <= self.opts.epsilon**2;
    self.num_converged = self.converged_mask.sum() - self.padded_slots;
    self.delta_max = numpy.sqrt(self.delta_sq.max());
    self.gain = gain1;
    # print norm (maybe generate norm as a diagnostic?)
    return (self.num_converged >= self.convergence_target),self.delta_max,self.delta_sq,num_flagged;
    
  def get_converged_mask (self):
    """Returns mask (of shape datashape) showing which t/f slots have been converged""";
    cmask = numpy.zeros(self.datashape,bool);
    self.tile_data(cmask)[...] = self.tile_subshape(self.converged_mask);
    return cmask;
    
  def get_gainflags (self,flagnull=False,lower=None,upper=None):
    """works out gain flags based on upper and lower bounds (also taking into account previously raised flags)
    Returns array of gain-based flags of shape datashape.""";
    for i,g0 in enumerate(self.gain):
      # accumulate mask
      mask = False;
      if flagnull:
        mask |= (g0==0);
      if lower:
        mask |= (abs(g0)<lower);
      if upper:
        mask |= (abs(g0)>upper);
    # got anything?
    if numpy.any(mask):
      self.gainflags |= mask;
    # reformat entries to datashape
    fl = numpy.zeros(self.datashape,bool);
    self.tile_data(fl)[...] = self.tile_subshape(self.gainflags);
    outflags = {};
    for p in self._antennas:
      outflags[p] = fl;
    return outflags;

  def gpgq (self,pq,i,j):
    """Returns Gp*conj(Gq), tiled into subtile shape.
    Computes it on-demand, if not already cached""";
    g = self._gpgq.get((i,j));
    if g is None:
      g = self._gpgq[i,j] = self.tile_subshape( self.gain[i]*numpy.conj(self.gain[j]) );
    return g;

  def gpgq_inv (self,pq,i,j,regularize=0):
    """Returns 1/((Gp+reg)*conj(Gq+reg)), tiled into subtile shape.
    Computes it on-demand, if not already cached""";
    if (i,j) not in self._gpgq_inv:
      gpi = [0,0];
      for i in range(2):
        x = gpi[i] = self._gp_inv[i] = 1/(self.gain[i]+regularize);
        x[~numpy.isfinite(x)] = 0; 
      for i in range(2):
        for j in range(2):
          self._gpgq_inv[i,j] = self.tile_subshape(gpi[i]*numpy.conj(gpi[j]));
    return self._gpgq_inv[i,j];

  def get_last_timeslot (self):
    """Returns dict of p->g, where p is a parm ID, and g is the gain solution for the last timeslot.
    This can be used to initialize new gain objects (for init_value)"""
    return dict([(key,value[-1,...]) for key,value in enumerate(self.gain) ]);

  def get_2x2_gains (self,datashape,expanded_dataslice,tiler=None):
    tiler = tiler or self;
    gainsdict = {};
    gx = tiler.expand_subshape(self.gain[0],datashape,expanded_dataslice)
    gy = tiler.expand_subshape(self.gain[1],datashape,expanded_dataslice)
    for p in self._antennas:
      gainsdict[p] = [ gx,meq.sca_vells(0),meq.sca_vells(0),gy ];
    return gainsdict;


