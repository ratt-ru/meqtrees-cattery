# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

import numpy
import math
import operator
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

from .DataTiler import DataTiler    

square = lambda x:(x*numpy.conj(x)).real;

class GainDiagPhase (DataTiler):
  """Support class to handle a set of subtiled gains in the form of diagonal G matrices with phase-only solutions
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
    self._parms = parms = set();
    for ifr in self._solve_ifrs:
      self._antennas.update(ifr);
      self._parms.update([(p,i) for p in ifr for i in range(2)]);
    self._unity = numpy.ones(self.subshape,dtype=complex if not opts.real_only else float);
    # init_value=1: init each parm with the _unity array
    # subsequent steps create new arrays, so sufficient to use the same initial value object for all antennas
    if init_value == 1:
      self.gain = dict([ (pp,self._unity) for pp in parms ]);
    # if init_value is a dict, use it to initialize each array with a different value
    # presumably this happens when going to the next tile -- we use the previous tile's solution
    # as a starting point
    elif isinstance(init_value,dict):
      self.gain = dict([ (pp,self._unity.copy()) for pp in parms ]);
      for p,value in init_value.items():
        g = self.gain.get(p);
        if g is not None:
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
      self.gain = dict([ (pp,default) for pp in parms ]);
    # setup gain flags
    self.gainflags = {}
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
    return self.gain.copy();

  def set_values (self,gain):
    self._reset();
    self.gain = gain;

  def iterate (self,lhs,rhs,bitflags,bounds=None,verbose=0,niter=0,weight=None):
    """Does one iteration of Gp*lhs*Gq^H -> rhs""";
    self._reset();
    gain0 = {};  # dict of gain parm updates from steps 0 and 1
    gain1 = {};
    # pre-averaged differences
    gaindiff2 = dict([ (p,0) for p in self._antennas ]);
    #
#    if verbose:
#      print [ (rhs[key].min(),lhs[key].min()) for key in self._solve_ifrs ];
#    for step,g0,g1 in [(0,self.gain,gain1)]:
    # this does two iterations at a time              
    for step,g0,g1 in (0,self.gain,gain0),(1,gain0,gain1):
      # converge from one side (solve for Gp)
      for p,i in list(self.gain.keys()):
        pmask = self.gainflags.get(p,False);
        # build up sums
        sum_reim = 0;
        sum_sq = 0;
        for q,j in list(self.gain.keys()):
          pq,direct,conjugate = pq_direct_conjugate(p,q,rhs);
          if pq in self._solve_ifrs:
            # get weights, data, model
            ww = weight.get(pq,None) if weight else 1;
            if ww is None:
              continue;
            m,d = conjugate(lhs[pq][i*2+j]),direct(rhs[pq][i*2+j]);
            if is_null(m) or is_null(d):
              continue;
            # get gain flag mask -- this will be the same shape as the gains
            pqmask = pmask|self.gainflags.get(q,False);
            # get bitflag mask -- same shape as the data
            bfmask = bitflags.get(pq,0)!=0;
            # apply bitflag mask to make zero model/data
            m,d = m*ww,d*ww;
            if numpy.any(bfmask):
              m[bfmask] = 0;
              d[bfmask] = 0;
            # compute update
            mh = self.tile_data(m)*self.tile_subshape(g0[q,j]);
            dmh = self.tile_data(d)*mh;
            mh2 = abs(mh)**2;
            if (p,q) in verbose_baselines and i==j:
              print("S%d %s%s:%s%s"%(step,p,q,i,j),"D",d[verbose_element],"MH",m[verbose_element], \
                  "Gq",g0[q,j][verbose_element],"V",mh[verbose_element],"DV",dmh[verbose_element],"VHV",mh2[verbose_element]);
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
        if p in verbose_stations:
          print("S%d %s:%s"%(step,p,i),"sum DV",sum_reim[verbose_element],"sum VHV",sum_sq[verbose_element]);
          print("S%d %s:%s"%(step,p,i),"G'",(sum_reim/sum_sq),(sum_reim/sum_sq)[verbose_element]);
        gold = g0[p,i];
        # null sumsq means no valid data (or gain flagged)
        if is_null(sum_sq):
          gnew = g1[p,i] = gold;
          mask = True;
        else:
          # null sumsq in some slot means null model, so keep the gain constant there
          # (this is ok -- null model means simply that the slot was flagged)
          gnew = g1[p,i] = sum_reim/sum_sq;
          # reset amplitude to 1
          absval = abs(gnew);
          absnull = absval==0;
          if absnull.any():
            dprint(2,"%s: %d slots have null gain solutions, keeping previous value"%(p,absnull.sum()));
          gnew /= absval;
          mask = (sum_sq==0)|absnull;
          gnew[mask] = gold[mask];
        # inf/nan gains means something else is very wrong, better print a diagnostic
        mask = (~mask)&(~numpy.isfinite(gnew));
        if mask.any():
          gnew[mask] = gold[mask];
          dprint(2,"%s: %d values reset due to INF/NAN"%(p,mask.sum()));
        if p in verbose_stations:
          print("S%d %s:%s"%(step,p,i),"G''",g1[p,i],g1[p,i][verbose_element]);
        # take difference at second step
        gaindiff2[p] += square(gnew - gold)
        gaindiff2[p][mask] = 0;
        # apply solution averaging
        if self.opts.average == 1 or (self.opts.average == 2 and step):
          gnew += gold;
          gnew /= 2;
          gnew /= abs(gnew);
        # apply feed-forward
        if self.opts.feed_forward:
          g0[p,i] = gnew;
    # apply gain flags based on bounds
    num_flagged = 0;
    if bounds:
      lower,upper = bounds;
      for p in self._antennas:
        # accumulate mask
        mask = False;
        for i in range(2):
          g0 = gain1.get((p,i),None);
          if g0 is not None:
            if lower:
              mask |= (abs(g0)<lower);
            if upper:
              mask |= (abs(g0)>upper);
        # got anything?
        if numpy.any(mask):
          for i in range(2):
            gain1[p,i][mask] = 1;
            self.gain[p,i][mask] = 1; 
          num_flagged += mask.sum();
          gaindiff2[p][mask] = 0;
          if p in self.gainflags:
            self.gainflags[p] |= mask;
          else:
            self.gainflags[p] = mask;
    # now sum the ||delta-G||^2 over all gains
    deltanorm_sq = sum(gaindiff2.values());
    # norm-squared of new gain solution, per each t/f slot
    gainnorm_sq = sum([ square(g1) for pp,g1 in gain1.items() ]);
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
    for (p,i),g0 in self.gain.items():
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
        if p in self.gainflags:
          self.gainflags[p] |= mask;
        else:
          self.gainflags[p] = mask;
    # reformat entries to datashape
    outflags = {};
    for p,gf in self.gainflags.items():
      fl = numpy.zeros(self.datashape,bool);
      self.tile_data(fl)[...] = self.tile_subshape(gf);
      outflags[p] = fl;
    return outflags;

  def gpgq (self,pq,i,j):
    """Returns Gp*conj(Gq), tiled into subtile shape.
    Computes it on-demand, if not already cached""";
    g = self._gpgq.get((pq,i,j));
    if g is None:
      g = self._gpgq[pq,i,j] = self.tile_subshape( self.gain.get((pq[0],i),self._unity)*
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
          gp = self.gain.get(pi,self._unity);
          gpi = self._gp_inv[pi] =  1/(gp+regularize);
          gpi[~numpy.isfinite(gpi)] = 0; 
          if pi[0] in verbose_stations_corr:
            print("G:%s:%s"%pi,self.gain[pi][verbose_element]);
            print("Ginv:%s:%s"%pi,gpi[verbose_element]);
        gpgqi.append(gpi);
      g = self._gpgq_inv[pq,i,j] = self.tile_subshape(gpgqi[0]*numpy.conj(gpgqi[1]));
    return g;

  def get_last_timeslot (self):
    """Returns dict of p->g, where p is a parm ID, and g is the gain solution for the last timeslot.
    This can be used to initialize new gain objects (for init_value)"""
    return dict([(key,value[-1,...]) for key,value in self.gain.items() ]);

  def reset_residuals (self):
    self._residual = {};

  def residual (self,lhs,rhs,pq,tiler=None,cache=False):
    """Returns residual R = apply(lhs) - rhs, tiled into subtile shape.
    Computes it on-demand, if not already cached""";
    cache = cache and tiler;
    res = self._residual_cache.get(pq);
    if not cache or res is None:
      corr = self.apply(lhs,pq,index=True,cache=cache,tiler=tiler);
      res = matrix_sub(corr,rhs[pq]);
      if cache:
        self._residual_cache[pq] = res;
      if pq in verbose_baselines_corr:
        print(pq,"D",[ 0 if is_null(g) else g[verbose_element] for g in lhs[pq] ]);
        print(pq,"M",[ 0 if is_null(g) else g[verbose_element] for g in rhs[pq] ]);
        print(pq,"C",[ 0 if is_null(g) else g[verbose_element] for g in corr ]);
        print(pq,"R",[ 0 if is_null(g) else g[verbose_element] for g in res ]);
    return res;

  def residual_inverse (self,lhs,rhs,pq,regularize=0,tiler=None,cache=False):
    """Returns residual R = apply_inverse(lhs) - rhs, tiled into subtile shape.
    Computes it on-demand, if not already cached""";
    cache = False and cache and tiler;
    res = self._residual_inverse_cache.get(pq);
    if not cache or res is None:
      corr = self.apply_inverse(lhs,pq,cache=cache,regularize=regularize,tiler=tiler);
      res = matrix_sub(corr,rhs[pq]);
      if cache:
        self._residual_inverse_cache[pq] = res;
      if pq in verbose_baselines_corr:
        print(pq,"D",[ 0 if is_null(g) else g[verbose_element] for g in lhs[pq] ]);
        print(pq,"M",[ 0 if is_null(g) else g[verbose_element] for g in rhs[pq] ]);
        print(pq,"C",[ 0 if is_null(g) else g[verbose_element] for g in corr ]);
        print(pq,"R",[ 0 if is_null(g) else g[verbose_element] for g in res ]);
    return res;

  def apply (self,lhs,pq,index=True,cache=False,tiler=None):
    """Returns lhs with gains applied: Gp*lhs*Gq^H."""
    cache = False and cache and tiler;
    tiler = tiler or self;
    if index:
      lhs = lhs[pq];
    appl = self._apply_cache.get(pq) if cache else None;
    if appl is None:
#      print [ (getattr(tiler.untile_data(tiler.tile_data(d)),'shape',()),getattr(self.gpgq(pq,i,j),'shape',())) for d,(i,j) in zip(lhs,IJ2x2) ];
      appl = [ 0 if is_null(d) else tiler.untile_data(tiler.tile_data(d)*self.gpgq(pq,i,j))
                                   for d,(i,j) in zip(lhs,IJ2x2) ];
      if cache:
        self._apply_cache[pq] = appl;
      if pq in verbose_baselines_corr:
        print(pq,"LHS",[ 0 if is_null(g) else g[verbose_element] for g in lhs ]);
        print(pq,"APPL(LHS)",[ 0 if is_null(g) else g[verbose_element] for g in appl ]);
    return appl;

  def apply_inverse (self,rhs,pq,cache=False,regularize=0,tiler=None):
    """Returns rhs with inverse gains applied: Gp^{-1}*rhs*Gq^{H-1}."""
    cache = False and cache and tiler;
    tiler = tiler or self;
    appl = self._apply_inverse_cache.get(pq) if cache else None;
    if appl is None:
      mod = rhs[pq];
      appl = [ 0 if is_null(m) else tiler.untile_data(tiler.tile_data(m)*
                  self.gpgq_inv(pq,i,j,regularize=regularize)) for m,(i,j) in zip(mod,IJ2x2) ];
      if cache:
        self._apply_inverse_cache[pq] = appl;
    return appl;

  def get_2x2_gains (self,datashape,expanded_dataslice,tiler=None):
    tiler = tiler or self;
    gainsdict = {};
    for p in self._antennas:
      gx = self.gain.get((p,0),self._unity);
      gy = self.gain.get((p,1),self._unity);
      gainsdict[p] = [ tiler.expand_subshape(gx,datashape,expanded_dataslice),meq.sca_vells(0),meq.sca_vells(0),
                       tiler.expand_subshape(gy,datashape,expanded_dataslice) ];
    return gainsdict;


