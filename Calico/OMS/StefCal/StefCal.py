# -*- coding: utf-8 -*-

from Timba import pynode
from Timba.Meq import meq
import numpy
import sys
import math
import operator
import Kittens.utils
import time
import cPickle
import os.path
import traceback

from MatrixOps import *

_verbosity = Kittens.utils.verbosity(name="stefcal");
dprint = _verbosity.dprint;
dprintf = _verbosity.dprintf;

def GCD (a,b):
  """Return greatest common divisor using Euclid's Algorithm."""
  while b:
    a,b = b,a%b
  return a;

def LCM (a,b,*args):
  """Return lowest common multiple of two arguments."""
  if not args:
    return a*b//GCD(a,b);
  else:
    return reduce(LCM,[a,b]+list(args));

def print_variance (variance):
  """Given a dictionary of per-baseline variances, computes per-station and mean overall variance""";
  meanvar = [];
  meanvar_sta = {};
  for pq in variance.keys():
    dprint(1,"variance on %s-%s is"%pq,[ v if numpy.isscalar(v) else v.flat[0] for v in variance[pq] ]);
    v1 = [ v if numpy.isscalar(v) else v.flat[0] for v in variance[pq] if v != 0 ];
    meanvar += v1;
    meanvar_sta[pq[0]] = meanvar_sta.get(pq[0],[]) + v1;
    meanvar_sta[pq[1]] = meanvar_sta.get(pq[1],[]) + v1;
  # compute mean
  meanvar = math.sqrt((numpy.array(meanvar)**2).mean());
  meanvar_sta = dict([ (p,math.sqrt((numpy.array(x)**2).mean())) for p,x in meanvar_sta.iteritems() ]);
  for p,x in meanvar_sta.iteritems():
    dprint(1,"variance on %s is %f"%(p,x));
  dprint(1,"overall mean variance is",meanvar);


def dump_data_model (data,model,ifrs,filename="dump.txt"):
  ff = file(filename,"w");
  print "dumping data/model to",filename;
  numpy.set_printoptions(threshold=1000000000);
  for pq in ifrs:
    for i,(d,m) in enumerate(zip(data[pq],model[pq])):
      xy = ("xx","xy","yx","yy")[i];
      if not numpy.isscalar(d) and not numpy.isscalar(m):
        ff.write("# data %s-%s %d (%s)\n"%(pq[0],pq[1],i,xy));
        ff.write(numpy.array_str(d)+"\n");
        ff.write("# model %s-%s %d (%s)\n"%(pq[0],pq[1],i,xy));
        ff.write(numpy.array_str(m)+"\n");
  ff.close();


global_gains = {};

class StefCalVisualizer (pynode.PyNode):
  def __init__ (self,*args):
    pynode.PyNode.__init__(self,*args);

  def update_state (self,mystate):
    mystate('freq_average',False);
    mystate('label','G');
    mystate('index',[]);
    self.set_symdeps("Domain");

  def get_result (self,request,*children):
    vellsets = [];
    gainsets = global_gains.get(self.label);
    if not gainsets:
      return meq.result();
    nsets = len(gainsets);
    keys = sorted(gainsets[0].keys());
    self.set_state('plot_label',keys);
    for pp in keys:
      for gains in (gainsets if self.index == [] else [gainsets[self.index]]):
        if self.freq_average:
          vellsets += [ meq.vellset(array_to_vells(x.mean(1)) if x.ndim>1 else x) for x in gains[pp] ];
        else:
          vellsets += [ meq.vellset(x) for x in gains[pp] ];
    res = meq.result(cells=request.cells);
    res.vellsets = vellsets;
    if self.index == []:
      res.dims = [ len(keys),nsets,2,2 ];
    else:
      res.dims = [ len(keys),2,2 ];
    return res;


class StefCalNode (pynode.PyNode):
  def __init__ (self,*args):
    pynode.PyNode.__init__(self,*args);
    self._dataset_id = None;
    self.ifr_gain = {};

  def update_state (self,mystate):
    """Standard function to update our state""";
    # list of ifr pairs (as p:q strings) corresponding to first axis of tensor input
    mystate('ifrs',[]);
    # list of ifr pairs which we use to solve for gains. If empty, all ifrs are used
    mystate('solve_ifrs',[]);
    # correlation names
    mystate('corr_names',["x","y"]);
    # labels for gain, ifr gain and differential gain parameters
    mystate('gain_parm_label',"G");
    mystate('ifr_gain_parm_label',"IG");
    mystate('diffgain_parm_label',"dE");
    # which implementation to use
    # given X.Y.Z, import module X.Y, and use symbol Z from that
    # given just X, use X.X
    mystate('implementation',"GainDiag");
    path = self.implementation.split('.');
    modname,classname = '.'.join(['Calico.OMS.StefCal']+(path[:-1] or path[-1:])),path[-1];
    __import__(modname);
    module = sys.modules[modname];
    self._impl_class = getattr(module,classname);
    self._polarized = getattr(self._impl_class,'polarized',True);
    # convergence criteria
    mystate('epsilon',1e-5);            # updates <epsilon are considered converged
    mystate('diffgain_epsilon',1e-5);            # updates <epsilon are considered converged
    mystate('max_iter',50);             # max gain iter in first major cycle
    mystate('max_iter_1',20);            # max gain iter in middle of major cycle
    mystate('max_iter_2',20);            # max gain iter in last major cycle
    mystate('diffgain_max_iter',5);     # max diffgain iters
    mystate('max_major',10);
    mystate('convergence_quota',0.9);   # what percentage of parms should converge
    mystate('diffgain_convergence_quota',1);   # what percentage of parms should converge
    # subtiling for gains
    mystate('gain_subtiling',[1,1]);
    mystate('diffgain_subtiling',[]);
    # smoothing for gains and differential gains
    mystate('gain_smoothing',[]);
    mystate('diffgain_smoothing',[]);
    # use stored solution (if available) as starting guess
    mystate('init_from_table',True);
    # use previous tile (timeslot) as starting guess -- if table not available
    mystate('init_from_previous',True);
    # use this value as starting guess -- if previous two not available
    mystate('init_value',1);
    # regularization factor applied to gain solutions for correction
    mystate('regularization_factor',0);
    # regularize intermediate corrections (when solving for dEs)?
    mystate('regularize_intermediate',False);
    mystate('iter_twosided',False);
    mystate('iter_averaging',2);
    # return residuals (else data)
    mystate('residuals',True);
    # return corrected residuals/data (else uncorrected)
    mystate('correct',True);
    # flagmask to aassign to flagged visiblities
    mystate('flagmask',1);
    # If True, equation is GDG^H=M. Else use D=GMG^H.
    mystate('gain_bounds',[]);
    mystate('gains_on_data',True);
    # solve for ifr gains as we go along
    mystate('solve_ifr_gains',True);
    # apply previous ifr gain solution, if available
    mystate('apply_ifr_gains',True);
    # name of ifr gain tables
    mystate('ifr_gain_table','ifrgains.cp');
    # visualize G gains by saving them in the global_gains dict. If >1, then all intermediate iterations will also be saved.
    mystate('visualize_gains',0);
    # visualize dE gains by saving them in the gloabl_gains dict
    mystate('visualize_diffgains',0);
    # verbosity level
    mystate('verbose',0);
    # enables dumping of intermediates to text file, if >=0
    mystate('dump_diffgain',[]);
    mystate('dump_domain',[]);
    # print the per-baseline variance of incoming data
    mystate('print_variance',False);
    # lis of all ifrs, as p,q pairs
    self._ifrs = [ tuple(x.split(':')) for x in self.ifrs ];
    # parse set of solvable ifrs
    self._solvable_ifrs = set([ tuple(x.split(":")) for x in (self.solve_ifrs or self.ifrs) ]);
    # other init
    _verbosity.set_verbose(self.verbose);
    _verbosity.enable_timestamps(True,modulo=6000);
    # initial value from which to start iterating
    self._init_value_gain = self.init_value;
    self._init_value_dg = {};

  def get_result (self,request,*children):
    timestamp0 = time.time();
    # get dataset ID from request
    dataset_id,domain_id = meq.split_request_id(request.request_id);
    # get domain ID from request
    time0,time1,timestep,numtime,freq0,freq1,freqstep,numfreq = request.cells.domain.domain_id;

    # if new dataset ID, do setup for start of new dataset
    if dataset_id != self._dataset_id:
      self._dataset_id = dataset_id;
      self._init_value_gain = self.init_value;
      self._init_value_dg = {};
      dprint(1,"new dataset id",dataset_id);
      # if asked to solve for IFR gains, set up dicts for collecting stats
      if self.solve_ifr_gains:
        self.ig_sum_reim = dict([ (pq,[0j]*4) for pq in self._ifrs ]);
        self.ig_sum_sq   = dict([ (pq,[0.]*4) for pq in self._ifrs ]);
        self.ifr_gain_update = dict([ (pq,[0.]*4) for pq in self._ifrs ]);
      # read previous IFR gains from table, if asked to apply them
      self.ifr_gain = {};
      if self.apply_ifr_gains and os.path.exists(self.ifr_gain_table):
        try:
          self.ifr_gain = cPickle.load(file(self.ifr_gain_table));
          dprint(1,"loaded %d ifr gains from %s"%(len(self.ifr_gain),self.ifr_gain_table));
        except:
          traceback.print_exc();
          dprint(1,"error loading gains from",self.ifr_gain_table);

    # child 0 is data
    # child 1 is direction-independent model
    # children 2 and on are models subject to dE terms
    num_diffgains = len(children)-2;
    if num_diffgains < 0:
      raise TypeError,"StefCalNode: at least 2 children (data, model) must be provided";

    # check inputs and populate mappings
    pqij_all = [];      # list of all (p,q),i,j tuples
    pqij_data = [];     # subset of (p,q),i,j tuples for which we have non-null input
    pqij_solvable = []; # subset of (p,q),i,j tuples for which we solve for gains
    data  = {};         # mapping from (p,q) to four data time-freq planes
    model0 = {};        # mapping from (p,q) to four model (M0) time-freq planes
    dgmodel = [ {} for i in range(num_diffgains) ];
                        # for each diff gain, mapping from (p,q) to M1,M2,... model time-freq planes (4 each)
    model = {};         # this is the full model, M0+M1+M2+...
    flags = {};         # this is a bool array of flagged values

    antennas = set([p for p,q in self._ifrs]) | set([q for p,q in self._ifrs]);

    # per-baseline noise
    variance = {};
    # antennas for which we have non-trivial data
    solvable_antennas = set();

    datares = children[0]
    modelres = children[1];
    if any( [ ch.dims != datares.dims for ch in children[1:] ] ):
      raise TypeError,"tensor dimensions of data and model(s) must match";
    # expecting Nx2x2 matrices
    if len(datares.dims) == 3:
      if datares.dims[1] != 2 or datares.dims[2] != 2:
        raise TypeError,"data and model must be of rank Nx2x2";
      nifrs = datares.dims[0];
      # setup antenna names
      if nifrs != len(self.ifrs):
        raise TypeError,"first dimension of data and model must match the number of interferometers in the ifrs field";
      # setup list of data, values and parameter names
      nvells = 0;
      for pq in self._ifrs:
        # get IFR gain for this p,q
        ifrgain = self.ifr_gain.get(pq,[1,1,1,1]);
        flag4 = False;
        # now loop over the 4 matrix elements
        for num,(i,j) in enumerate(IJ2x2):
          # get data, and apply ifr gains if we have them
          d = getattr(datares.vellsets[nvells],'value',0);
          if not is_null(d):
            d *= ifrgain[num];
          # get model
          m = getattr(modelres.vellsets[nvells],'value',0);
          if hasattr(datares.vellsets[nvells],'flags'):
            flag4 += (datares.vellsets[nvells].flags != 0);
          # if model and/or data is null, then we're unpolarized, so skip this from the matrix entirely
          skip = is_null(d) if self._polarized else (is_null(m) or is_null(d));
          if skip:
            pass;
          else:
            # if this is the first datum, then check shape, and prepare subtilings etc.
            if pq in self._solvable_ifrs:
              solvable_antennas.update(pq);
            if not nvells:
              # this is the basic time-frequency shape
              datashape = tuple(d.shape);
              # figure out subtiling
              # if not specified, use whole tile as solution interval
              if self.gain_subtiling:
                # replace nulls in subtiling with solution interval
                lcm_subtiling = gain_subtiling = [ min(gs,ds) or ds for gs,ds in zip(self.gain_subtiling,datashape) ];
              else:
                lcm_subtiling = gain_subtiling = datashape;
              if len(gain_subtiling) != len(datashape):
                raise ValueError,"gain_subtiling vector must have the same length as the data shape";
              if min(gain_subtiling) < 1:
                raise ValueError,"invalid gain_subtiling %s"%self.gain_subtiling;
              # if diffgains are also present, then work out the least-common-multiple subtiling
              if num_diffgains:
                dg_subtiling = self.diffgain_subtiling or datashape;
                dg_subtiling = [ min(gs,ds) or ds for gs,ds in zip(dg_subtiling,datashape) ];
                if len(dg_subtiling) != len(datashape):
                  raise ValueError,"diffgain_subtiling vector must have the same length as the data shape";
                if min(dg_subtiling) < 1:
                  raise ValueError,"invalid diffgain_subtiling %s"%dg_subtiling;
                lcm_subtiling = [ LCM(a,b) for a,b in zip(gain_subtiling,dg_subtiling) ];
              # data must be expanded to match the LCM subtiling
              expanded_datashape = tuple([ (nd/np+(1 if nd%np else 0))*np for nd,np in zip(datashape,lcm_subtiling) ]);
              dprint(1,"gain parm LCM subtiling is",lcm_subtiling);
              # this counts how many valid visibilities we have per each antenna, per each time/freq slot
              vis_per_antenna = dict([(p,numpy.zeros(expanded_datashape,dtype=int)) for p in antennas ]);
              # if tiling does not tile the data shape perfectly, we'll need to expand the input arrays
              # Define pad_array() as a function for this: it will be identity if no expansion is needed
              if datashape != expanded_datashape:
                expanded_dataslice = tuple([ slice(0,nd) for nd in datashape ]);
                def pad_array (x):
                  if is_null(x):
                    return 0;
                  x1 = numpy.zeros(expanded_datashape,dtype=x.dtype);
                  x1[expanded_dataslice] = x;
                  return x1;
                self._expanded_size = reduce(operator.mul,expanded_datashape);
                self._expansion_ratio = self._expanded_size/float(reduce(operator.mul,datashape));
                dprint(1,"input arrays will be expanded to shape",expanded_datashape,"ratio %.2f"%self._expansion_ratio);
              else:
                self._expanded_size = reduce(operator.mul,datashape);
                self._expansion_ratio = 1;
                expanded_dataslice = None;
                pad_array = identity_function;
            # now check inputs and add them to data and model dicts
            if d.shape != datashape:
              print d.shape,datashape,d;
              raise TypeError,"data shape mismatch at %s:%s:%s:%s"%(pq[0],pq[1],self.corr_names[i],self.corr_names[j]);
            if not is_null(m) and m.shape != datashape:
              raise TypeError,"model shape mismatch at %s:%s:%s:%s"%(pq[0],pq[1],self.corr_names[i],self.corr_names[j]);
            # add to data/model matrices, applying the padding function defined above
            m0 = model0.setdefault(pq,[0,0,0,0])[num] = pad_array(m);
            data.setdefault(pq,[0,0,0,0])[num] = pad_array(d);
            # add the noise variance
            if self.print_variance:
              v = d[1:,...]-d[:-1,...];
              v = v[numpy.isfinite(v)];
              variance.setdefault(pq,[0,0,0,0])[num] = (v.real.std(0)+v.imag.std(0))/2;
            # also accumulate initial model, as M0+M1+M2
            # if max_major==0, then we don't solve for diff
            if num_diffgains:
              m0 = model.setdefault(pq,[0,0,0,0])[num] = 0 if is_null(m0) else m0.copy();
              for k in range(num_diffgains):
                m1 = children[2+k].vellsets[nvells].value
                m1 = pad_array(m1);
                dgmodel[k].setdefault(pq,[0,0,0,0])[num] = m1;
                m0 += m1;
            else:
              model[pq] = model0[pq];
          # increment vells #
          nvells += 1;
        # do we have any flags? apply them to all 4 corrs
        if flag4 is not False:
          dprint(4,pq,"has",flag4.sum(),"flags");
          for dataset in [data,model,model0] + dgmodel:
            for x in dataset[pq]:
              if not is_null(x):
                x[pad_array(flag4)] = 0;
        # now count valid elements
        # nel will be a time/freq array with values from 0 to 4, depending on how many data/model visiblities are non-null
        nel = sum([0 if is_null(d) else ((d!=0)|(m!=0)).astype(int) for d,m in zip(data[pq],model[pq])]);
        valid = (nel>1).astype(int)
        vis_per_antenna[pq[0]] += valid;
        vis_per_antenna[pq[1]] += valid;
    else:
      # in principle could also handle [N], but let's not bother for now
      raise TypeError,"data and model must be of rank Nx2x2";

    GainClass = self._impl_class;
    dprintf(0,"solving with %s, using %d of %d inteferometers, with %d solvable antennas\n",
      GainClass.__name__,len(self._solvable_ifrs),len(self.ifrs),len(solvable_antennas));

    # flags will contain a flag array per each baseline (but not per correlation), of shape expanded_datashape.
    # These are additional flags raised in the processing below. Note that data and model will be set to 0 at each flagged point.
    flags = {};
    # print [ (p,vis_per_antenna[p][0:10]) for p in antennas ];
    ## this tells us how many elements there are in correlation matrix pq (per time/freq point)
    ## sum([0 if is_null(d) else ((d!=0)|(m!=0)).astype(int) for d,m in zip(data[pq],model[pq])])
    # this gives us the number of valid visibilities (in the solvable IFR set) per each time/freq slot
    nel_pq = dict([(pq,sum([0 if is_null(d) else ((d!=0)|(m!=0)).astype(int)
                    for d,m in zip(data[pq],model[pq])])) for pq in list(self._solvable_ifrs) ]);
    # this is how many valid visibility matrices we have total per slot
    num_valid_vis = sum([ int(n>1) if numpy.isscalar(n) else (n>1).astype(int) for pq,n in nel_pq.iteritems() ]);
    # we need as many as there are parameters per antenna, times ~4, to constrain the problem fully
    threshold = len(solvable_antennas)*6;
    dprintf(1,"%d valid visibility matrices per time/freq slot are required for a solution\n",threshold);
    invalid_slots = (num_valid_vis>0)&(num_valid_vis<threshold);
    nfl = invalid_slots.sum();
    if nfl:
      dprintf(1,"%d of %d slots have insufficient visibilities and will be flagged\n",nfl,self._expanded_size);
      flags = dict([(pq,invalid_slots.copy()) for pq in self._ifrs ]);
      for dataset in [data,model,model0] + dgmodel:
        for dd in dataset.itervalues():
          for d in dd:
            if not is_null(d):
              d[invalid_slots] = 0;
      # init flags array
      flags = dict([(pq,invalid_slots.copy()) for pq in self._ifrs]);
#    mnorm = sum([sum([d*numpy.conj(d) for d in dd]) for dd in model.itervalues() ]);

    if self.gain_smoothing:
      dprint(0,"a Gaussian smoothing kernel of size",self.gain_smoothing,"will be applied");
    if self.gain_bounds:
      dprint(0,"gains will be flagged on amplitudes outside of",self.gain_bounds);
    dprint(1,"initial gain value is",self._init_value_gain.values()[0].flat[0] if isinstance(self._init_value_gain,dict) else
      self._init_value_gain);
    # init gain parms object
    gain = GainClass(datashape,expanded_datashape,gain_subtiling,self._solvable_ifrs,
              self.epsilon,self.convergence_quota,smoothing=self.gain_smoothing,
              bounds=self.gain_bounds,
              twosided=self.iter_twosided,averaging=self.iter_averaging,
              init_value=self._init_value_gain,
              verbose=self.verbose);

    if self.print_variance:
      print_variance(variance);

    # init diffgains
    if num_diffgains:
      dprintf(0,"also solving for %d differential gains\n",num_diffgains);
      if self.diffgain_smoothing:
        dprint(0,"a Gaussian smoothing kernel of size",self.diffgain_smoothing,"will be applied");
      for i in range(num_diffgains):
        initval = self._init_value_dg.get(i,self.init_value);
        dprint(1,"initial gain value #%d is"%i,initval.items()[0] if isinstance(initval,dict) else initval);
      diffgains = [
        GainClass(datashape,expanded_datashape,dg_subtiling,self._solvable_ifrs,
          self.diffgain_epsilon,self.diffgain_convergence_quota,
          smoothing=self.diffgain_smoothing,
          twosided=self.iter_twosided,averaging=self.iter_averaging,
          init_value=self._init_value_dg.get(i,self.init_value),
          verbose=self.verbose)
        for i in range(num_diffgains) ];
    else:
      diffgains = [];

    gain_iterate = (data,model) if self.gains_on_data else (model,data);
    # start major loop -- alternates over gains and diffgains
    for nmajor in range(self.max_major+1):
      if (domain_id == self.dump_domain or self.dump_domain == -1):
        dump_data_model(data,model,self._solvable_ifrs,"dump_G%d.txt"%nmajor);
      # first iterate normal gains to convergence
      gain_maxdiffs = [];
      maxiter = (self.max_iter_2 if nmajor == self.max_major else self.max_iter_1) if nmajor else self.max_iter;
#      print maxiter,self.max_iter,self.max_iter_1,self.max_iter_2;
      t0 = time.time();
      chisq0 = 0;
      for niter in range(maxiter):
        # iterate over normal gains
        converged,maxdiff,deltas,nfl = gain.iterate(niter=niter,*gain_iterate);
        dprint(3,"%d slots converged, max delta is %f"%(gain.num_converged,gain.delta_max));
        gain_maxdiffs.append(maxdiff);
        if self.visualize_gains > 1:
          global_gains.setdefault('G',[]).append(gain.get_2x2_gains(datashape,expanded_dataslice));
        if nfl:
          dprint(2,"%d out-of-bound gains flagged at iteration %d"%(nfl,niter));
#        print "value",gain.gain.values()[0][0][0,0];
        # check chi-square
        if 1: # ( niter and not niter%100 ) or niter >= maxiter-1 or converged:
          chisq = self.compute_chisq(gain,data,model,self.gains_on_data);
          dprint(2,"iter %d max gain update is %g chisq is %.12g diff %.3g"%(niter+1,gain.delta_max,chisq,(chisq0-chisq)/chisq));
          chisq0 = chisq;
        # break out if converged
        if converged:
          break;
      dprint(1,"gains %sconverged, chisq %.12g (last G update %g) after %d iterations and %.2fs"%(
            ("" if converged else "not "),chisq,gain.delta_max,niter+1,time.time()-t0));
      dprint(2,"  convergence was"," ".join(["%.2g"%x for x in gain_maxdiffs]));
      # break out if no diffgains to iterate over, or if we're on the last major cycle
      if not num_diffgains or nmajor >= self.max_major:
        break;
      else:
        # we have solved for G=inverse of G-Jones essentially, thus minimzing
        # G*D*G^H <- M0+corrupt(M1)+corrupt(M2)+...
        # model is the full model, M0+corrupt(M1)+corrupt(M2)+....
        # subtract this from corrected data: D1 = G*D*G^H - M0 - corrupt(M1) - corrupt(M2) - ... to obtain residuals
        if self.gains_on_data:
          data1 = dict([ (pq,gain.residual(data,model,pq)) for pq in self._solvable_ifrs ]);
        else:
          data1 = dict([ (pq,gain.residual_inverse(data,model,pq,
              regularize=self.regularization_factor if self.regularize_intermediate else 0))
              for pq in self._solvable_ifrs ]);
#        for pq in list(self._solvable_ifrs)[:1]:
#          print [ (pq,[ d1 if is_null(d1) else abs(d1).max() for d1 in data1[pq] ]) for pq in self._solvable_ifrs ];
        # now loop over all diffgains and iterate each set once
        for idg,dg in enumerate(diffgains):
          # add current estimate of corrupt(Mi) back into data1, and subtract from current model
          for pq in self._solvable_ifrs:
            corr = dg.apply(dgmodel[idg],pq,cache=True);
            mm = model[pq];
            for n,(d,m,c) in enumerate(zip(data1[pq],mm,corr)):
              d += c;
              mm[n] = -c if is_null(m) else m-c;
#          print data1['0','A'][0][0,0],model0['0','A'][0][0,0],model['0','A'][0][0,0],dgmodel[idg]['0','A'][0][0,0];
#          print "d1",data1['0','A'][0],"m1",dgmodel[idg]['0','A'][0];
          if ( idg == self.dump_diffgain or self.dump_diffgain == -1 ) and \
             ( domain_id == self.dump_domain or self.dump_domain == -1 ):
            dump_data_model(dgmodel[idg],data1,self._solvable_ifrs,"dump_E%d-%d.txt"%(idg,nmajor));
          # iterate this diffgain solution
          t0 = time.time();
          chisq0 = -1;
          gain_maxdiffs = [];
          for niter in range(self.diffgain_max_iter):
            converged,maxdiff,deltas,nfl = dg.iterate(dgmodel[idg],data1,niter=niter);
            dprint(3,"%d slots converged, max delta is %f"%(dg.num_converged,dg.delta_max));
            gain_maxdiffs.append(maxdiff);
            if 1:# (niter and not niter%100) or niter >= self.diffgain_max_iter-1 or converged:
              chisq = self.compute_chisq(dg,data1,dgmodel[idg],False);
              dprint(2,"dg%d iter %d max gain update is %g chisq is %.12g diff %.3g"%(idg,niter+1,gain.delta_max,chisq,(chisq0-chisq)/chisq));
              chisq0 = chisq;
            if converged:
              break;
          dprint(2,"diffgain #%d converged after %d iterations and %.2fs"%(idg,niter,time.time()-t0));
          dprint(2,"  convergence was"," ".join(["%.2g"%x for x in gain_maxdiffs]));
#          print dg.gain.keys();
#          print "dE(0)",dg.gain['0',0];
#          print "dgcorrupt",dg.corrupt(dgmodel[idg],('0','A'),cache=True)[0];
          # add back to model, and subtract from data1 if needed
          for pq in self._solvable_ifrs:
            corr = dg.apply(dgmodel[idg],pq,cache=True);
            for d,m,c in zip(data1[pq],model[pq],corr):
              m += c;
              if idg<num_diffgains-1:
                d -= c;

    # if we were solving for diffgains, then model is not completely up-to-date, since the non-solvable baselines
    # have been ignored. Fill them in here. Also, reset residuals
    if num_diffgains:
      for pq in set(self._ifrs)-set(self._solvable_ifrs):
        model[pq] = model0[pq];
        for i,dg in enumerate(diffgains):
          for m,c in zip(model[pq],dg.apply(dgmodel[i],pq,cache=True)):
            m += c;

    # visualize gains
    if self.visualize_gains == 1:
      global_gains['G'] = [ gain.get_2x2_gains(datashape,expanded_dataslice) ];
    if self.visualize_diffgains and num_diffgains:
      for i,dg in enumerate(diffgains):
        global_gains['dE:%d'%i] = [ dg.get_2x2_gains(datashape,expanded_dataslice) ];

    # remember init value for next tile
    if self.init_from_previous:
      self._init_value_gain = gain.get_last_timeslot();
      for i,dg in enumerate(diffgains):
        self._init_value_dg[i] = dg.get_last_timeslot();

    # update IFR gain solutions, if asked to
    if self.solve_ifr_gains:
      for pq in self._ifrs:
        dd = data[pq];
        if self.gains_on_data:
          mm = gain.apply_inverse(model,pq,cache=True,regularize=self.regularization_factor);
        else:
          mm = gain.apply(model,pq,cache=True);
        for num,(d,m) in enumerate(zip(dd,mm)):
          # work out update to ifr gains
          if numpy.isscalar(m):
            m = numpy.array(m);
          if numpy.isscalar(d):
            d = numpy.array(d);
          dh = numpy.conj(d);
          sri = self.ig_sum_reim[pq][num] = self.ig_sum_reim[pq][num] + (m*dh).sum();
          ssq = self.ig_sum_sq[pq][num]   = self.ig_sum_sq[pq][num] + (d*dh).sum();
          if ssq != 0:
            self.ifr_gain_update[pq][num] = sri/ssq;
#          if num == 0 and pq[0] == '0':
#           print m[0,0],d[0,0],dh[0,0]
#            print pq,(m*dh).sum(),(d*dh).sum(),sri/ssq;

    # work out result -- residual or corrected visibilities, depending on our state
    variance = {};
    nvells = maxres = 0;
    for pq in self._ifrs:
      m = model.get(pq);
      if self.gains_on_data:
        out = res = gain.residual(data,model,pq);
        if not self.residuals:
          out = gain.apply(data,pq);
      else:
        out = res = gain.residual_inverse(data,model,pq,regularize=self.regularization_factor);
        if not self.residuals:
          out = gain.apply_inverse(data,pq,regularize=self.regularization_factor);
      fl4 = flags.get(pq);
      for n,x in enumerate(out):
        vs = datares.vellsets[nvells];
        val = getattr(vs,'value',None);
        if val is not None:
          try:
            val[...] = x[expanded_dataslice] if expanded_dataslice \
              and not is_null(x) else x;
          except:
            print x,getattr(x,'shape',None);
          if self.print_variance:
            v = val[1:,...] - val[:-1,...];
            v = v[numpy.isfinite(v)];
            variance.setdefault(pq,[0,0,0,0])[n] = (v.real.std(0)+v.imag.std(0))/2;
        if fl4 is not None:
          fl = getattr(vs,'flags',None);
          if fl is None:
            fl = vs.flags = meq.flags(expanded_datashape);
          fl[fl4] |=  self.flagmask;
        # compute stats
        maxres = max(maxres,abs(res[n]).max() if not numpy.isscalar(res[n]) else abs(res[n]));
        nvells += 1;

    if self.print_variance:
      print_variance(variance);

    # if last domain, then write ifr gains to file
    if self.solve_ifr_gains and time1 >= numtime:
      # apply updates
      for pq in self._ifrs:
        self.ifr_gain[pq] = [ g*g1 for g,g1 in zip(self.ifr_gain.get(pq,[1,1,1,1]),self.ifr_gain_update[pq]) ];
      dprint(2,"IFR gain solutions update: ",", ".join(
            ["%s%s:%s%s %s"%(p,self.corr_names[i],q,self.corr_names[j],
            self.ifr_gain_update[(p,q),i,j])
            for (p,q),i,j in pqij_data[0:3]]));
      dprint(2,"IFR gain solutions: ",", ".join(
            ["%s%s:%s%s %s"%(p,self.corr_names[i],q,self.corr_names[j],
            self.ifr_gain[(p,q),i,j])
            for (p,q),i,j in pqij_data[0:3]]));
      # save
      try:
        cPickle.dump(self.ifr_gain,file(self.ifr_gain_table,'w'));
        dprint(1,"saved %d ifr gains to %s"%(len(self.ifr_gain_update),self.ifr_gain_table));
      except:
        traceback.print_exc();
        dprint(0,"error saving ifr gains to",self.ifr_gain_table);

    dt = time.time()-timestamp0;
    m,s = divmod(dt,60);
    dprint(0,"%s residual max %g last chisq %g (last G update %g), elapsed time %dm%0.2fs"%(
              request.request_id,maxres,chisq,gain.delta_max,m,s));

    return datares;

  def compute_chisq (self,gain,data,model,gains_on_data):
    chisq = 0;
    nterms = 0;
    for pq in self._solvable_ifrs:
      if gains_on_data:
        res = gain.residual(data,model,pq);
      else:
        res = gain.residual(model,data,pq);
      for r in res:
        if numpy.isscalar(r):
          chisq1 = (r*numpy.conj(r));
          nterms += 1; # self._expanded.size;
        else:
          chisq1 = (r*numpy.conj(r)).sum();
          nterms += r.size;
#          print pq,"chi-sqare contribution is",chisq1,"r",r if numpy.isscalar(r) else r[5,0];
        chisq += chisq1;
    # nterms may be exagerrated due to data padding, so correct by the expansion ratio
    return chisq/nterms*self._expansion_ratio;


