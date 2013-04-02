
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
import scipy.ndimage.measurements

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
    mystate('flag_unity',True);
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
    # define processing function
    if self.freq_average:
      if self.flag_unity:
        def gain_to_vellset (ix,x):
          x = numpy.ma.masked_array(x,x==(1 if ix in (0,3) else 0));
          if x.ndim > 1:
            x = x.mean(1);
          return meq.vellset(array_to_vells(x),flags=mask_to_flags(x.mask)) if x.mask.any() else meq.vellset(array_to_vells(x));
      else:
        def gain_to_vellset (ix,x):
          return meq.vellset(array_to_vells(x.mean(1)) if x.ndim>1 else x);
    else:
      if self.flag_unity:
        def gain_to_vellset (ix,x):
          x = numpy.ma.masked_array(x,x==(1 if ix in (0,3) else 0));
          return meq.vellset(array_to_vells(x),flags=mask_to_flags(x.mask)) if x.mask.any() else meq.vellset(array_to_vells(x));
      else:
        def gain_to_vellset (ix,x):
          return meq.vellset(array_to_vells(x));
    
    for pp in keys:
      for gains in (gainsets if self.index == [] else [gainsets[self.index]]):
        xx,xy,yx,yy = gains[pp];
        # xy /= xx
        # yx /= yy
        vellsets += [ gain_to_vellset(ix,x) for ix,x in enumerate((xx,xy,yx,yy)) ];

    res = meq.result(cells=request.cells);
    res.vellsets = vellsets;
    if self.index == []:
      res.dims = [ len(keys),nsets,2,2 ] if nsets>1 else [ len(keys),2,2 ];
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
    # list of baseline lengths
    mystate('baselines',[]);
    if not self.baselines:
      self.baselines = [0]*len(self.ifrs);
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
    mystate('epsilon',1e-5);            # when the update is ||G-G'||<epsilon, we are converged
    mystate('delta',1e-6);              # when chisq changes by less than delta, we are converged. First loop of major cycle
    mystate('delta_1',1e-6);            # middle loops of major cycle
    mystate('delta_2',1e-6);            # last loop of major cycle
    mystate('diffgain_epsilon',1e-5);   # epsilon for diffgain solutions
    mystate('diffgain_delta',1e-6);     # delta for diffgain solutions
    mystate('max_iter',50);             # max gain iter in first loop of major cycle
    mystate('max_iter_1',20);           # max gain iter in middle loops of major cycle
    mystate('max_iter_2',20);           # max gain iter in last loop of major cycle
    mystate('max_diverge',2);           # give up if we diverge for that many iterations in a row
    mystate('diffgain_max_iter',5);     # max diffgain iters
    mystate('diffgain_max_diverge',2);  # give up if we diverge for that many iterations in a row
    mystate('max_major',10);
    mystate('convergence_quota',0.9);   # what percentage of parms should converge
    mystate('diffgain_convergence_quota',1);   # what percentage of parms should converge
    # weigh G solutions using noise estimates
    mystate('weigh_gains',True);
    # weigh dE solutions using noise estimates
    mystate('weigh_diffgains',True);
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
    # rescale data to model
    mystate('rescale',True);
    # use this value as starting guess -- if previous two not available
    mystate('init_value',1);
    # if set, then updated solutions are forward-fed during iteration
    mystate('omega',.5);
    mystate('omega_1',.5);
    mystate('omega_de',.5);
    # avergaing modes
    mystate('average',2);
    mystate('average_1',2);
    mystate('average_de',2);
    mystate('feed_forward',False);
    mystate('feed_forward_1',False);
    mystate('feed_forward_de',False);
    # regularization factor applied to gain solutions for correction
    mystate('regularization_factor',0);
    # regularize intermediate corrections (when solving for dEs)?
    mystate('regularize_intermediate',False);
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
    # IFR gains are per-frequency, or one value across entire band
    mystate('per_chan_ifr_gains',False);
    # IFR gains are for diagonal elements only if True, for full 2x2 if False
    mystate('diag_ifr_gains',False);
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
    # make list of ifrs sorted by baselines
    self.ifr_by_baseline = zip(self._ifrs,self.baselines);
    self.ifr_by_baseline.sort(lambda x,y:cmp(x[1],y[1]));
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
        self.ifr_gain_update = dict([ (pq,[1.]*4) for pq in self._ifrs ]);
      # read previous IFR gains from table, if asked to apply them
      self.ifr_gain = {};
      if self.apply_ifr_gains and os.path.exists(self.ifr_gain_table):
        try:
          self.ifr_gain = cPickle.load(file(self.ifr_gain_table));
          dprint(1,"loaded %d ifr gains from %s"%(len(self.ifr_gain),self.ifr_gain_table));
          # reset off-diagonals to 1
          if self.diag_ifr_gains:
            dprint(1,"resetting off-diagonal elements to 1");
            for gg in self.ifr_gain.itervalues():
              gg[1] = gg[2] = 1;
            print self.ifr_gain[self.ifr_gain.keys()[0]];
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

    # This will contain a mask of flagged values. I would use masked arrays,
    # but they seem to slow something down, so no no. Instead, we'll use 0.0 for missing values in model and
    # data (since zeroes do not upset the equations), and maintain a flagmask for computing things like
    # stats. For each p,q, flagmask contains a boolean array (same shape as model/data).
    flagmask = {};

    antennas = set([p for p,q in self._ifrs]) | set([q for p,q in self._ifrs]);

    # per-baseline noise
    variance = {};
    # antennas for which we have non-trivial data
    solvable_antennas = set();
    # this will count the valid visibilities per each antenna, per each time/freq slot
    vis_per_antenna = None;

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
      nvells = -1;
      for pq in self._ifrs:
        # get IFR gain for this p,q
        ifrgain = self.ifr_gain.get(pq,[1,1,1,1]);
#        print pq,"IFR gain is",ifrgain;
        # now loop over the 4 matrix elements
        for num,(i,j) in enumerate(IJ2x2):
          # increment vells count upfront (this is why we start at -1)
          nvells += 1;
          # get data, and apply ifr gains if we have them
          d = getattr(datares.vellsets[nvells],'value',0);
          if not is_null(d):
            d *= ifrgain[num];
          # get model
          m = getattr(modelres.vellsets[nvells],'value',0);
          if hasattr(datares.vellsets[nvells],'flags'):
            flags = (datares.vellsets[nvells].flags != 0);
          else:
            flags = None;
          # does this need to be skipped? only process data otherwise
          if not ( is_null(d) if self._polarized else (is_null(m) or is_null(d)) ):
            # if this is the first datum, then check shape, and prepare subtilings etc.
            if pq in self._solvable_ifrs:
              solvable_antennas.update(pq);
            # for the first valid result, setup shapes and stuff
            if not model0:
              # this is the basic time-frequency shape
              datashape = tuple(d.shape);
              self._datasize = reduce(operator.mul,datashape);
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
                self._expansion_ratio = self._expanded_size/float(self._datasize);
                dprint(1,"input arrays will be expanded to shape",expanded_datashape,"ratio %.2f"%self._expansion_ratio);
              else:
                self._expanded_size = self._datasize;
                self._expansion_ratio = 1;
                expanded_dataslice = None;
                pad_array = identity_function;
              # this counts how many valid visibilities we have per each antenna, per each time/freq slot
              vis_per_antenna = dict([(p,numpy.zeros(expanded_datashape,dtype=int)) for p in antennas ]);
            # now check inputs and add them to data and model dicts
            if d.shape != datashape:
              print d.shape,datashape,d;
              raise TypeError,"data shape mismatch at %s:%s:%s:%s"%(pq[0],pq[1],self.corr_names[i],self.corr_names[j]);
            if not is_null(m) and m.shape != datashape:
              raise TypeError,"model shape mismatch at %s:%s:%s:%s"%(pq[0],pq[1],self.corr_names[i],self.corr_names[j]);
            # add to data/model matrices, applying the padding function defined above
            m0 = model0.setdefault(pq,[0,0,0,0])[num] = pad_array(m);
            d0 = data.setdefault(pq,[0,0,0,0])[num] = pad_array(d);
            # apply flags
            if flags is not None:
              flags = pad_array(flags);
              if not is_null(m0):
                m0[flags] = 0;
              if not is_null(d0):
                d0[flags] = 0;
              invalid = (d0==0)&(m0==0);
              if pq in flagmask:
                flagmask[pq] |= invalid;
              else:
                flagmask[pq] = invalid;
  #              print pq,validmask[pq];
            # also accumulate initial model, as M0+M1+M2
            # if max_major==0, then we don't solve for diff
            if num_diffgains:
              m0 = model.setdefault(pq,[0,0,0,0])[num] = 0 if is_null(m0) else m0.copy();
              for k in range(num_diffgains):
                m1 = children[2+k].vellsets[nvells].value
                m1 = pad_array(m1);
                if flags is not None and not is_null(m1):
                  m1[flags] = 0;
                dgmodel[k].setdefault(pq,[0,0,0,0])[num] = m1;
                m0 += m1;
        # ok, done looping over the 2x2 visibility matrix elements. 
        # If we have found anything valid at all, finalize flagmasks etc.
        if pq in model0:
          # copy M0 to model, if no diffgains were present
          if not num_diffgains:
            model[pq] = model0[pq];
          # look at flagmask to see how many valid correlations we have, and zero the flagged ones
          fmask = flagmask.get(pq);
          if fmask is not None:
            valid = self._expanded_size - fmask.sum();
            if valid > 0:
              dprint(4,pq,"has %d of %d unflagged correlation matrices"%(valid,self._datasize));
              for dataset in [data,model,model0] + dgmodel:
                for x in dataset.get(pq,[]):
                  if not is_null(x):
                    x[fmask] = 0;
              validmask = (~fmask).astype(int);
              vis_per_antenna[pq[0]] += validmask;
              vis_per_antenna[pq[1]] += validmask;
            else:
              # if nothing is valid, remove baseline from dicts
              dprint(4,pq,"is completely flagged, skipping");
              for dataset in [data,model,model0] + dgmodel:
                del dataset[pq];
          else:
            dprint(4,pq,"has no flagged correlation matrices, all data is valid");
            vis_per_antenna[pq[0]] += 1;
            vis_per_antenna[pq[1]] += 1;
    else:
      # in principle could also handle [N], but let's not bother for now
      raise TypeError,"data and model must be of rank Nx2x2";
    
    valid_ifrs = data.keys();
    solvable_ifrs = set(self._solvable_ifrs)&set(valid_ifrs);
    if not solvable_ifrs:
      dprint(1,"no valid data found for solvable IFRs  -- nothing to stefcal!");
      return datares;

    # compute the noise estimate, and weights based on this
    if self.weigh_gains:
      nw = nm = nnw = 0;
      noise,weight = self.compute_noise(data,flagmask);
      # print and total up stats, check for funny situations
      for pq,bl in self.ifr_by_baseline:
        if pq in data:
          w = weight.get(pq,0);
          rms = noise.get(pq,0);
          if rms is not None:
            dprint(4,"%20s %.2fm"%("-".join(pq),bl),"noise %.2g weight %.2f"%(rms,w));
          if w:
            nw += 1;
          else:
            if is_null(data[pq][0]):
              nm += 1;
            else:
              fm = flagmask.get(pq);
              np = self._datasize if fm is None else (self._expanded_size - fm.sum());
              if np:
                dprint(1,"funny,",pq,"has zero weights, yet %d valid points"%np);
                nnw += 1;
              else:
                nm += 1;
      dprint(1,"%d baselines with non-zero weights (%d with null weight, %d have no valid data)"%(nw,nnw,nm));
    else:
      weight = noise = None;

    GainClass = self._impl_class;
    dprintf(0,"solving with %s, using %d of %d inteferometers (%d have valid data), with %d solvable antennas\n",
      GainClass.__name__,len(self._solvable_ifrs),len(self.ifrs),len(solvable_ifrs),len(solvable_antennas));

    # flags will contain a flag array per each baseline (but not per correlation), of shape expanded_datashape.
    # These are additional flags raised in the processing below. Note that data and model will be set to 0 at each flagged point.
    flags = {};
    ## this tells us how many elements there are in correlation matrix pq (per time/freq point)
    ## sum([0 if is_null(d) else ((d!=0)|(m!=0)).astype(int) for d,m in zip(data[pq],model[pq])])
    # this gives us the number of valid visibilities (in the solvable IFR set) per each time/freq slot
    nel_pq = dict([(pq,sum([0 if is_null(d) else ((d!=0)|(m!=0)).astype(int)
                    for d,m in zip(data[pq],model[pq])])) for pq in data.iterkeys() ]);
    # this is how many valid visibility matrices we have total per slot
    num_valid_vis = sum([ int(n>1) if numpy.isscalar(n) else (n>1).astype(int) for pq,n in nel_pq.iteritems() ]);
    # we need as many as there are parameters per antenna, times ~4, to constrain the problem fully
    threshold = len(solvable_antennas)*4;
    dprintf(1,"%d valid visibility matrices per time/freq slot are required for a solution\n",threshold);
    invalid_slots = (num_valid_vis>0)&(num_valid_vis<threshold);
    nfl = invalid_slots.sum();
    if nfl:
      dprintf(1,"%d of %d slots have insufficient visibilities and will be flagged\n",nfl,self._expanded_size);
      flags = dict([(pq,invalid_slots.copy()) for pq in data.iterkeys() ]);
      for pq,fl in flagmask.iteritems():
        fl[invalid_slots] = True;
      for dataset in [data,model,model0] + dgmodel:
        for dd in dataset.itervalues():
          for d in dd:
            if not is_null(d):
              d[invalid_slots] = 0;


## -------------------- rescale model/data
    if self.rescale:
      scale = {};
      finite = {};
      # compute scales as s(p) = ||sum_q Mpq||/||sum_q Dpq||
      for p in antennas:
        dsum = msum = 0;
        for q in antennas:
          d = m = None;
          if (p,q) in data:
            d,m = data.get((p,q)),model.get((p,q));
          elif (q,p) in data:
            d,m = data.get((q,p)),model.get((q,p));
          if d is None or m is None:
            continue;
          dsum += sum([x*numpy.conj(x) for x in d]);
          msum += sum([x*numpy.conj(x) for x in m]);
        if dsum is not 0:
          s = numpy.sqrt(msum.real/dsum.real);
          f = numpy.isfinite(s);
        if dsum is 0 or (~f).all():
          dprint(2,"no valid data (and thus no scale) for",p);
        else:
          s[~f] = 0;
          scale[p] = s,f;
      # apply scales
      dprint(1,"min/max scaling factors are",min([s[f].min() for s,f in scale.itervalues()]),
                                             max([s.max() for s,f in scale.itervalues()]));
      dprint(2,"per-antenna data scaling factors are "," ".join(["%s %.3g,"%(p,s.max()) for p,(s,f) in scale.iteritems() ]));
      for (p,q),dd in data.iteritems():
        s,f = scale.get(p,(None,None));
        if s is not None:
          matrix_scale1(dd,s);

## -------------------- other init
    if self.gain_smoothing:
      dprint(0,"a Gaussian smoothing kernel of size",self.gain_smoothing,"will be applied");
    if self.gain_bounds:
      dprint(0,"gains will be flagged on amplitudes outside of",self.gain_bounds);
    dprint(1,"initial gain value is",self._init_value_gain.values()[0].flat[0] if isinstance(self._init_value_gain,dict) else
      self._init_value_gain);
    # init gain parms object
    gain = GainClass(datashape,expanded_datashape,gain_subtiling,solvable_ifrs,
              self.epsilon,self.convergence_quota,smoothing=self.gain_smoothing,
              bounds=self.gain_bounds,
              init_value=self._init_value_gain,
              verbose=self.verbose);

    if self.print_variance:
      print_variance(variance);

## -------------------- init diffgains
    if num_diffgains:
      dprintf(0,"also solving for %d differential gains\n",num_diffgains);
      if self.diffgain_smoothing:
        dprint(0,"a Gaussian smoothing kernel of size",self.diffgain_smoothing,"will be applied");
      for i in range(num_diffgains):
        initval = self._init_value_dg.get(i,self.init_value);
        dprint(1,"initial gain value #%d is"%i,initval.items()[0] if isinstance(initval,dict) else initval);
      diffgains = [
        GainClass(datashape,expanded_datashape,dg_subtiling,solvable_ifrs,
          self.diffgain_epsilon,self.diffgain_convergence_quota,
          smoothing=self.diffgain_smoothing,
          init_value=self._init_value_dg.get(i,self.init_value),
          verbose=self.verbose)
        for i in range(num_diffgains) ];
    else:
      diffgains = [];

    #model_unscaled,data_unscaled = model,data;
    ## send to phase center
    #corr = {};
    #for p,q in data.keys():
      #mm = model[p,q];
      #phase = numpy.angle(mm[0]*mm[3])/2;
      #corr = numpy.exp(-1j*phase);
      #model[p,q] = matrix_scale(mm,corr);
      #data[p,q] = matrix_scale(data[p,q],corr);

    ## apply depolarisation scaling
#    dps = self.compute_depol_scaling(solvable_antennas,model);
#    model = self.apply_scaling(dps,model);
#    data = self.apply_scaling(dps,data);

## -------------------- start of major loop
    gain_iterate = (data,model) if self.gains_on_data else (model,data);
    # start major loop -- alternates over gains and diffgains
    for nmajor in range(self.max_major+1):
      if (domain_id == self.dump_domain or self.dump_domain == -1):
        dump_data_model(data,model,solvable_ifrs,"dump_G%d.txt"%nmajor);
      ## -------------------- first, iterate normal gains to convergence
      # setup convergence criteria
      if nmajor == 0:
        maxiter,delta = self.max_iter,self.delta;
      elif nmajor < self.max_major:
        maxiter,delta = self.max_iter_1,self.delta_1;
      else:
        maxiter,delta = self.max_iter_2,self.delta_2;
      # run solution loop
      self.solve_gains("G",gain,gain_iterate,(data,model,self.gains_on_data),
                      noise,weight,flagmask,flags,
                      maxiter,self.max_diverge,delta,
                      averaging=(self.average_1 if nmajor else self.average),
                      omega=(self.omega_1 if nmajor else self.omega),
                      feed_forward=(self.feed_forward_1 if nmajor else self.feed_forward));
      ## -------------------- now iterate over diffgains
      # break out if no diffgains to iterate over, or if we're on the last major cycle
      if not num_diffgains or nmajor >= self.max_major:
        break;
      else:
        # subtract this from corrected data: D1 = G*D*G^H - M0 - corrupt(M1) - corrupt(M2) - ... to obtain residuals
        # (G may need to be inverted depending on mode)
        if self.gains_on_data:
          data1 = dict([ (pq,gain.residual(data,model,pq)) for pq in solvable_ifrs ]);
        else:
          data1 = dict([ (pq,gain.residual_inverse(data,model,pq,
              regularize=self.regularization_factor if self.regularize_intermediate else 0))
              for pq in solvable_ifrs ]);
          ## check for NANs in the data
          self.check_finiteness(data1,"corrected data");
              
        if self.weigh_diffgains:
          resnoise,resweight = self.compute_noise(data1,flagmask);
          if _verbosity.verbose > 3:
            for pq,bl in self.ifr_by_baseline:
              rms = resnoise.get(pq);
              if rms is not None:
                dprint(4,"%20s %.2fm"%("-".join(pq),bl),"residual std %.2g weight %.2f"%
                      (rms,resweight.get(pq,0)));
        else:
          resweight = resnoise = None;
#        for pq in list(self._solvable_ifrs)[:1]:
#          print [ (pq,[ d1 if is_null(d1) else abs(d1).max() for d1 in data1[pq] ]) for pq in self._solvable_ifrs ];
        # reset current model to a _copy_ of model0 -- each DG term will be added to it (in place) at the end of
        # each DG loop iteration
        for pq in solvable_ifrs:
          model[pq] = matrix_copy(model0[pq]);
        # now loop over all diffgains and iterate each set once
        for idg,dg in enumerate(diffgains):
          # add current estimate of corrupt(Mi) back into data1, and subtract from current model
          for pq in solvable_ifrs:
            corr = dg.apply(dgmodel[idg],pq,cache=True);
            dd = data1[pq];
            for n,c in enumerate(corr):
              dd[n] += c;
#          print data1['0','A'][0][0,0],model0['0','A'][0][0,0],model['0','A'][0][0,0],dgmodel[idg]['0','A'][0][0,0];
#          print "d1",data1['0','A'][0],"m1",dgmodel[idg]['0','A'][0];
          if ( idg == self.dump_diffgain or self.dump_diffgain == -1 ) and \
             ( domain_id == self.dump_domain or self.dump_domain == -1 ):
            dump_data_model(dgmodel[idg],data1,solvable_ifrs,"dump_E%d-%d.txt"%(idg,nmajor));
          # iterate this diffgain solution
          self.check_finiteness(data1,"data1 after DG%d added in"%idg);
          self.solve_gains("dE:%d"%idg,dg,(dgmodel[idg],data1),(data1,dgmodel[idg],False),
                          resnoise,resweight,flagmask,
                          flags=None,
                          maxiter=self.diffgain_max_iter,max_diverge=self.diffgain_max_diverge,
                          delta=self.diffgain_delta,averaging=self.average_de,omega=self.omega_de,
                          feed_forward=self.feed_forward_de);
          # add to model, and subtract back from data1 if needed
          for pq in solvable_ifrs:
            corr = dg.apply(dgmodel[idg],pq,cache=True);
            d1 = data1[pq];
            mm = model[pq];
            # serious bug here -- things weren't added back to the model properly. Or were they?
            for i,c in enumerate(corr):
              if is_null(mm[i]):
                mm[i] = c;
              else:
                mm[i] += c;
              if idg<num_diffgains-1:
                d1[i] -= c;
        # just in case, reset flagged model values back to 0
        for pq,fm in flagmask.iteritems():
          if pq in model:
            for m in model[pq]:
              if not is_null(m):
                if m[fm].max() != 0:
                  print "Ooops,",pq,"ended up with a non-zero flagged model or data!";
                m[fm] = 0;

    # if we were solving for diffgains, then model is not completely up-to-date, since the non-solvable baselines
    # have been ignored. Fill them in here. Also, reset residuals
    if num_diffgains:
      for pq in set(valid_ifrs)-set(solvable_ifrs):
        mm = model[pq] = model0[pq];
        for idg,dg in enumerate(diffgains):
          for i,c in enumerate(dg.apply(dgmodel[idg],pq,cache=True)):
            mm[i] += c;

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
        dd = data.get(pq);
        if dd is None:
          continue;
        if self.gains_on_data:
          mm = gain.apply_inverse(model,pq,cache=True,regularize=self.regularization_factor);
        else:
          mm = gain.apply(model,pq,cache=True);
        for num,(d,m) in enumerate(zip(dd,mm)):
          # skip off-diagonal elements when in diagonal mode
          if self.diag_ifr_gains and num in (1,2):
            continue;
          # work out update to ifr gains
          if numpy.isscalar(m):
            m = numpy.array(m);
          if numpy.isscalar(d):
            d = numpy.array(d);
          dh = numpy.conj(d);
          if self.per_chan_ifr_gains:
            mdh = (m*dh).sum(0);
            ddh = (d*dh).sum(0);
          else:
            mdh = (m*dh).sum();
            ddh = (d*dh).sum();
          sri = self.ig_sum_reim[pq][num] = self.ig_sum_reim[pq][num] + mdh;
          ssq = self.ig_sum_sq[pq][num]   = self.ig_sum_sq[pq][num] + ddh;
          if numpy.isscalar(ssq):
            if ssq != 0:
              self.ifr_gain_update[pq][num] = sri/ssq;
          else:
            if (ssq!=0).any():
              self.ifr_gain_update[pq][num] = sri/ssq;
              self.ifr_gain_update[pq][num][ssq==0] = 1; 
#          if num == 0 and pq[0] == '0':
#           print m[0,0],d[0,0],dh[0,0]
#            print pq,(m*dh).sum(),(d*dh).sum(),sri/ssq;

    # work out result -- residual or corrected visibilities, depending on our state
    variance = {};
    nvells = maxres = 0;
    for pq in self._ifrs:
      m = model.get(pq);
      if m is None:
        for i in range(4):
          vs = datares.vellsets[nvells];
          if getattr(vs,'value',None) is not None:
            fl = getattr(vs,'flags',None);
            if fl is None:
              fl = vs.flags = meq.flags(datashape);
            fl[...] |= self.flagmask;
          nvells += 1;
        continue;
      else:
        if self.gains_on_data:
          out = res = gain.residual(data,model,pq);
          if not self.residuals:
            out = gain.apply(data,pq);
          elif not self.correct:
            out = matrix_negate(gain.residual_inverse(model,data,pq));
        else:
          out = res = gain.residual_inverse(data,model,pq,regularize=self.regularization_factor);
          if not self.residuals:
            out = gain.apply_inverse(data,pq,regularize=self.regularization_factor);
          elif not self.correct:
            out = matrix_negate(gain.residual(model,data,pq));
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
          if fl4 is not None:
            fl = getattr(vs,'flags',None);
            if fl is None:
              fl = vs.flags = meq.flags(datashape);
            fl[fl4 if expanded_dataslice is None else fl4[expanded_dataslice]] |=  self.flagmask;
          # compute stats
          maxres = max(maxres,abs(res[n]).max() if not numpy.isscalar(res[n]) else abs(res[n]));
          nvells += 1;

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
    dprint(0,"%s residual max %g, elapsed time %dm%0.2fs"%(
              request.request_id,maxres,m,s));

    return datares;

  def compute_noise (self,data,flagmask):
    """Computes delta-std and weights of data""";
    noise = {};
    weight = {};
    for pq,dd in data.iteritems():
      flag = flagmask.get(pq);
      dvalid = Ellipsis if flag is None else ~(flag[1:,...]|flag[:-1,...]);
      vv = [];
      for d in dd:
        if is_null(d):
          vv.append(0);
        else:
          delta = d[1:,...] - d[:-1,...];
          v = (delta.real[dvalid].std()+delta.imag[dvalid].std())/(2*math.sqrt(2));
          vv.append(float(v));
      # convert to weight
      if vv[1] or vv[2]:
        noise[pq] = (vv[1]+vv[2])/2;
        weight[pq] = (1/noise[pq])**2;
      elif vv[0] or vv[3]:
        noise[pq] = (vv[0]+vv[3])/2;
        weight[pq] = (1/noise[pq])**2;
    # normalize weights
    if weight:
      meanweight = sum(weight.itervalues())/len(weight);
      for pq,w in weight.iteritems():
        weight[pq] = w/meanweight;
    return noise,weight;

  def compute_chisq (self,gain,data,model,gains_on_data,noise=None,flagmask={}):
    chisq0 = chisq1 = 0;
    nterms = 0;
    for pq in self._solvable_ifrs:
      if pq not in data:
        continue;
      if gains_on_data:
        res = gain.residual(data,model,pq);
      else:
        res = gain.residual(model,data,pq);
      rms = noise.get(pq,0) if noise else 1;
      if not rms or not numpy.isfinite(rms):
        w = 1;
      else:
        w = rms**(-2);
#      fmask = reduce(operator.or_,[(d==0)&(m==0) for (d,m) in zip(data[pq],model[pq])]);
      fmask = flagmask.get(pq);
      n = self._datasize if fmask is None else (self._expanded_size - fmask.sum());
      for ir,r in enumerate(res):
        fin = numpy.isfinite(r);
        if not fin.all():
          dprintf(4,"%s element %d: %d slots are INF/NAN, omitting from chisq sum\n",pq,ir,fin.size-fin.sum());
        elif not is_null(r):
          chisq0 += (r*numpy.conj(r)*w).real.sum();
          chisq1 += (r*numpy.conj(r)).real.sum();
          nterms += n;
    return float(chisq0)/nterms,float(chisq1)/nterms;

  def compute_depol_scaling (self,antennas,model):
    scaling = {};
    #for p in antennas:
      #msum = [0j,0j,0j,0j];
      #for q in antennas:
        #if (p,q) in model:
          #m = model[p,q];
        #elif (q,p) in model:
          #m = matrix_conj(model[q,p]);
        #else:
          #continue;
        #for i in range(4):
          #msum[i] += m[i];
      #if msum[0] is not 0j and msum[3] is not 0j:
        #scaling[p] = [1,-msum[1]/msum[3],-msum[2]/msum[0],1];
    for p in antennas:
      modsum = NULL_MATRIX();
      modsqsum = NULL_MATRIX();
      have_scale = False;
      for q in antennas:
        if (p,q) in model:
          m = model[p,q];
          mh = matrix_conj(m);
        elif (q,p) in model:
          mh = model[q,p];
          m = matrix_conj(mh);
        else:
          continue;
        matrix_add1(modsum,mh);
        matrix_add1(modsqsum,matrix_multiply(mh,m));
        have_scale = True;
      if have_scale:
        scaling[p] = matrix_multiply(modsum,matrix_invert(modsqsum));
    return scaling;

  def apply_scaling (self,scaling,data):
    result = {};
    for (p,q),d in data.iteritems():
      scale = scaling.get(p);
      result[p,q] = matrix_multiply(scale,d) if scale is not None else d;
      dh = matrix_conj(d);
      scale = scaling.get(q);
      result[q,p] = matrix_multiply(scale,dh) if scale is not None else dh;
    return result;
    
  def check_finiteness (self,data,label,complete=False):
    # check for NANs in the data
    for pq,dd in data.iteritems():
      for i,d in enumerate(dd):
        if type(d) is numpy.ndarray:
          fin = numpy.isfinite(d);
          if not fin.all():
            dprintf(0,"%s %s element %d: %d slots are INF/NAN\n",label,pq,i,d.size-fin.sum());
            if not complete:
              break;


  def solve_gains (self,label,gain,iter_args,chisq_args,noise,weight,flagmask,flags=None,
                  maxiter=10,max_diverge=2,delta=1e-6,averaging=2,omega=.5,feed_forward=False):
    """Runs a single gain solution loop to completion"""
    gain_dchi = [];
    gain_maxdiffs = [];
    # initial chi-sq, and fallback chisq for divergence
    gain._reset();
    init_chisq,init_chisq_unnorm = self.compute_chisq(gain,noise=noise,flagmask=flagmask,*chisq_args);
    lowest_chisq = (init_chisq,gain.get_values());
    num_diverged = 0;
    dprint(2,"solving for %s, initial chisq is %.12g"%(label,init_chisq));
    t0 = time.time();
    chisq0 = init_chisq;
    # iterate
    for niter in range(maxiter):
      # iterate over normal gains
      converged,maxdiff,deltas,nfl = gain.iterate(niter=niter,weight=weight,flagmask=flagmask,
        averaging=averaging,omega=omega,feed_forward=feed_forward,
        *iter_args);
      dprint(4,"%d slots converged, max delta is %f"%(gain.num_converged,float(gain.delta_max)));
      gain_maxdiffs.append(float(maxdiff));
      if self.visualize_gains > 1:
        global_gains.setdefault(label,[]).append(gain.get_2x2_gains(datashape,expanded_dataslice));
      if nfl:
        dprint(2,"%d out-of-bound gains flagged at iteration %d"%(nfl,niter));
      ## check chi-square
      chisq,chisq_unnorm = self.compute_chisq(gain,noise=noise,flagmask=flagmask,*chisq_args);
      dchi = (chisq0-chisq)/chisq;
      gain_dchi.append(dchi);
      dprint(3,"iter %d ||%s||=%g max update is %g chisq is %.8g (%.8g) diff %.3g"%(niter+1,label,float(gain.gainnorm),float(gain.delta_max),chisq,chisq_unnorm,dchi));
      chisq0 = chisq;
      # check convergence
      if converged:
        dprint(1,"%s converged at chisq %.12g (last gain update %g) after %d iterations and %.2fs"%(
                label,chisq,float(gain.delta_max),niter+1,time.time()-t0));
        break;
      # if chi-sq decreased, remember this
      if dchi >= 0:
        lowest_chisq = (chisq,gain.get_values());
        num_diverged = 0;
        # and check for chisq-convergence
        if niter > 1 and dchi < delta:
          dprint(1,"%s chisq converged at %.12g (last gain update %g) after %d iterations and %.2fs"%(
                label,chisq,float(gain.delta_max),niter+1,time.time()-t0));
          break;
      # else chisq is increasing, check if we need to stop
      else:
        num_diverged += 1;
        if num_diverged >= max_diverge:
          chisq,gainvals = lowest_chisq;
          gain.set_values(gainvals);
          dprint(1,"%s chisq diverging, stopping at %.12g after %d iterations and %.2fs"%(
                    label,chisq,niter+1,time.time()-t0));
          break;
    else:
      dprint(1,"%s max iterations (%d) reached at chisq %.12g (last gain update %g) after %.2fs"%(
              label,maxiter,chisq,float(gain.delta_max),time.time()-t0));
    # check if we have a lower chisq to roll back to
    if chisq > init_chisq:
      dprint(1,"DANGER WILL ROBINSON! Final chisq is higher than initial chisq. Rolling back to initial values.");
      gain.set_values(gainvals);
    elif chisq > lowest_chisq[0]:
      chisq,gainvals = lowest_chisq;
      gain.set_values(gainvals);
    # print other stats
    dprint(2,"  delta-chisq were"," ".join(["%.4g"%x for x in gain_dchi]));
    dprint(2,"  convergence criteria were"," ".join(["%.2g"%x for x in gain_maxdiffs]));
    if flags is not None:
      gainflags = getattr(gain,'gainflags',None);
      if gainflags:
        dprint(2,"total gain-flags: "," ".join(["%s:%d"%(p,gf.sum()) for p,gf in sorted(gainflags.iteritems())]));
        for p in solvable_antennas:
          flag4 = gainflags.get(p);
          if flag4 is None:
            continue;
          for q in antennas:
            pq = (p,q) if (p,q) in data else ((q,p) if (q,p) in data else None);
            if pq:
              if pq in flags:
                flags[pq] |= flag4;
              else:
                flags[pq] = flag4;
