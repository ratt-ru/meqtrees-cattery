#from memory_profiler import profile

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
import gc
import scipy.ndimage.measurements

from MatrixOps import *
import DataTiler

_verbosity = Kittens.utils.verbosity(name="stefcal");
dprint = _verbosity.dprint;
dprintf = _verbosity.dprintf;

## bitflag constants used below
FPRIOR  = 1;  # prior flags
FINSUFF = 2;  # insufficient data for solution
FNOCONV = 4;  # no convergence
FCHISQ  = 8;  # chisq too high
FSOLOOB = 16; # solution out of bounds

DEBUG_SLICE = (slice(0,10), 0)

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
    mystate('norm_offdiag',False);
    mystate('label','G');
    mystate('index',[]);
    self.set_symdeps("Domain");

#  @profile
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
        if self.norm_offdiag:
          xy,yx = xy/xx,yx/yy
        vellsets += [ gain_to_vellset(ix,x) for ix,x in enumerate((xx,xy,yx,yy)) ];

    res = meq.result(cells=request.cells);
    res.vellsets = vellsets;
    if self.index == []:
      res.dims = [ len(keys),nsets,2,2 ] if nsets>1 else [ len(keys),2,2 ];
    else:
      res.dims = [ len(keys),2,2 ];
    return res;

    
from GainOpts import *
    
# TERMINOLOGY
# self._datashape:           original shape of input data, e.g. (NT,NF)
#       self._datasize       corresponding # of slots (=NT*NF)
# self._expanded_datashape:  datashape padded out to nearest multiple of solution tiling (NT0,NF0)
#       self._expanded_size: corresponding # of slots (=NT0*NF0)
# self._expanded_dataslice: slice applied to array of shape expanded_datashape to extract the original datashape,
#                           e.g. (slice(0,NT),slice(0,NF))

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
    # use cross-hand data to estimate noise
    mystate('use_polarizations_for_noise',False);
    # compute noise estimates per-channel
    mystate('noise_per_chan',True);
    # verbosity level
    mystate('verbose',0);
    # verbosity level
    mystate('critical_flag_threshold',20);
    # number of diffgains
    mystate('diffgain_labels',[]);
    # init gain objects
    self.gain  = GainOpts("","gain","G");
    self.bgain  = GainOpts("","gain1","B");
    self.gainopts_all = [ self.gain,self.bgain ];
    # update state of all diffgain objects, and make list of enabled ones
    for opt in self.gainopts_all:
      dprint(2,"initializing options for",opt.label);
      opt.update_state(self,mystate,verbose=self.verbose);
    # self.gainopts: list of enabled DI gain objects
    self.gainopts = [ opt for opt in self.gainopts_all if opt.enable ];
    # init diffgain objects, one per source
    self.dgopts = [];
    for i,label in enumerate(self.diffgain_labels):
      dg = GainOpts("","diffgain","dE");
      # init each diffgain first with the diffgain_xxx vars, then with diffgain_xxx_LABEL
      dprint(2,"initializing options for",dg.label);
      dg.update_state(self,mystate,verbose=self.verbose);
      # this will have overwritten the label, so reset
      dg.label += ":%s"%label;
      dprint(2,"initializing specific options for",dg.label);
      dg.update_state(self,option_suffix=label);
      if dg.enable:
        self.dgopts.append(dg);
    self.use_float_di = all([gg.use_float for gg in self.gainopts])
    self.use_float_dd = all([dg.use_float for dg in self.dgopts])
    dprintf(2,"using float di %s dd %s\n",self.use_float_di,self.use_float_dd)
    # we're polarized if at least one gain is polarized
    self.polarized = any([ opt.polarized for opt in self.gainopts+self.dgopts ]);
    # roll back solutions if final chisq exceeds initial chi-sq
    mystate('chisq_rollback',True);
#    # labels for gain, ifr gain and differential gain parameters
#    mystate('gain_parm_label',"G");
#    mystate('ifr_gain_parm_label',"IG");
#    mystate('diffgain_parm_label',"dE");
    mystate('num_major_loops',10);
    # downsampling settings
    mystate('downsample_output',True);
    mystate('downsample_subtiling',[]);
    # use stored solution (if available) as starting guess
    mystate('init_from_table',True);
    # use previous tile (timeslot) as starting guess -- if table not available
    mystate('init_from_previous',True);
    # rescale data to model
    mystate('rescale',True);
    # use this value as starting guess -- if previous two not available
    mystate('init_value',1);
    # regularization factor applied to gain solutions for correction
    mystate('regularization_factor',0);
    # regularize intermediate corrections (when solving for dEs)?
    mystate('regularize_intermediate',False);
    # name of gain tables to which solutions are saved (or from which they are loaded)
    mystate('ifr_gain_table','ifrgains.cp');
    # filenames for solutions
    # return residuals (else corrected data)
    mystate('residuals',True);
    # subtract sources with diffgains, if returning corrected data
    mystate('subtract_dgsrc',False);
    # flagbit to assign to flagged visiblities
    mystate('output_flag_bit',0);
    # apply IFR gains (if False, then none of the other IFR gain-related options apply)
    mystate('apply_ifr_gains',False);
    # solve for IFR gains (if False, then simply load & apply)
    mystate('solve_ifr_gains',True);
    # save solutions
    mystate('save_ifr_gains',True);
    # ignore loaded solutions 
    mystate('reset_ifr_gains',True);
    if not self.apply_ifr_gains:
      self.solve_ifr_gains = self.save_ifr_gains = self.reset_ifr_gains = False;
    # IFR gains are per-frequency, or one value across entire band
    mystate('per_chan_ifr_gains',False);
    # IFR gains are for diagonal elements only if True, for full 2x2 if False
    mystate('diag_ifr_gains',False);
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
    self._init_value_gain = self._init_value_bgain = self.init_value;
    self._init_value_dg = {};

#  @profile
  def get_result (self,request,*children):
    dprint(1,"get_result entry");
    timestamp0 = time.time();
    # get dataset ID from request
    dataset_id,domain_id = meq.split_request_id(request.request_id);
    # get domain ID from request
    time0,time1,timestep,numtime,freq0,freq1,freqstep,numfreq = request.cells.domain.domain_id;
    # child 0 is data
    # child 1 is direction-independent model
    # children 2 and on are models subject to dE terms
    num_diffgains = len(self.dgopts);
    if num_diffgains != len(children)-2:
      raise TypeError,"StefCalNode: %d children (data,model,model1,...) must be provided"%(num_diffgains+2);
    
    # if new dataset ID, do setup for start of new dataset
    if dataset_id != self._dataset_id:
      self._dataset_id = dataset_id;
 
      # try to load gain solutions if available
      for opt in self.gainopts + self.dgopts:
        opt.load_initval(self.init_value);
      GainOpts.flush_tables();
      
      dprint(1,"new dataset id",dataset_id);
      # if asked to solve for IFR gains, set up dicts for collecting stats
      if self.solve_ifr_gains:
        self.ig_sum_reim = dict([ (pq,[0j]*4) for pq in self._ifrs ]);
        self.ig_sum_sq   = dict([ (pq,[0.]*4) for pq in self._ifrs ]);
        self.ifr_gain_update = dict([ (pq,[1.]*4) for pq in self._ifrs ]);
      # read previous IFR gains from table, if asked to apply them
      self.ifr_gain = {};
      if self.apply_ifr_gains and not self.reset_ifr_gains and os.path.exists(self.ifr_gain_table):
        try:
          self.ifr_gain = cPickle.load(file(self.ifr_gain_table));
          dprint(1,"loaded %d ifr gains from %s"%(len(self.ifr_gain),self.ifr_gain_table));
          # reset off-diagonals to 1
          if self.diag_ifr_gains:
            dprint(1,"resetting off-diagonal elements to 1");
            for gg in self.ifr_gain.itervalues():
              gg[1] = gg[2] = 1;
        except:
          traceback.print_exc();
          dprint(0,"error loading ifr gains from",self.ifr_gain_table);

    # child 0 is data
    # child 1 is direction-independent model
    # children 2 and on are models subject to dE terms

    # check inputs and populate mappings
    pqij_all = [];      # list of all (p,q),i,j tuples
    pqij_data = [];     # subset of (p,q),i,j tuples for which we have non-null input
    pqij_solvable = []; # subset of (p,q),i,j tuples for which we solve for gains
    data  = {};         # mapping from (p,q) to four data time-freq planes
    model0 = {};        # mapping from (p,q) to four model (M0) time-freq planes
    dgmodel = [ {} for i in range(num_diffgains) ];
                        # for each diff gain, mapping from (p,q) to M1,M2,... model time-freq planes (4 each)
    dgmodel_corr = [ {} for i in range(num_diffgains) ];
                        # diffgain model, corrupted with appropriate diffgain term
    model = {};         # this is the full model, M0+M1+M2+...

    # This will contain a mask of flagged values. I would use masked arrays, but they seem to slow something down, so no no. 
    # Instead, we'll use 0.0 for missing values in model and data (since zeroes do not upset the equations), and maintain a 
    # bitflags array for masking things out.
    # Need bitflags rather than a single flag so that we can distinguish what the origin of the flag was 
    # (and also because we may clear flags in the first cycle of the major loop)
    bitflags = {};

    antennas = set([p for p,q in self._ifrs]) | set([q for p,q in self._ifrs]);

    # per-baseline noise
    variance = {};
    # antennas for which we have non-trivial data
    solvable_antennas = set();
    # this will count the valid visibilities per each antenna, per each time/freq slot
    vis_per_antenna = None;

    datares = children[0]
    data_dims = getattr(datares,'dims',None);
    if data_dims is None:
      raise TypeError,"No data dimensions. Have you specified a valid input column?";
    modelres = children[1];
    if any( [ getattr(ch,'dims',None) != data_dims for ch in children[1:] ] ):
      raise TypeError,"Dimensions of data and model(s) do not match. Have you specified a sky model?";
    # expecting Nx2x2 matrices
    if len(datares.dims) == 3:
      if datares.dims[1] != 2 or datares.dims[2] != 2:
        raise TypeError,"Data and model must be of rank Nx2x2";
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
          # get data
          d = getattr(datares.vellsets[nvells],'value',0);
          # apply ifr gains if we have them
          if not is_null(d):
            # convert if needed
            if self.use_float_di:
              d = d.astype(numpy.complex64)  # make copy by default
            if d.flags['WRITEABLE']:
              d *= ifrgain[num]
            else:
              d = d*ifrgain[num]
          # get model
          m = getattr(modelres.vellsets[nvells],'value',0);
          if hasattr(datares.vellsets[nvells],'flags'):
            flags = (datares.vellsets[nvells].flags != 0);
          else:
            flags = None;
          # does this need to be skipped? only process data otherwise
          if not ( is_null(d) if self.polarized else (is_null(m) or is_null(d)) ):
            # if this is the first datum, then check shape, and prepare subtilings etc.
            # for the first valid result, setup shapes and stuff
            if not model0:
              def get_dtype (dd):
                if dd:
                  return numpy.complex64 if self.use_float_dd else numpy.complex128
                else:
                  return numpy.complex64 if self.use_float_di else numpy.complex128
              # this is the basic time-frequency shape
              self._datashape = datashape = tuple(d.shape);
              self._datasize = reduce(operator.mul,datashape);
              # figure out subtiling
              self._expanded_datashape = expanded_datashape = GainOpts.resolve_tilings(datashape,*(self.gainopts+self.dgopts));
              # if tiling does not tile the data shape perfectly, we'll need to expand the input arrays
              # Define pad_array() as a function for this: it will set to be identity if no expansion is needed
              if datashape != expanded_datashape:
                self._expanded_dataslice = expanded_dataslice = tuple([ slice(0,nd) for nd in datashape ]);
                def pad_array (x,initval=0,dd=False):
                  if is_null(x):
                    return 0;
                  x1 = numpy.empty(expanded_datashape, dtype=bool if type(initval) is bool else get_dtype(dd))
                  x1[...] = initval;
                  x1[expanded_dataslice] = x;
                  return x1;
                self._expanded_size = reduce(operator.mul,expanded_datashape);
                self._expansion_ratio = self._expanded_size/float(self._datasize);
                self._expansion_mask = numpy.zeros(expanded_datashape,bool);
                self._expansion_mask[expanded_dataslice] = True;
                dprint(1,"input arrays will be expanded to shape",expanded_datashape,"ratio %.2f"%self._expansion_ratio);
              else:
                self._expanded_size = self._datasize;
                self._expansion_ratio = 1;
                self._expanded_dataslice = expanded_dataslice = None;
                self._expansion_mask = numpy.ones(datashape,bool);
                def pad_array (x,initval=0,dd=False):
                  if is_null(x):
                    return x
                  elif type(initval) is bool:
                    return x
                  else:
                    return x.astype(get_dtype(dd))  # copy=True implicitly
              # this counts how many valid visibilities we have per each antenna, per each time/freq slot
              vis_per_antenna = dict([(p,numpy.zeros(expanded_datashape,dtype=int)) for p in antennas ]);
            # now check inputs and add them to data and model dicts
            if d.shape != datashape:
              raise TypeError,"data shape mismatch at %s:%s:%s:%s, %s vs %s" % (pq[0], pq[1],
                self.corr_names[i], self.corr_names[j], d.shape, datashape )
            if not is_null(m) and m.shape != datashape:
              raise TypeError,"model shape mismatch at %s:%s:%s:%s, %s vs %s" % (pq[0], pq[1],
                self.corr_names[i], self.corr_names[j], m.shape, datashape )
            # add to data/model matrices, applying the padding function defined above
            m0 = model0.setdefault(pq,[0,0,0,0])[num] = pad_array(m);
            d0 = data.setdefault(pq,[0,0,0,0])[num] = pad_array(d);
            # apply flags
            if flags is not None:
              flags = pad_array(flags,True);
              if not is_null(m0):
                m0[flags] = 0;
              if not is_null(d0):
                d0[flags] = 0;
              invalid = (d0==0)&(m0==0);
              self.add_flags(bitflags,pq,invalid*FPRIOR);
            # if there are at least some valid points on this baseline, add it to solvable antennas
            if (bitflags.get(pq) is None) or not bitflags[pq].all():
              if pq in self._solvable_ifrs:
                solvable_antennas.update(pq);
  #              print pq,validmask[pq];
            # get models for dE-subjected terms
            if num_diffgains:
              for k in range(num_diffgains):
                m1 = children[2+k].vellsets[nvells].value
                m1 = pad_array(m1,dd=True);
                if flags is not None and not is_null(m1):
                  m1[flags] = 0;
                dgmodel[k].setdefault(pq,[0,0,0,0])[num] = m1;
        # ok, done looping over the 2x2 visibility matrix elements. 
        # If we have found anything valid at all, finalize flagmasks etc.
        if pq in model0:
          # look at bitflags to see how many valid correlations we have, and zero the flagged ones
          fmask = bitflags.get(pq);
          if fmask is not None:
            fmask = fmask!=0;
            valid = self._expanded_size - fmask.sum();
            if valid > 0:
              dprint(4,"%s-%s"%pq,"has %d of %d unflagged correlation matrices"%(valid,self._datasize));
              for dataset in [data,model0] + dgmodel:
                for x in dataset.get(pq,[]):
                  if not is_null(x):
                    x[fmask] = 0;
              validmask = (~fmask).astype(int);
              vis_per_antenna[pq[0]] += validmask;
              vis_per_antenna[pq[1]] += validmask;
            else:
              # if nothing is valid, remove baseline from dicts
              dprint(4,"%s-%s"%pq,"is completely flagged, skipping");
              for dataset in [data,model0] + dgmodel:
                del dataset[pq];
          else:
            dprint(4,"%s-%s"%pq,"has no flagged correlation matrices, all data is valid");
            vis_per_antenna[pq[0]] += 1;
            vis_per_antenna[pq[1]] += 1;
    else:
      # in principle could also handle [N], but let's not bother for now
      raise TypeError,"data and model must be of rank Nx2x2";

    # hang onto datares record since we'll be putting the results into it
    # release modelres and all the other child results (they're already held in model and dgmodel)
    # use resize to explicitly release the memory since we KNOW nobody else is using it
    dprint(1,"constructed internal arrays, trying to release array memory");
    modelres = children = None
    gc.collect()
    dprint(1,"released memory");
    
    valid_ifrs = data.keys();
    solvable_ifrs = set(self._solvable_ifrs)&set(valid_ifrs);
    if not solvable_ifrs:
      dprint(1,"no valid data found for solvable IFRs  -- nothing to stefcal!");
      return datares;
    dprint(1,"Found %d solvable antennas"%len(solvable_antennas));
    dprint(2,"  valid ifrs outside the solvable set:"," ".join(["%s-%s"%pq for pq in set(valid_ifrs)-set(self._solvable_ifrs)]));
    dprint(2,"  ifrs with no data:"," ".join(["%s-%s"%pq for pq in set(self._ifrs)-set(valid_ifrs)]));
    min_baselines_for_solution = len(solvable_antennas)+1;
    
    # pq00 set for debugging purposes, first valid IFR
    pq00 = sorted(valid_ifrs)[0]


## -------------------- downsample data and model, if needed
    downsample_subtiling = self.downsample_subtiling;
    downsampler = None;
    if downsample_subtiling:
      downsample_subtiling = [ max(d,1) for d in downsample_subtiling ];
      for iaxis,ds in enumerate(downsample_subtiling):
        for opt in self.gainopts+self.dgopts:
          if opt.subtiling[iaxis]%ds != 0:
            raise RuntimeError,"axis %d: %s solution interval must be a multiple of downsample interval"%(iaxis,opt.name);
      if max(downsample_subtiling) == 1:
        downsample_subtiling = None;
      else:
        # create retiler for going from downsampled to full resolution
        downsampler = DataTiler.DataTiler(expanded_datashape,downsample_subtiling,original_datashape=datashape);
        # create retilers for going from gain tiling to full resolution
        for opt in self.gainopts+self.dgopts:
          opt.vis_tiler = DataTiler.DataTiler(expanded_datashape,opt.subtiling,original_datashape=datashape,force_subtiling=True);
          if not self.downsample_output:
            opt.tiler = opt.vis_tiler;
          # update option settings
          opt.smoothing = [ ds/float(st) for ds,st in zip(opt.smoothing,downsample_subtiling) ];
          opt.subtiling = [ gs//st for gs,st in zip(opt.subtiling,downsample_subtiling) ];
        # keep copies of model and data at original sampling
        orig_sampled_data = data.copy();
        orig_sampled_bitflags = bitflags.copy();
        orig_sampled_model0 = model0.copy();
        orig_sampled_dgmodel = [ dg.copy() for dg in dgmodel ];
        # resample
        downsample_factor = reduce(operator.mul,downsample_subtiling);
        dprint(1,"resampling data by a factor of %d=%s"%(downsample_factor,"x".join(map(str,downsample_subtiling))));
        for pq in data.keys():
          flags = bitflags.get(pq);
          # resample flags, and compute number of valid slots per resampled interval, and a norm based on this
          if flags is not None and not numpy.isscalar(flags):
            nv = downsample_factor - downsampler.reduce_tiles(downsampler.tile_data(flags!=0));
            fl = bitflags[pq] = FPRIOR*(nv==0);
            norm = numpy.where(fl,0,1./nv);
          else:
            norm = 1./downsample_factor;
          # resample data
          for vissets in [data,model0] + dgmodel:
            dd = vissets.get(pq);
            if dd is not None:
              vissets[pq] = [ downsampler.reduce_tiles(downsampler.tile_data(d))*norm if 
                              d is not None and not numpy.isscalar(d) else d for d in dd ];
        # change other settings
        orig_sampled_expanded_datashape = expanded_datashape;
        orig_sampled_datashape = datashape;
        orig_expansion_mask = self._expansion_mask;
        self._expansion_mask = numpy.zeros(expanded_datashape);
        self._datashape = datashape = [ int(math.ceil(ds/float(st))) for ds,st in zip(datashape,downsample_subtiling) ];
        self._expansion_mask[tuple([ slice(0,nd) for nd in datashape ])] =True;
        self._expanded_datashape = expanded_datashape = [ ds/st for ds,st in zip(expanded_datashape,downsample_subtiling) ];
        self._datasize /= downsample_factor;
        self._expanded_size /= downsample_factor;
    
## -------------------- rescale data to model if asked to
    if self.rescale and self.rescale != "no":
      scale = {};
      finite = {};
      # compute scales as s(p) = ||sum_q Mpq||/||sum_q Dpq||
      for p in antennas:
        dsum = msum = 0;
        for q in antennas:
          d = m = None;
          if (p,q) in data:
            d,m = data.get((p,q)),model0.get((p,q));
          elif (q,p) in data:
            d,m = data.get((q,p)),model0.get((q,p));
          if d is None or m is None:
            continue;
          dsum += sum([x*numpy.conj(x) for x in d]);
          msum += sum([x*numpy.conj(x) for x in m]);
        if self.rescale == "scalar":
          dsum = dsum.sum() if dsum is not 0 else 0+0j;
          msum = msum.sum() if msum is not 0 else 0+0j;
          if dsum:
            scale[p] = numpy.power(msum.real/dsum.real,0.25),True;
          else:
            dprint(2,"no valid data (and thus no scale) for",p);
        else:
          if dsum is not 0:
            s = numpy.power(msum.real/dsum.real,0.25);
            f = numpy.isfinite(s);
          if dsum is 0 or (~f).all():
            dprint(2,"no valid data (and thus no scale) for",p);
          else:
            s[~f] = 0;
            scale[p] = s,f;
      # apply scales
      if scale:
        if self.rescale == "scalar":
          dprint(2,"per-antenna data scaling factors are ",", ".join(["%s %.3g"%(p,s) for p,(s,f) in scale.iteritems() ]));
        else:
          dprint(1,"min/max scaling factors are",min([s[f].min() for s,f in scale.itervalues()]),
                                                max([s[f].max() for s,f in scale.itervalues()]));
          dprint(2,"per-antenna data scaling factors are ",", ".join(["%s %.3g"%(p,s.max()) for p,(s,f) in scale.iteritems() ]));
      else:
        dprint(1,"rescaling not done, as none of the antennas appear to have any valid data");
      for (p,q),dd in data.iteritems():
        s1,f = scale.get(p,(None,None));
        s2,f = scale.get(q,(None,None));
        if s1 is not None and s2 is not None:
          matrix_scale1(dd,s1*s2);
          
## -------------------- compute the noise estimate, and weights based on this
    noise,weight = self.compute_noise(data,bitflags);
    nw = nm = nnw = 0;
    # print and total up stats, check for funny situations
    for pq,bl in self.ifr_by_baseline:
      if pq in data:
        w = weight.get(pq,0);
        rms = noise.get(pq,0);
        if not is_null(w):
          nw += 1;
        else:
          if is_null(data[pq][0]):
            nm += 1;
          else:
            fm = bitflags.get(pq);
            np = self._datasize if fm is None else (self._expanded_size - (fm!=0).sum());
            if np:
              dprint(1,"oops,","%s-%s"%pq,"has zero weights, yet %d valid points"%np);
              nnw += 1;
            else:
              nm += 1;
    dprint(1,"%d baselines with non-zero weights (%d with null weight, %d have no valid data)"%(nw,nnw,nm));
    if nw < min_baselines_for_solution:
      dprint(0,"At least %d valid baselines required for a solution! Aborting."%min_baselines_for_solution);
      raise RuntimeError,"not enough valid baselines";
      
    skip_solve = False;
    
  
## ----------------------- if solving for gains, check for required number of baselines per each t/f slot
## ----------------------- flag those that are missing
    solve_any = any([opt.solve for opt in self.gainopts+self.dgopts]);
        
    if solve_any:
      ##************** NB: move this to GAIN CLASS, as this needs to be done per subtile!
      ## this tells us how many elements there are in correlation matrix pq (per time/freq point)
      ## sum([0 if is_null(d) else ((d!=0)|(m!=0)).astype(int) for d,m in zip(data[pq],model[pq])])
      # this gives us the number of valid visibilities (in the solvable IFR set) per each time/freq slot
      nel_pq = dict([(pq,sum([0 if is_null(d) else ((d!=0)|(m!=0)).astype(int)
                      for d,m in zip(data[pq],model0[pq])])) for pq in data.iterkeys() ]);
      # this is how many valid visibility matrices we have total per slot
      num_valid_vis = sum([ int(n>1) if numpy.isscalar(n) else (n>1).astype(int) for pq,n in nel_pq.iteritems() ]);
      # we need as many as there are parameters per antenna, times ~4, to constrain the problem fully
      dprintf(1,"Found %d solvable antennas\n"%len(solvable_antennas));
      invalid_slots = (num_valid_vis>0)&(num_valid_vis<min_baselines_for_solution);
      dprintf(1,"%d valid visibility matrices per time/freq slot are required for a solution\n",min_baselines_for_solution);
      dprintf(1,"   for this data chunk, Nvalid ranges from %d to %d\n"%(num_valid_vis.min(),num_valid_vis.max()));
      nfl = nflagged_due_to_insufficient = invalid_slots.sum();
      if nfl:
        dprintf(1,"   %d of %d slots flagged due to insufficient visibilities\n",nfl,self._expanded_size);
        for pq in data.iterkeys():
          self.add_flags(bitflags,pq,invalid_slots*FINSUFF);
        for dataset in [data,model0] + dgmodel:
          for dd in dataset.itervalues():
            for d in dd:
              if not is_null(d):
                d[invalid_slots] = 0;
        if nfl == self._expanded_size:
          dprintf(0,"***WARNING: insufficient valid visibilities, no solutions will be attempted\n",nfl,self._expanded_size);
          skip_solve = True;
    else:
      skip_solve = True;
      nflagged_due_to_insufficient = 0;
      
        
## -------------------- init gain solvers
    if not skip_solve:
      dprintf(0,"Solvable: %d of %d inteferometers (%d have valid data), with %d solvable antennas\n",
        len(self._solvable_ifrs),len(self.ifrs),len(solvable_ifrs),len(solvable_antennas));
    for opt in self.gainopts+self.dgopts:
      opt.init_solver(datashape,expanded_datashape,solvable_ifrs,downsample_subtiling);

    if self.print_variance:
      print_variance(variance);

    dprint(2,"***DEBUG*** data",pq00,data[pq00][0][DEBUG_SLICE])
    dprint(2,"***DEBUG*** model0",pq00,model0[pq00][0][DEBUG_SLICE])
    if dgmodel:
      x = 0
      for dgm in dgmodel:
        x = x + dgm[pq00][0][DEBUG_SLICE] 
      dprint(2,"***DEBUG*** dgms",pq00,x)

    # now compute model is M0+M1+M2... 
    # where terms M1 etc. have dE's on them, if already initialized
    if len(self.dgopts):
      dprint(1,"adding dE-enabled terms into model");
      for pq,mod0 in model0.iteritems():
        mm = model[pq] = matrix_copy(mod0);
        for dg,dgm in zip(self.dgopts, dgmodel):
          if opt.has_init_value:
            cc = dg.solver.apply(dgm,pq,cache=False);
          else:
            cc = dgm[pq]
          for i,c in enumerate(cc):
            mm[i] += c;
    else:
      model = model0

    dprint(0,"***DEBUG*** model",pq00,model[pq00][0][DEBUG_SLICE])
    dprint(0,"***DEBUG*** model0",pq00,model0[pq00][0][DEBUG_SLICE])
            
    ## last check for NANs in the data and model
    self.check_finiteness(data,"data",bitflags);
    self.check_finiteness(model,"model",bitflags);
#    cPickle.dump(model,file("dump-model.cp","w"),2);

## -------------------- start of major loop
    initdata,initmodel,initweight = data,model,weight;
    if not skip_solve:
      # at beginning of each major loop, for an equation of the form D = G.B.(M0+dE1.M1+dE1^H+...)B^H.G^H
      #
      # initdata: contains the original data D
      # model0: contains the original M0
      # initmodel: contains M = M0+M1+M2+... i.e. without dE values, or with initial guess for dE values
      
      dprint(1,"resetting data and model to initial values");
      
      num_solvable = sum([ opt.solve for opt in self.gainopts ]);
      num_solvable_dd = sum([ opt.solve for opt in self.dgopts ]);
      
      if num_solvable<2 and not num_solvable_dd:
        max_major = 0;
      else:
        max_major = self.num_major_loops;
      
      # start major loop
      # do as many loops as specified, plus one more, to finalize the final solvable gain term
      for nmajor in range(max_major+1):
        last_loop = (nmajor == max_major);
        if last_loop:
          looptype = 2;
        else:
          looptype = 1 if nmajor else 0;
        
        if (domain_id == self.dump_domain or self.dump_domain == -1):
          dump_data_model(data,model,solvable_ifrs,"dump_G%d.txt"%nmajor);

        # at beginning of loop, model: contains M = M0+dE1.M1+dE1^H+... with the latest dE values
        # reset data and weight to initial values
        data,weight = initdata,initweight;
        flagged = False;  # will be True if new flags arise in DI terms
        pq0 = sorted(solvable_ifrs)[0]
        dprintf(2,"%s data type of data is %s\n",pq0,data[pq0][0].dtype)
        dprintf(2,"%s data type of model is %s\n",pq0,model[pq0][0].dtype)
        
        # at beginning of major loop 1, reset flags raised in loop 0.
        # This means that loop-0 flags are treated as preliminary-only
        if nmajor == 1:
          for bf in bitflags.itervalues():
            bf &= ~(FCHISQ|FNOCONV|FSOLOOB);
        
        ## make list of models that will be used at each stage of the direction-independent solve,
        ## by applying each DI term (from the innermost to the outermost)
        dimodels = [ model ];
        for opt in self.gainopts[-1:0:-1]:
          # modify model, if in second loop, or gain is already initialized
          if nmajor or opt.has_init_value:
            dprint(1,"applying %s %s to model"%(("current" if nmajor else "prior"),opt.label));
            model = dict([ (pq,opt.solver.apply(model,pq,cache=True)) for pq in solvable_ifrs ]);
          dimodels.append(model);

        ## loop over all DI terms (unless we're in the last major cycle, in which case only do the first one)
        di_solved = False;
        for i,opt in enumerate(self.gainopts):
          model = dimodels.pop();
          ## do we need to solve for this DI term? 
          ## On the first loop, optionally skip the first term
          ## On the last loop, we solve for the outermost term only
          dprintf(1,"%s: solvable %d from major loop %d (current %d)\n",opt.label,opt.solve,opt.nmajor_start,nmajor);
          if opt.solve and (nmajor >= opt.nmajor_start) and not (last_loop and di_solved):
            # recompute noise, unless this is the outer term, since it will have been rescaled
            if i:
              newnoise,weight = self.compute_noise(data,bitflags);
            ## iterate to convergence
            flagged |= self.run_gain_solution(opt,model,data,weight,bitflags,flag_null_gains=True,looptype=looptype);
            di_solved = True;
          ## apply correction to data
          dprint(1,"applying %s-inverse to data"%opt.label);
          data = dict([ (pq,opt.solver.apply_inverse(data,pq,
              regularize=self.regularization_factor if self.regularize_intermediate or last_loop else 0))
              for pq in solvable_ifrs ]);
          dprint(1,"done");
          ## check for NANs in the data
          self.check_finiteness(data,"corrected data",bitflags);
          
        ## ok, now: data is corrected data B^{-1}.G^{-1}.D.G^{-H}.B^{-H}
        ## model is simply M
          
        ## -------------------- now iterate over all diffgains
        if num_diffgains and not last_loop:
          ## if something was flagged in the DI solutions above, apply flags to dgmodels
          ## by resetting appropriate values to 0. data and model should already contain nulls at flagged slots.
          #if flagged:
            #for pq,fm in bitflags.iteritems():
              #for visset in [data,model,model0]+dgmodels:
                #for m in visset.get(pq,()):
                  #if not is_null(m):
                    #m[fm] = 0;
          # subtract model from corrected data to make residuals
          for pq in solvable_ifrs:
            for n,(dd,mm) in enumerate(zip(data[pq],model[pq])):
              dd -= mm;
          # recompute noise on residuals -- this will have been rescaled by the DI solutions      
          noise,resweight = self.compute_noise(data,bitflags);
          # reset current model to a _copy_ of model0 -- each DG term will be added to it (in place) at the end of
          # each DG loop iteration
          for pq in solvable_ifrs:
            model[pq] = matrix_copy(model0[pq]);
          # now loop over all diffgains and iterate each set once
          for idg,dgopt in enumerate(self.dgopts):
            # we want to solve for dE1 minimizing D=G.(M0+dE1.M1.dE1^H+dE2.M2.dE2^H+...).G^H
            # which is the same as G^{-1}.D.G^{-H} = M0 + dE1.M1.dE1^H + dE2.M2.dE2^H + ...
            # which is the same as D_corr - (M0 + dE1.M1.dE1^H + dE2.M2.dE2^H + ...) + dE1.M1.dE1^H = dE1.M1.dE1^H
            # which is the same as D_corr - full_model_old                           + model1_old   = model1
            # at this point, data is D_corr - full_model_old
            # so, add model1_old to data, and fit model1 to it
            # dgmodel_corr always contains the modelN_old values
            report_prec = False
            for pq in solvable_ifrs:
              if not report_prec:
                dprintf(2,"%s %s data type of dgmodel is %s\n",dgopt.label,pq,dgmodel[idg][pq][0].dtype)
              corr = dgopt.solver.apply(dgmodel[idg],pq);
              ##dgm: corr = dgmodel_corr[idg][pq];
              if not report_prec:
                dprintf(2,"%s %s data type of corrected is %s\n",dgopt.label,pq,corr[0].dtype)
                report_prec = True
              for (c,dd) in zip(corr,data[pq]):
                dd += c;
            if ( idg == self.dump_diffgain or self.dump_diffgain == -1 ) and \
              ( domain_id == self.dump_domain or self.dump_domain == -1 ):
              dump_data_model(dgmodel[idg],data,solvable_ifrs,"dump_E%d-%d.txt"%(idg,nmajor));
            # iterate this diffgain solution
            self.check_finiteness(data,"data after DG%d added in"%idg,bitflags);
            flagged = self.run_gain_solution(dgopt,dgmodel[idg],data,weight,bitflags,flag_null_gains=False,looptype=looptype);
            # now, add to model1 to model, and subtract back from data if needed
            for pq in solvable_ifrs:
              corr = dgopt.solver.apply(dgmodel[idg],pq);
              #dgm: corr = dgmodel_corr[idg][pq] = dgopt.solver.apply(dgmodel[idg],pq,cache=True);
              for i,(c,mm,dd) in enumerate(zip(corr,model[pq],data[pq])):
                if is_null(mm):
                  model[i] = c;
                else:
                  mm += c;
                if idg<num_diffgains-1:
                  dd -= c;
          # at end of loop over diffgains:
          # model already contains an up-to-date model with dEs applied
          
      # end major loop:
      #         data contains the corrected data (corrected by all DI terms)
      #         model contains an up-to-date model with dEs applied
      # However, this is only the case for solvable_ifrs.
      # We still need to update both for non-solvable IFRs
      dprint(1,"updating non-solvable IFRs");
      missing_ifrs = set(valid_ifrs) - set(solvable_ifrs);
      for pq in missing_ifrs:
        # fix up model
        mm = model[pq] = model0[pq];
        for idg,dg in enumerate(self.dgopts):
          #dgm: corr = dgmodel_corr[idg][pq] = dg.solver.apply(dgmodel[idg],pq,cache=True);
          corr = dg.solver.apply(dgmodel[idg],pq);
          for i,c in enumerate(corr):
            mm[i] += c;
      # apply corrections to missing baselines in data
      data1 = initdata;
      for opt in self.gainopts:
        data1 = dict([ (pq,opt.solver.apply_inverse(data1,pq,
                      regularize=self.regularization_factor))
                    for pq in missing_ifrs ]);
      data.update(data1);
      dprint(1,"saving solutions");        
      for opt in self.gainopts+self.dgopts:
        opt.save_values();
      GainOpts.flush_tables();
    # endif not skip_solve
    else:
      # no solve -- simply apply corrections to data
      data = initdata;
      for opt in self.gainopts:
        data = dict([ (pq,opt.solver.apply_inverse(data,pq,
                      regularize=self.regularization_factor))
                      for pq in data.iterkeys() ]);

#    dprint(0,"***DEBUG*** data",pq00,data[pq00][0][DEBUG_SLICE])
#    dprint(0,"***DEBUG*** model",pq00,model[pq00][0][DEBUG_SLICE])
    
    # corrdata will contain the corrected data (with all DI terms applied)
    # data will contain the original data
    # model already contains an up-to-date model with dEs applied
    corrdata,data = data,initdata;
    
    # remember init value for next tile
    if self.init_from_previous:
      for opt in self.gainopts+self.dgopts:
        opt.update_initval();
    
    dprint(1,"checking flagging");  
    # check for excessive flagging
    nfl = ndata = 0;
    for pq in data.iterkeys():
      ndata += self._datasize;
      fmask = bitflags.get(pq);
      # count flagged but unpadded points
      if fmask is not None:
        nfl += (((fmask&~FPRIOR)!=0)&self._expansion_mask).sum();
    flagpc = (nfl - nflagged_due_to_insufficient*len(data))*100./ndata;
    if flagpc > self.critical_flag_threshold:
      dprint(0,"***CRITICAL*** %.2f%% (%d/%d) data points were flagged in the stefcal process. Can not take. Stopping!"%(flagpc,nfl,ndata));
      raise RuntimeError,"Too many data points (%.2f%%) flagged. Check your stefcal settings?"%flagpc;
    dprint(1,"%.2f%% (%d/%d) data points were flagged in the stefcal process. Can take."%(flagpc,nfl,ndata));
    
    # go back to full-resolution data, if we were downsampling
    if downsample_subtiling:
      self._expanded_datashape = expanded_datashape = orig_sampled_expanded_datashape;
      self._datashape = datashape = orig_sampled_datashape;
      if not self.downsample_output:
        gain._reset();
        data = orig_sampled_data;
        model0 = model = orig_sampled_model0;
        dgmodel = orig_sampled_dgmodel;
        self._expansion_mask = orig_expansion_mask;
        self._datasize = reduce(operator.mul,datashape);
        for pq,bf in bitflags.items():
          fl1 = bitflags[pq] = numpy.zeros(expanded_datashape,int);
          fl1[...] = downsampler.expand_subshape(fl);
          #if num_diffgains:
            #model = {};
            #for pq,mod0 in model0.iteritems():
              #mm = model[pq] = mod0;
              #for idg,dg in enumerate(diffgains):
                #for i,c in enumerate(dg.apply(dgmodel[idg],pq,tiler=diffgaintiler)):
                  #mm[i] += c;

    # visualize gains
    for opt in self.gainopts+self.dgopts:
      if opt.visualize:
        print "saving global gain",opt.label;
        global_gains[opt.label] = [ opt.solver.get_2x2_gains(expanded_datashape,expanded_dataslice,
                                    tiler=opt.vis_tiler) ];
                                    
    # update IFR gain solutions, if asked to
    corrupt_model = None;
    if self.solve_ifr_gains:
      dprint(1,"generating corrupt model for IFR gain update");
      # make corrupted model
      corrupt_model = model;
      for opt in self.gainopts[-1::-1]:
        corrupt_model = dict([ (pq,opt.solver.apply(corrupt_model,pq,tiler=opt.tiler)) for pq in valid_ifrs ]);
#      dprint(1,"checking finiteness");
#      self.check_finiteness(corrupt_model,"corrupt model",bitflags);
#      self.check_finiteness(data,"data",bitflags);
#      cPickle.dump((corrupt_model,model,data,bitflags,[ (o.solver._gmat,o.solver._ginv) for o in self.gainopts]),file("dump.cp","wb"),2);
      dprint(1,"done, updating IFR gains");
      # now update IFR solutions
      for pq in self._ifrs:
        dd = data.get(pq);
        if dd is None:
          continue;
        mm = corrupt_model.get(pq);
        flagmask = bitflags.get(pq);
        for num,(d,m) in enumerate(zip(dd,mm)):
          # skip off-diagonal elements when in diagonal mode
          if self.diag_ifr_gains and num in (1,2):
            continue;
          # zero flagged data/model (ok to destroy now!)
          if not is_null(flagmask):
            d[flagmask] = 0;
            m[flagmask] = 0;
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
          dprint(4,"ifr gains: summed m.d* and d.d* for",pq,num);
          sri = self.ig_sum_reim[pq][num] = self.ig_sum_reim[pq][num] + mdh;
          ssq = self.ig_sum_sq[pq][num]   = self.ig_sum_sq[pq][num] + ddh;
          dprint(4,"updated sums for",pq,num);
          if numpy.isscalar(ssq):
            if ssq != 0:
              self.ifr_gain_update[pq][num] = sri/ssq;
          else:
            if (ssq!=0).any():
              self.ifr_gain_update[pq][num] = sri/ssq;
              self.ifr_gain_update[pq][num][ssq==0] = 1; 
          dprint(4,"updated ifr gain for",pq,num);
#          if num == 0 and pq[0] == '0':
#           printm[0,0],d[0,0],dh[0,0]
#            print pq,(m*dh).sum(),(d*dh).sum(),sri/ssq;
      dprint(1,"IFR gains updated");

    # work out result -- residual or corrected visibilities, depending on our state
    variance = {};
    nvells = 0;
    dprint(1,"computing result");
    for pq in self._ifrs:
      dd = corrdata.get(pq);
      mm = model.get(pq);
      if mm is None:
        for i in range(4):
          vs = datares.vellsets[nvells];
          if getattr(vs,'value',None) is not None:
            fl = getattr(vs,'flags',None);
            if fl is None:
              fl = vs.flags = meq.flags(datashape);
            else:
              fl = vs.flags = fl.copy()
            fl[...] |= self.output_flag_bit;
          nvells += 1;
        continue;
      else:
        if self.residuals:
          out = [ d-m for d,m in zip(dd,mm) ];
#          out = mm  ### write model!
#          if pq == pq00:
#            dprint(0,"***DEBUG*** residuals:",pq00,out[0][DEBUG_SLICE])
        else:
          out = dd;
          # subtract dE'd sources, if so specified
          if self.subtract_dgsrc:
            for idg,dg in enumerate(self.dgopts):
              corr = dg.solver.apply(dgmodel[idg],pq);
              for d,m in zip(out,corr):
                d -= m;
            #dgm: for idg,dgcorr in enumerate(dgmodel_corr):
            #dgm:   for d,m in zip(out,dgcorr[pq]):
            #dgm:     d -= m;
        # get flagmask, clear prior flags
        flagmask = bitflags.get(pq);
        # clear prior flags
        if flagmask is not None:
          flagmask &= ~FPRIOR;
          if self.downsample_output and downsampler and not numpy.isscalar(flagmask):
            flagmask = downsampler.expand_subshape(flagmask);
          if not flagmask.any():
            flagmask = None;
        for n,x in enumerate(out):
          if self.downsample_output and downsampler and not numpy.isscalar(x):
            x = downsampler.expand_subshape(x);
          vs = datares.vellsets[nvells];
          val = getattr(vs,'value',None);
          if val is not None:
            vs.value = val = val.copy()
            try:
              val[...] = x[expanded_dataslice] if expanded_dataslice \
                and not is_null(x) else x;
            except:
              print x,getattr(x,'shape',None);
          if not is_null(flagmask) and self.output_flag_bit:
            newflags = (flagmask!=0);
            nnew = newflags.sum();
            if nnew:
              counts = [];
              for bit,label in (FINSUFF,"n/d"),(FNOCONV,"n/c"),(FCHISQ,"chi2"),(FSOLOOB,"oob"):
                nf = ((flagmask&bit)!=0).sum();
                if nf:
                  counts.append("%s: %d"%(label,nf));
              dprint(3,"generated %d new flags in baseline %s-%s (%s)"%(nnew,pq[0],pq[1]," ".join(counts)));
              fl = getattr(vs,'flags',None);
              if fl is None:
                fl = vs.flags = meq.flags(datashape);
              else:
                fl = vs.flags = fl.copy()
              fl[newflags if expanded_dataslice is None else newflags[expanded_dataslice]] |=  self.output_flag_bit;
          # compute stats
          nvells += 1;
    dprint(1,"computing result: done");

    # if last domain, then write ifr gains to file
    if self.solve_ifr_gains and time1 >= numtime:
      dprint(1,"saving IFR gains");
      # get freq slicing (to go back from expanded shape to true data shape)
      if self.per_chan_ifr_gains and self._expanded_dataslice:
        slc = self._expanded_dataslice[1];
      else:
        slc = None;
#      dprint(0,"slice is",slc);
      # apply updates
      for pq in self._ifrs:
        self.ifr_gain[pq] = [ g*(g1 if numpy.isscalar(g1) or not slc else g1[slc]) 
            for g,g1 in zip(self.ifr_gain.get(pq,[1,1,1,1]),self.ifr_gain_update[pq]) ];
#        dprint(0,pq,"shape is",getattr(self.ifr_gain[pq][0],'shape',[1]));
      dprint(2,"IFR gain solutions update: ",", ".join(
            ["%s%s:%s%s %s"%(p,self.corr_names[i],q,self.corr_names[j],
            self.ifr_gain_update[(p,q),i,j])
            for (p,q),i,j in pqij_data[0:3]]));
      dprint(2,"IFR gain solutions: ",", ".join(
            ["%s%s:%s%s %s"%(p,self.corr_names[i],q,self.corr_names[j],
            self.ifr_gain[(p,q),i,j])
            for (p,q),i,j in pqij_data[0:3]]));
      # save
      if self.save_ifr_gains:
        try:
          cPickle.dump(self.ifr_gain,file(self.ifr_gain_table,'w'),2);
          dprint(1,"saved %d ifr gains to %s"%(len(self.ifr_gain),self.ifr_gain_table));
        except:
          traceback.print_exc();
          dprint(0,"error saving ifr gains to",self.ifr_gain_table);

    dt = time.time()-timestamp0;
    m,s = divmod(dt,60);
    dprint(0,"%s elapsed time %dm%0.2fs"%(
              request.request_id,m,s));

    return datares;

  def compute_noise (self,data,bitflags):
    """Computes delta-std and weights of data""";
    noise = {};
    weight = {};
    sumweight = 0;
    nweight = 0;
    # make null-flag array
    dfshape = list(self._expanded_datashape);
    dfshape[0] -= 1;
    alltrue = numpy.ones(dfshape,bool);
    allfalse = numpy.zeros(dfshape,bool);
    # loop over data
    for pq,dd in data.iteritems():
      flag = bitflags.get(pq);
      if not is_null(flag):
        flag = (flag!=0);
        dflag = flag[1:,...]|flag[:-1,...];
        dvalid = ~dflag;
      else:
        dflag,dvalid = allfalse,alltrue;
      # number of valid slots (per channel, or in total)
      num_valid = dvalid.sum(0) if self.noise_per_chan else dvalid.sum();
      vv2 = [];
      for d in dd:
        if is_null(d):
          vv2.append(0);
        else:
          # take forward difference, null at flagged points
          delta = d[1:,...]-d[:-1,...];
          delta[dflag] = 0;
          # take squared real and imaginary parts of this
          # sum them, since taking the difference reduces the noise by sqrt(2); so the
          # squared-mean-diff is a factor of 2 higher
          if self.noise_per_chan:
            v2 = (numpy.square(delta.real).sum(0) + numpy.square(delta.imag).sum(0))/num_valid;
            v2[num_valid==0] = 0;
            vv2.append(v2[numpy.newaxis,:]);
          else:
            v2 = (numpy.square(delta.real).sum() + numpy.square(delta.imag).sum())/num_valid if num_valid else 0;
            vv2.append(v2);
      # convert to weight
      # if XY/YX is well-defined, use it, else use the XX/YY estimates
      if self.use_polarizations_for_noise and not is_null(vv2[1]) and not is_null(vv2[2]):
        n = noise[pq] = numpy.sqrt((vv2[1]+vv2[2])/2);
        w = weight[pq] = 1/n;
        w[n==0] = 0; 
      elif not is_null(vv2[0]) and not is_null(vv2[3]):
        n = noise[pq] = numpy.sqrt((vv2[0]+vv2[3])/2);
        w = weight[pq] = 1/n;
        if numpy.isscalar(n):
          if n == 0:
            w = 0;
        else:
          w[n==0] = 0; 
      else:
        n,w = 0,0;
    ## normalize weights ## NB why? what was I thinking?
    #if nweight:
      #meanweight = sumweight/nweight;
      #for pq,w in weight.iteritems():
        #weight[pq] = w/meanweight;
    
    if _verbosity.verbose>3:
      dprint(4,"noise estimates by baseline:");
      npq = sorted([ (n,pq) for pq,n in noise.iteritems() ],cmp=lambda x,y:cmp(numpy.mean(x[0]),numpy.mean(y[0])));
      for n,pq in npq:
        dprint(4,"  %s-%s"%pq," ".join(["%.2g"%float(x) for x in n.ravel()[:20]]));
        
    return noise,weight;

  def compute_chisq (self,model,data,gain,weight=None,bitflags={}):
    # per-slot normalized and unnormalized chisq
    chisq0 = numpy.zeros(self._expanded_datashape);
    chisq1 = numpy.zeros(self._expanded_datashape);
    # nterms: per-slot number of terms in chi-sq sum
    nterms = numpy.zeros(self._expanded_datashape,int);
    # loop over all IFRS
    for pq in self._solvable_ifrs:
      if pq not in data:
        continue;
      res = gain.residual(model,data,pq);
      w = weight.get(pq,1)**2 if weight else 1;
#      fmask = reduce(operator.or_,[(d==0)&(m==0) for (d,m) in zip(data[pq],model[pq])]);
      fmask = bitflags.get(pq);
      fmask = fmask is not None and (fmask!=0);
      # n0 is the nominal number of t/f slots for which we expect to have a residual
      n0 = self._datasize - ((fmask&self._expansion_mask).sum() if not is_null(fmask) else 0);
      for ir,r in enumerate(res):
        if not is_null(r):
          # fin is a mask of finite residuals, in unflagged slots
          # in principle all unflagged residuals ought to be finite, but I'm covering
          # my ass here in case of some pathologies/bugs
          fin = numpy.isfinite(r);
          if not is_null(fmask):
            fin &= ~fmask;
          # n is the number of slots for which we have a finite, unflagged residual
          n = fin.sum();
          if n < n0:
            dprintf(4,"%s element %d: %d/%d slots are unexpectedly INF/NAN, omitting from chisq sum\n",pq,ir,n0-n,n0);
          # add residual to chisq sum
          if n:
            rsq = (r*numpy.conj(r)).real;
            chisq0[fin] += (rsq*w)[fin];
            chisq1[fin] += (rsq)[fin];
            nterms[fin] += 2;  # each slot contributes two terms (real and imag)
    # ok chisq0 and chisq1 contain the per-slot chi-squares. Take their mean
    tot_terms = nterms.sum();
    if tot_terms:
      chisq0sum = float(chisq0.sum())/tot_terms;
      chisq1sum = float(chisq1.sum())/tot_terms;
      mask = nterms>0;
      norm = nterms;
      chisq0[mask] /= norm[mask];
      chisq1[mask] /= norm[mask];
    else:
      chisq0sum = chisq1sum = 0;
    return chisq0sum,chisq1sum,chisq0,chisq1;

  def check_finiteness (self,data,label,bitflags,complete=False):
    # only do this in high verbosity mode
    if _verbosity.verbose<4:
      return;
    # check for NANs in the data
    for pq,dd in data.iteritems():
      fmask = bitflags.get(pq);
      fmask = fmask is not None and (fmask!=0);
      for i,d in enumerate(dd):
        if type(d) is numpy.ndarray:
          fin = numpy.isfinite(d);
          if not is_null(fmask):
            fin |= fmask;
          if not fin.all():
            dprintf(0,"%s %s element %d: %d slots are INF/NAN\n",label,pq,i,d.size-fin.sum());
            if not complete:
              break;
              
  def add_flags (self,bitflags,pq,fmask):
    if pq in bitflags:
      bitflags[pq] |= fmask;
    else:
      bitflags[pq] = fmask;

  def run_gain_solution (self,gopt,model,data,weight,bitflags,flag_null_gains=False,looptype=0):
    """Runs a single gain solution loop to completion"""
    flagged = False;
    gain_dchi = [];
    gain_maxdiffs = [];
    # initial chi-sq, and fallback chisq for divergence
    gopt.solver._reset();
    init_chisq,init_chisq_unnorm,init_chisq_arr,init_chisq_unnorm_arr = \
      self.compute_chisq(model,data,gopt.solver,weight=weight,bitflags=bitflags);
    self._set_ds_array('$init_chisq',init_chisq_arr);
    lowest_chisq = (init_chisq,gopt.solver.get_values());
    lowest_chisq_iter = 0;
    num_diverged = 0;
    dprint(2,"solving for %s, initial chisq is %.12g"%(gopt.label,init_chisq));
    t0 = time.time();
    chisq0 = init_chisq;
    # iterate
    for niter in range(gopt.max_iter):
      # iterate over normal gains
      # bounds-flagging is enabled after iteration 3
      converged,maxdiff,deltas,nflag = gopt.solver.iterate(model,data,bitflags,
                                          niter=niter,weight=weight if gopt.weigh else None,
                                          bounds=gopt.bounds if niter>2 else None);
      dprint(3,"iter %d: %.2f%% (%d/%d) conv, %d gfs, max update %g"%(
          niter+1,gopt.solver.num_converged*100./gopt.solver.real_slots,gopt.solver.num_converged,gopt.solver.real_slots,nflag,float(gopt.solver.delta_max)));
      gain_maxdiffs.append(float(maxdiff));
      if gopt.visualize > 1:
        global_gains.setdefault(gopt.label,[]).append(gopt.solver.get_2x2_gains(expanded_datashape,expanded_dataslice));
      ## compute chisq if converged, or stopping, or using chi-sq convergence (delta!=0)
      delta = gopt.delta;
      if looptype:
        delta1 = getattr(gopt,"delta_loop%d"%looptype);
        if delta1 == "same":
          delta1 = gopt.delta;
      if delta != 0 or gopt.max_diverge or converged or niter == gopt.max_iter-1:
        chisq,chisq_unnorm,chisq_arr,chisq_unnorm_arr = self.compute_chisq(model,data,gopt.solver,weight=weight,bitflags=bitflags);
      if delta != 0:
        dchi = (chisq0-chisq)/chisq;
        gain_dchi.append(dchi);
        dprint(3,"iter %d: ||%s||=%g max update is %g chisq is %.8g (%.8g) diff %.3g, %d/%d conv"%(
          niter+1,gopt.label,float(gopt.solver.gainnorm),float(gopt.solver.delta_max),
          chisq,chisq_unnorm,dchi,gopt.solver.num_converged,gopt.solver.real_slots));
        if niter:
          # check chisq convergence
          dchi_arr = (chisq0_arr-chisq_arr)/chisq0_arr;
          diverging_mask = (dchi_arr<0);       # mask of diverging chi-squares
          diverging_sig_mask = (dchi_arr<-1e-2*chisq0_arr);       # mask of significantly diverging chi-squares
          nd = diverging_mask.sum();
          nds = diverging_sig_mask.sum();
          converging_mask = ((dchi_arr<delta) & ~diverging_mask);
          # make list of (num_converged,thr) pairs
          conv_stats = [(converging_mask.sum(),delta)] + [ 
              (((dchi_arr<d) & ~diverging_mask).sum(),d) for d in delta*10,delta*100,delta*1000,delta*10000 ];
          dprint(3,"iter %d: chisq converged %s; diverging %d; significantly %d"%(niter+1,
                  ", ".join(["%d to %g"%x for x in conv_stats ]),nd,nds));
        chisq0 = chisq;
        chisq0_arr = chisq_arr;
      # check solutions convergence
      if converged:
        dprint(1,"%s converged at chisq %.12g (last gain update %g) after %d iterations and %.2fs"%(
                gopt.label,chisq,float(gopt.solver.delta_max),niter+1,time.time()-t0));
        break;
      # else check for chi-sq convergence
      elif delta != 0:
        # if we got here (avoiding the break above
        # if chi-sq decreased, remember this
        if dchi >= 0:
          if chisq < lowest_chisq[0]:
            lowest_chisq = (chisq,gopt.solver.get_values());
            lowest_chisq_iter = niter+1;
          num_diverged = 0;
          # and check for chisq-convergence
          if niter > 1 and dchi < delta:
            dprint(1,"%s chisq converged at %.12g (last gain update %g) after %d iterations and %.2fs"%(
                  gopt.label,chisq,float(gopt.solver.delta_max),niter+1,time.time()-t0));
            break;
        # else chisq is increasing, check if we need to stop
        else:
          self._set_ds_array('$diverging_chisq',chisq_arr);
          if chisq > init_chisq*10:
            dprint(1,"%s: mad jump in chisq, aborting"%gopt.label);
            break;
          num_diverged += 1;
          if num_diverged >= gopt.max_diverge:
            dprint(1,"%s chisq diverging, stopping at %.12g after %d iterations and %.2fs"%(
                      gopt.label,chisq,niter+1,time.time()-t0));
            break;
    else:
      dprint(1,"%s max iterations (%d) reached at chisq %.12g (last gain update %g) after %.2fs"%(
              gopt.label,gopt.max_iter,chisq,float(gopt.solver.delta_max),time.time()-t0));
    # check if we have a lower chisq to roll back to
    rolled_back = False;
    if self.chisq_rollback:
      if chisq > lowest_chisq[0]:
        if not lowest_chisq_iter:
          dprint(1,"GONE NOWHERE FAST! Final chisq %g is higher than initial chisq %g. Rolling back!"%(chisq,lowest_chisq[0]));
        else:
          dprint(2,"Final chi-sq %g higher than minimum %g achieved at iteration %d, rolling back."%(
            chisq,lowest_chisq[0],lowest_chisq_iter));
        if chisq > lowest_chisq[0]*1.01 and gopt.flag_chisq:
          dprint(2,"Recomputing rolled-back chisq for flagging purposes");
          chisq,chisq_unnorm,chisq_arr,chisq_unnorm_arr = self.compute_chisq(model,data,gopt.solver,weight=weight,bitflags=bitflags);
        rolled_back = True;
        self._set_ds_array('$high_discarded_chisq',chisq_arr);
        chisq,gainvals = lowest_chisq;
        gopt.solver.set_values(gainvals);
    dprint(2,"  delta-chisq were"," ".join(["%.4g"%x for x in gain_dchi]));
    dprint(2,"  convergence criteria were"," ".join(["%.2g"%x for x in gain_maxdiffs]));
    #
    # flag non-converged solutions, unless rollback was in effect
    #
    if not rolled_back and gopt.flag_non_converged:
      fnonconv = (~gopt.solver.get_converged_mask())&self._expansion_mask;
      nfl = fnonconv.sum();
      if nfl:
        dprint(2,"flagging %f%% of data (%d/%d slots) due to non-convergence"%(nfl*100./self._datasize,nfl,self._datasize));
        fmask = fnonconv*FNOCONV;
        for pq in data.iterkeys():
          self.add_flags(bitflags,pq,fmask);
        flagged = True;
        chisq_arr[fnonconv] = 0;
    else:
      fnonconv = False;
    #
    # implement gain-clip flagging
    #
    if flag_null_gains or gopt.bounds:
      lower,upper = gopt.bounds or (None,None);
      gainflags = gopt.solver.get_gainflags(flagnull=flag_null_gains,lower=lower,upper=upper);
      if gainflags:
        dprint(2,"flagged gains per antenna:",", ".join(["%s %.2f%%"%(p,gf.sum()*100./gopt.solver.real_slots) for p,gf in sorted(gainflags.iteritems())]));
        for pq in data.iterkeys():
          fmask = gainflags.get(pq[0],0);
          fmask |= gainflags.get(pq[1],0);
          if not is_null(fmask):
            chisq_arr[fmask] = 0;
            self.add_flags(bitflags,pq,(fmask)*FSOLOOB);
            flagged = True;
    #
    # implement chisq-based flagging
    #
    if gopt.visualize:
      self._set_ds_array('$final_chisq',chisq_arr);
    if gopt.flag_chisq:
      if False:  # flag on histogram and mean
        # make histogram of chisq values
        chisq_nonzero = (chisq_arr!=0);
        # min bin must NOT be zero, we don't want to be picking up the flagged values
        chisq_histbins = [ chisq*10**x for x in [-99]+list(numpy.arange(-5,5,.5)) ];
        chisq_hist,dum = numpy.histogram(chisq_arr,bins=chisq_histbins);
        if _verbosity.verbose>2:
          for i,ch in enumerate(chisq_hist):
            dprint(3,"  chisq histogram: %d values in [%g,%g)"%(ch,chisq_histbins[i],chisq_histbins[i+1]));
        # find max bin, take this to be the nominal chi-sq value
        imaxbin = numpy.argmax(chisq_hist);
        dprint(3,"max chisq bin is %g"%chisq_histbins[imaxbin]);
        b0 = chisq_histbins[max(imaxbin-1,0)]
        b1 = chisq_histbins[min(imaxbin+1,len(chisq_histbins)-1)]
        mm = chisq_arr[(chisq_arr>=b0)&(chisq_arr<=b1)].mean();
        dprint(3,"M (mean chisq value in bin range %g:%g) = %g"%(b0,b1,mm));
        cPickle.dump((chisq_hist,chisq_histbins),file("stefcal.chisq.dump","w"),2);
      else:
        chisq_masked = numpy.ma.masked_array(chisq_arr,chisq_arr==0,fill_value=0);
        mode = getattr(gopt,"flag_chisq_loop%d"%looptype);
        mm = numpy.ma.median(chisq_masked,axis={FREQMED:0,TIMEMED:1}.get(mode,None));
        if numpy.ma.is_masked(mm):
          mm = mm.filled(fill_value=1e+999);
        dprint(3,"M (median chisq value):",mm);
      # get flagmask based on theshold  
      chisq_flagmask = (chisq_arr > mm*gopt.flag_chisq_threshold);
      nfl,nsl = chisq_flagmask.sum(),(chisq_arr!=0).sum();
      dprint(2,"flagging %f%% of data (%d/%d slots) due to chi-sq > M*%.1f"%(
          nfl*100./nsl,nfl,nsl,gopt.flag_chisq_threshold));
      # make trial stats for other thresholds
      for f in 2,3,5,10,20:
        nfl1 = ((chisq_arr>mm*f)).sum(); 
        dprint(2,"  for threshold M*%.1f this would have been %f%%"%(f,nfl1*100./nsl));
      # apply flags
      if chisq_flagmask.any():
        flagged = True;
        fmask = FCHISQ*chisq_flagmask;
        for pq in data.iterkeys():
          self.add_flags(bitflags,pq,fmask);
      # dump chisq in output record
      if gopt.visualize:
        chisq_arr[chisq_flagmask] = 0; 
        self._set_ds_array('$final_chisq_flagged',chisq_arr);
#    self.set_state('$final_chisq_histbins',chisq_histbins); 
#    self.set_state('$final_chisq_hist',chisq_hist); 
    return flagged;
    
  def _set_ds_array (self,field,array):
    vv = meq.vells(array.shape);
    vv[...] = array;
    self.set_state(field,array); 
