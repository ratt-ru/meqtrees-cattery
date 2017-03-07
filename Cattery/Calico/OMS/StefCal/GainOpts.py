#from memory_profiler import profile

import sys
import Kittens.utils
import os.path
import os
import cPickle
import numpy
import traceback

MODE_SOLVE_SAVE = "solve-save";
MODE_SOLVE_NOSAVE = "solve-nosave"
MODE_SOLVE_APPLY = "apply"

TIMEMED = "time"; 
FREQMED = "freq";
TOTMED  = "overall";

_verbosity = Kittens.utils.verbosity(name="gainopts");
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



class GainOpts (object):
  """Encapsulates a set of options for a gains object"""
  def __init__ (self,desc,name,label,tdl_basespace=None,node=None,mystate=None,pre_opts=[],post_opts=[]):
    self.desc,self.name,self.label = desc,name,label;
    ### init options on TDL side of things
    if tdl_basespace:
      from Timba.TDL import *
      self.tdloption_namespace = "%s_%s"%(tdl_basespace,name.lower());
      timeint = TDLOption("timeint","Solution interval, time axis (0 for full axis)",[0,1],more=int,default=1,namespace=self);
      freqint = TDLOption("freqint","Solution interval, freq axis (0 for full axis)",[0,1],more=int,default=1,namespace=self);
      timesmooth = TDLOption("timesmooth","Smoothing kernel, time axis",[0],more=int,namespace=self);
      freqsmooth = TDLOption("freqsmooth","Smoothing kernel, freq axis",[0],more=int,namespace=self);
      timesmooth.when_changed(lambda x,opt=timeint:opt.show(not x));
      freqsmooth.when_changed(lambda x,opt=freqint:opt.show(not x));
      labelopt = TDLOption("label","Jones matrix label",[label],more=str,default=label,validator=bool,namespace=self);
      meddict = { TIMEMED:"per timeslot",FREQMED:"per channel",TOTMED:"overall" };
      menuopts = pre_opts + [
          labelopt,
          TDLOption("use_float","Use single precision",False,namespace=self),
          timeint,freqint,timesmooth,freqsmooth,
          TDLOption("flag_nonconv","Flag non-converging bins",False,namespace=self),
          TDLMenu("Flag using chi-square",
              TDLOption('flag_chisq_threshold',"Threshold, in N*median",[0,3,5,10],more=int,default=5,namespace=self),
              TDLOption('flag_chisq_loop0',"Thresholding type, first major cycle",meddict,default=TOTMED,namespace=self),
              TDLOption('flag_chisq_loop1',"Thresholding type, intermediate cycles",meddict,default=TOTMED,namespace=self),
              TDLOption('flag_chisq_loop2',"Thresholding type, final cycle",meddict,default=TOTMED,namespace=self),
            toggle='flag_chisq',namespace=self),
          TDLMenu("Flag using solution amplitude clipping",
              TDLOption("flag_ampl_low","Lower threshold (0 disables)",[0,.5],more=float,default=0,namespace=self),
              TDLOption("flag_ampl_high","Upper threshold (0 disables)",[0,1.5],more=float,default=0,namespace=self),
            toggle='flag_ampl',namespace=self),
          TDLOption("implementation","Jones matrix type",["GainDiag","Gain2x2","Gain2x2a","GainDiagCommon","GainDiagPhase" ] ,namespace=self),
          TDLOption("mode","Solution mode",
            {MODE_SOLVE_SAVE:"solve and save",MODE_SOLVE_NOSAVE:"solve, do not save",MODE_SOLVE_APPLY:"load and apply"},
            default=MODE_SOLVE_SAVE,namespace=self),
          TDLOption("visualize","Enable visualization",True,namespace=self),
          TDLMenu("Ninja options",
##              TDLOption("global_solution","One solution per all antennas",False,namespace=self),
              TDLOption("real_only","Constrain to real solutions only",False,namespace=self),
              TDLOption("nmajor_start","Start solving at major loop",[0,1],more=int,namespace=self),
              TDLOption("weigh","Apply noise-based weights",False,namespace=self),
              TDLOption("niter","Max iterations",[20,50,100],more=int,default=50,namespace=self),
              TDLOption("epsilon","Solution convergence criterion",[1e-4,1e-5,1e-6,1e-7,1e-8],more=float,default=1e-6,namespace=self),
              TDLOption("quota","Solution convergence quorum",[.95,.99,1],more=float,namespace=self),
              TDLOption("delta","Chi-sq convergence criterion (0 disables)",
                          [0,1e-3,1e-6],more=float,default=1e-6,namespace=self),
              TDLOption("delta_loop1","...for intermediate major loops",
                          ["same",0,1e-3,1e-6],more=float,default="same",namespace=self),
              TDLOption("delta_loop2","...for final major loop",
                          ["same",0,1e-3,1e-6],more=float,default="same",namespace=self),
              TDLOption("max_diverge","Iterations to allow chi-sq to diverge for",[0,1,5],more=int,default=1,namespace=self),
              TDLOption("omega","Averaging weight (omega)",[0.5,1.8],more=float,namespace=self),
              TDLOption("average","Averaging mode",[0,1,2],default=2,namespace=self),
              TDLOption("ff","Enable feed-forward averaging",True,namespace=self),
              TDLOption("table","Filename for solution table",["%s.cp"%name],more=str,namespace=self),
              TDLOption("intermediate_table","Filename for intermediate values table",[None,"intermediate-%s.cp"%name],more=str,namespace=self),
            )
        ] + post_opts;
      self._menuopt = TDLMenu("Use '%s' %s"%(label,desc),toggle='enabled',namespace=self,*menuopts)
      labelopt.when_changed(lambda label,opt=self._menuopt,desc=desc:opt.set_name("Use '%s' %s"%(label,desc)));
      self.tdl_options = [ self._menuopt ];
    # other stuff
    self.tiler = self.vis_tiler = None;
      
  def update_state (self,node,mystate=None,verbose=None,option_suffix=None):
    """Called from StefCal node to update state record""";
    if verbose is not None:
      _verbosity.set_verbose(verbose);
      _verbosity.enable_timestamps(True,modulo=6000);
    
    for option,default in [
            ('label','J'),
            ('enable',False),
            ('use_float',False),
            ('chisq_flag',5),
            ('epsilon',1e-5),            # when the update is ||G-G'||<epsilon, we are converged
            ('delta',1e-6),              # when chisq changes by less than delta, we are converged
            ('delta_loop1',"same"),              # when chisq changes by less than delta, we are converged
            ('delta_loop2',"same"),              # when chisq changes by less than delta, we are converged
            ('max_iter',50),             # max gain iter 
            ('max_diverge',10),          # max iterations to diverge allowed 
            ('convergence_quota',0.99),   # what percentage of parms should converge
            ('weigh',True),
            ('subtiling',[1,1]),
            ('smoothing',[]),
            ('omega',.5),
            ('average',2),
            ('feed_forward',False),
            ('solve',True),
            ('save',True),
            ('global',False),
            ('real_only',False),
            ('nmajor_start',0),
            ('flag_non_converged',False),
            ('flag_chisq',False),
            ('flag_chisq_threshold',5),
            ('flag_chisq_loop0',"overall"),
            ('flag_chisq_loop1',"overall"),
            ('flag_chisq_loop2',"overall"),
            ('visualize',True),
            ('bounds',[]),
            ('table','%s.cp'%self.name),
            ('intermediate_table',None),
            ('implementation','GainDiag'),
          ]:
      # for each OPTION, init node state field NAME_OPTION,
      optname = '_'.join([self.name,option_suffix,option] if option_suffix else [self.name,option]);
      # if mystate is given, call it to init the node field
      if mystate:
        mystate(optname,default);
      if hasattr(node,optname):
        value = getattr(node,optname);
        # ...and copy its value into our OPTION attribute
        setattr(self,option,value);
        # and print
        dprint(4,"%s: setting"%self.name,optname,"=",value);
    self.polarized = False;
    # get implementation class
    # given X.Y.Z, import module X.Y, and use symbol Z from that
    # given just X, use X.X
    if self.enable:
      path = self.implementation.split('.');
      modname,classname = '.'.join(['Calico.OMS.StefCal']+(path[:-1] or path[-1:])),path[-1];
      __import__(modname);
      module = sys.modules[modname];
      self.impl_class = getattr(module,classname);
      self.polarized = getattr(self.impl_class,'polarized',True);

  def set_stefcal_node_options (self,kw,visualize=False):
    """converts TDL options into node keyword options in the kw dict"""
    name = self.name;
    kw['%s_enable'%name]     = self.enabled;
    kw['%s_label'%name]      = self.label;
    kw['%s_use_float'%name]  = self.use_float;
    kw['%s_flag_non_converged'%name] = self.flag_nonconv;
    kw['%s_flag_chisq'%name] = self.flag_chisq;
    kw['%s_flag_chisq_threshold'%name] = self.flag_chisq_threshold;
    kw['%s_flag_chisq_loop0'%name] = self.flag_chisq_loop0;
    kw['%s_flag_chisq_loop1'%name] = self.flag_chisq_loop1;
    kw['%s_flag_chisq_loop2'%name] = self.flag_chisq_loop2;
    kw['%s_epsilon'%name]    = self.epsilon;
    kw['%s_delta'%name]      = self.delta;
    kw['%s_delta_loop1'%name]= self.delta_loop1;
    kw['%s_delta_loop2'%name]= self.delta_loop2;
    kw['%s_max_diverge'%name]= self.max_diverge;
    kw['%s_max_iter'%name]   = self.niter;
    kw['%s_convergence_quota'%name] = self.quota
    kw['%s_weigh'%name]      = self.weigh;
    kw['%s_omega'%name]      = self.omega;
    kw['%s_average'%name]    = self.average;
    kw['%s_feed_forward'%name] = self.ff;
    kw['%s_table'%name]      = self.table;
    kw['%s_intermediate_table'%name] = self.intermediate_table;
    kw['%s_solve'%name]      = (self.mode != MODE_SOLVE_APPLY);
    kw['%s_save'%name]       = (self.mode == MODE_SOLVE_SAVE);
    kw['%s_implementation'%name] = self.implementation;
    kw['%s_bounds'%name]     = self.flag_ampl and [self.flag_ampl_low,self.flag_ampl_high];
##   kw['%s_global_solution'%name] = self.global_solution;
    kw['%s_real_only'%name]  = self.real_only;
    kw['%s_nmajor_start'%name] = self.nmajor_start;
    kw['%s_visualize'%name]  = visualize and self.visualize;
    subtiling = [ self.timeint,self.freqint ];
    smoothing = [ self.timesmooth,self.freqsmooth ];
    for i in range(2):
      if smoothing[i]:
        subtiling[i] = 1;
    if all([not x for x in smoothing]):
      smoothing = [];
    kw['%s_subtiling'%name]  = subtiling;
    kw['%s_smoothing'%name]  = smoothing;
  
  _incoming_tables = {};
  
  def load_initval (self,default):
    """Loads initial values from table (if available)"""
    self.init_value = default;
    self.has_init_value = False;
    if not self.enable:
      return;
    if not os.path.exists(self.table):
      dprint(0,"not loading %s solutions: %s does not exist"%(self.label,self.table));
      return;
    try:
      struct = GainOpts._incoming_tables.get(self.table);
      if not struct:
        struct = GainOpts._incoming_tables[self.table] = cPickle.load(file(self.table));
      if not isinstance(struct,dict) or struct.get('version',0) < 2:
        dprint(0,"error loading %s solutions: %s format or version not known"%(self.label,self.table));
        return;
      gains = struct['gains'];
      if self.label not in gains:
        dprint(0,"no %s solutions found in %s (table contains: %s)"%(self.label,self.table,", ".join(sorted(gains.keys()))));
        return;
      gains = gains[self.label];
      if gains['implementation'] != self.implementation:
        dprint(0,"%s solutions in %s are for class %s, expected %s"%(self.label,self.table,gains['implementation'],self.implementation));
      self.init_value = gains['solutions'];
      self.has_init_value = True;
      dprint(1,"loaded %d %s solutions from %s"%(len(self.init_value),self.name,self.table));
    except:
      traceback.print_exc();
      dprint(0,"error loading %s solutions from"%self.label,self.table);

  _outgoing_tables = {};
        
  def save_values (self):
    if self.save:
      GainOpts._outgoing_tables.setdefault(self.table,{})[self.label] = \
        dict(solutions=self.solver.gain,implementation=self.implementation);

  def save_intermediate_values (self,niter):
    if self.intermediate_table:
      GainOpts._outgoing_tables.setdefault(self.intermediate_table,{})[self.label,niter] = \
        dict(solutions=self.solver.gain,implementation=self.implementation);

  @staticmethod 
  def flush_tables ():
    GainOpts._incoming_tables = {};
    for table,initval in GainOpts._outgoing_tables.iteritems():
      if initval:
        struct = dict(description="stefcal gain solutions table",version=2,gains=initval);
        try:
          cPickle.dump(struct,file(table,'w'),2);
          dprint(1,"saved %d gain set(s) to %s"%(len(initval),table));
        except:
          traceback.print_exc();
          dprint(0,"error saving gains to",table);
    GainOpts._outgoing_tables = {};

  @staticmethod
#  @profile
  def resolve_tilings (datashape,*opts):
    """Resolves a number of GainOpts into a common tiling"""
    lcm_tiling = [0]*len(datashape);
    
    for opt in opts:
      if opt.enable:
        if len(opt.subtiling) != len(datashape):
          raise ValueError,"%s tiling vector must have the same length as the data shape"%self.name;
        if min(opt.subtiling) < 0:
          raise ValueError,"invalid %s tiling %s"%(opt.name,opt.subtiling);
        opt.subtiling = opt.subtiling or [0]*len(datashape);
        # work out least-common-multiple subtile size for each axis on which a subtiling is defined.
        lcm_tiling = [ LCM(a,b) if a>0 and b>0 else max(a,b,0) for a,b in zip(lcm_tiling,opt.subtiling) ];
    
    # lcm_tiling along each axis now contains the LCM of each gain term's tiling, or 0 if all were 0
    # expand datashape as needed
    expanded_datashape = \
      tuple([ (nd/np+(1 if nd%np else 0))*np if np else nd for nd,np in zip(datashape,lcm_tiling) ]);
    dprint(1,"datashape",datashape,"expanded datashape is",expanded_datashape);
    
    # now, for any zeroes left in the tilings, replace with full solution interval,
    lcm_tiling = [ ls or ds for ls,ds in zip(lcm_tiling,expanded_datashape) ];
    for opt in opts:
      if opt.enable:
        opt.subtiling = [ gs or ds for gs,ds in zip(opt.subtiling,expanded_datashape) ];
        dprint(1,"%s tiling is"%opt.name,opt.subtiling,"smoothing is",opt.smoothing);
    dprint(1,"based on an LCM tiling of",lcm_tiling);

    return expanded_datashape;

  def init_solver (self,datashape,expanded_datashape,solvable_ifrs,downsample_subtiling):
    """Initializes gain solver object""";
    if not self.enable:
      return;
    dprintf(0,"stefcal %s solve=%d %s, using %d solvable inteferometers\n",
      self.label,self.solve,self.impl_class.__name__,
      len(solvable_ifrs));
    dprint(0,"  solution intervals:",self.subtiling,"smoothing kernel:",self.smoothing);
    if self.bounds:
      dprint(0,"  gains will be flagged on amplitudes outside of",self.bounds);
    if self.has_init_value:
      initval = self.init_value.values()[0][0];
      dprint(1,"  initial values loaded, first number is",initval.flat[0] if hasattr(initval,'flat') else initval);
    else:
      dprint(1,"  default initial value is",self.init_value);
    # init gain parms object
    self.solver = self.impl_class(datashape,expanded_datashape,
        self.subtiling,solvable_ifrs,opts=self,
        force_subtiling=bool(downsample_subtiling),
        init_value=self.init_value,
        verbose=_verbosity.verbose);
    dprint(1,"  subshape",self.solver.subshape,"tiled",self.solver.tiled_shape);

  def update_initval (self):
    self.init_value = self.solver.get_last_timeslot();
        