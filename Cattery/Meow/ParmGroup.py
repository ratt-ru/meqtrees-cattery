# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from Timba.TDL import *
from Meow import Context
from Meow import Bookmarks
from .SolverControl import SolverControl

import os.path
import os

# converts argument to a name
def namify (arg):
  if isinstance(arg,str):
    return arg;
  elif is_node(arg):
    return arg.name;
  else:
    raise TypeError("node name or node object expected, got '%s'"%str(type(arg)));

_all_parmgroups = [];

class Subgroup (object):
  def __init__ (self,name,nodes):
    self.name = name;
    self.nodes = ParmGroup._sort_members(nodes);
    self.solvable = True;

  def set_solvable (self,solvable):
    self.solvable = solvable;

class ParmGroup (object):
  class Controller (object):
    """A ParmGroup Controller implements a number of TDLOptions associated with the ParmGroup,
    and provides methods for changing the state of the parms"""
    def __init__ (self,pg,label,solvable=False,
                  **kw):
      self.pg = pg;
      self.label = label = "%s.%s"%(label,self.pg.label);
      self.tdloption_namespace = label.replace(":","_").replace(" ","_");
      # build up list of options
      self._option_list = [];
      # make individual solvability tdloptions for every parm, if so asked
      self._individual = pg._individual;
      sopts = [];  # these will need to be shown/hidden with the solvability control
      if self._individual:
        self._parm_opts = {};
        for parm in self.pg.nodes:
          opt = self._add_individual_parm(parm);
          sopts.append(opt);
          self._option_list.append(opt);
        self.solve_by_subgroup = False;
        self.solve_by_individual = True;
      # else provide a common equivalent applying to all parms
      else:
        # we do have individual solvability toggles
        # add controls for subgroups
        self._parm_solvability_toggles = dict();
        if self.pg.subgroups:
          subgroup_optlist = [];
          self._subgroup_opts = {};
          for subgroup in self.pg.subgroups:
            subgroup.controller = self;
            opts = self._subgroup_opts[subgroup.name] = record();
            opts.tdloption_namespace = (self.label+"."+subgroup.name).replace(":","_");
            opt = TDLOption('solvable',"Solve for %s"%subgroup.name,True,namespace=opts);
            opt.when_changed(subgroup.set_solvable);
            subgroup_optlist.append(opt);
          subgroup_toggles = TDLMenu('Select solvability by subgroup',
              toggle='solve_by_subgroup',default=True,namespace=self,
              *subgroup_optlist);
          self._option_list.append(subgroup_toggles);
          sopts.append(subgroup_toggles);
        else:
          self.solve_by_subgroup = False;
          subgroup_toggles = None;
        # add individual toggles for parms
        if not self.pg.individual_toggles:
          self.solve_by_individual = False;
        else:
          self._parm_solvability_toggles = dict();
          self._parm_opts = {};
          parm_sopts = [];
          for parm in self.pg.nodes:
            opt = self._add_nonindividual_parm(parm);
            parm_sopts.append(opt);
          if subgroup_toggles:
            parm_toggles = TDLMenu('Select solvability of individual parms',
              toggle='solve_by_individual',default=False,namespace=self,
              *parm_sopts);
            def parm_toggles_changed (enabled):
              if enabled:
                subgroup_toggles.set(False);
            def subgroup_toggles_changed (enabled):
              if enabled:
                parm_toggles.set(False);
            parm_toggles.when_changed(parm_toggles_changed);
            subgroup_toggles.when_changed(subgroup_toggles_changed);
          else:
            parm_toggles = TDLMenu('Toggle solvability of individual parms',*parm_sopts);
            self.solve_by_individual = True;
          self._option_list.append(parm_toggles);
          sopts.append(parm_toggles);
        self._option_list.append(TDLOption('initial_value',
                                  'Override initial value',[None,0.],more=float,namespace=self));
        so = [  TDLOption('time_deg','Polynomial degree, time',[0,1,2,3],more=int,namespace=self),
                TDLOption('freq_deg','Polynomial degree, freq',[0,1,2,3],more=int,namespace=self),
                TDLOption('subtile_time','Solution subinterval (subtile), time',[None],more=int,namespace=self),
                TDLOption('subtile_freq','Solution subinterval (subtile), freq',[None],more=int,namespace=self)
        ];
        sopts += so;
        self._option_list += so;
      # now, more common options which are the same for individual/non-individual parms
      self._option_list += [
        TDLMenu("Use non-default MEP table (default is %s%s)"%
            (self.pg.default_table_name," in MS dir" if pg._table_in_ms else ""),             TDLOption('nondefault_meptable','MEP table',TDLDirSelect("*mep",default=pg.table_name),namespace=self),
          toggle='use_nondefault_meptable',deafult=False,namespace=self),
        TDLOption('use_mep','Initialize with solution from MEP table',True,namespace=self),
        TDLOption('ignore_time',"...even if time domains don't match (e.g. in case of calibrators)",False,namespace=self)
      ];
      so = [  TDLOption('use_previous','Start from solution of previous time interval',True,namespace=self),
              TDLOption('save_all','Save solutions to MEP table even if not converged',True,namespace=self),
              TDLMenu('Parameter constraints (the HMS Suicidal Insanity submenu)',
                TDLOption('force_positive','Force parameters (or c00) to stay positive',False,namespace=self),
                TDLOption('constrain_min','Lower bound for parameter (or c00)',
                    [None],more=float,namespace=self),
                TDLOption('constrain_max','Upper bound for parameter (or c00)',
                    [None],more=float,namespace=self)
              )
      ];
      sopts += so;
      self._option_list += so;
      self._option_list.append(TDLJob(self._clear_mep_tables,"Clear out all previous solutions from MEP tables",job_id=self.label+"_clear_meptables"));
      # update option values from keyword arguments
      for opt in self._option_list:
        if hasattr(opt,'symbol') and opt.symbol in kw:
          opt.set_value(kw[opt.symbol]);
      # now make menu
      self._optmenu = TDLMenu("Solve for %s"%pg.name,
                        toggle='solvable',open=solvable,namespace=self,
                        *self._option_list);
      # show/hide solvability-related options
      def show_hide_sopts (show):
        for opt in sopts:
          opt.show(show);
      self._optmenu.when_changed(show_hide_sopts);

    def _add_individual_parm (self,parm):
      opts = self._parm_opts[parm.name] = record();
      opts.tdloption_namespace = (self.label+"."+parm.name).replace(":","_");
      # solvable parm options
      parm_sopts = [
            TDLOption('time_deg','Polynomial degree, time',[0,1,2,3],more=int,namespace=opts),
            TDLOption('freq_deg','Polynomial degree, freq',[0,1,2,3],more=int,namespace=opts),
            TDLOption('subtile_time','Solution subinterval (subtile), time',[None],more=int,namespace=opts),
            TDLOption('subtile_freq','Solution subinterval (subtile), freq',[None],more=int,namespace=opts)
      ];
      parm_solv = TDLMenu("Solve for %s"%parm.name,
            TDLOption('initial_value','Override initial value',[None,0.],more=float,namespace=opts),
            toggle='solvable',open=True,namespace=opts,
            *parm_sopts
          );
      # hide solvables
      def show_hide_sopts (show):
        # print 'showing',show,parm_sopts;
        for opt in parm_sopts:
          opt.show(show);
      parm_solv.when_changed(show_hide_sopts);
      return parm_solv;

    def _add_nonindividual_parm (self,parm):
      opts = self._parm_opts[parm.name] = record();
      opts.tdloption_namespace = (self.label+"."+parm.name).replace(":","_");
      return TDLOption("solvable","Solve for %s"%parm.name,True,namespace=opts);

    def runtime_options (self,submenu=True):
      return [ self._optmenu ];

    def _fill_state_record (self,state):
      """Fills a state record with common options""";
      state.table_name    = self.get_table_name();
      state.reset_funklet = not self.use_mep;
      state.ignore_time   = self.ignore_time;
      state.use_previous  = self.use_previous;
      state.save_all      = self.save_all;

    def make_cmdlist (self,solvable=True):
      """Makes a record suitbale for inclusion in a request command_by_list entry.
      If 'solvable' is True, solvability is determined by parm options.
      If 'solvable' is False, solvability is always False.
      """;
      solvable = solvable and self.solvable;
      if self._individual:
        cmdlist = [];
        deferred = [];
        for parm in self.pg.nodes:
          opts = self._parm_opts[parm.name];
          # if solvable, make a full record for this parm
          if solvable and opts.solvable:
            state = record(solvable=True);
            if opts.initial_value is not None:
              state.default_value = opts.initial_value;
            state.shape         = [opts.time_deg+1,opts.freq_deg+1];
            state.tiling        = record();
            if opts.subtile_time:
              state.tiling.time = opts.subtile_time;
            if opts.subtile_freq:
              state.tiling.freq = opts.subtile_freq;
            self._fill_state_record(state);
            cmdlist.append(record(name=[parm.name],state=state));
          # else if value is overridden, also make a full record
          elif opts.initial_value is not None:
            state = record(solvable=False);
            state.default_value = opts.initial_value;
            state.shape = [0];
            self._fill_state_record(state);
            cmdlist.append(record(name=[parm.name],state=state));
          # else include in list of other parms to be initialized with a single record later
          else:
            deferred.append(parm);
        # and now make an entry for the non-solvables
        if deferred:
          state = record(solvable=False,shape=[0]);
          self._fill_state_record(state);
          cmdlist.append(record(name=[parm.name for parm in deferred],state=state));
        return cmdlist;
      # else all parms are treated together, so make one state record for all
      # solvables
      else:
        cmdlist = [];
        if solvable:
          if self.solve_by_subgroup:
            solvables = [];
            nonsolvables = [];
            for parm in self.pg.nodes:
              solvable = True;
              for sg in self.pg.parm_subgroups.get(parm.name,[]):
                if not sg.solvable:
                  solvable = False;
                  break;
              if solvable:
                solvables.append(parm.name);
              else:
                nonsolvables.append(parm.name);
          elif self.solve_by_individual:
            solvables = [ parm.name for parm in self.pg.nodes if self._parm_opts[parm.name].solvable ];
            nonsolvables = [ parm.name for parm in self.pg.nodes if not self._parm_opts[parm.name].solvable ];
          else:
            solvables = [parm.name for parm in self.pg.nodes ];
            nonsolvables = [];
        else:
          solvables = [];
          nonsolvables = [ parm.name for parm in self.pg.nodes ];
        # command for solvables
        if solvables:
          state = record(solvable=True);
          if self.initial_value is not None:
            state.default_value = self.initial_value;
          state.shape         = [self.time_deg+1,self.freq_deg+1];
          state.tiling        = record();
          if self.subtile_time:
            state.tiling.time = self.subtile_time;
          if self.subtile_freq:
            state.tiling.freq = self.subtile_freq;
          self._fill_state_record(state);
          cmdlist.append(record(name=solvables,state=state));
        if nonsolvables:
          state = record(solvable=False);
          if self.initial_value is not None:
            state.default_value = self.initial_value;
          state.shape = [0];
          self._fill_state_record(state);
          cmdlist.append(record(name=nonsolvables,state=state));
        return cmdlist;

    def _clear_mep_tables (self,mqs,parent,**kw):
      tabname = self.get_table_name();
      if GUI.warning_box("Clearing solutions",
	 "This will clear out <b>all</b> previous solutions from table '%s'. Are you sure you want to do this?"%tabname,
          GUI.Button.Yes|GUI.Button.No,GUI.Button.No) != GUI.Button.Yes:
         return;
      try:    os.system("rm -fr "+tabname);
      except: pass;

    def get_table_name (self):
      if self.use_nondefault_meptable and self.nondefault_meptable:
        return self.nondefault_meptable;
      elif self.pg._table_in_ms and Context.mssel and Context.mssel.msname:
        return os.path.join(Context.mssel.msname,self.pg.default_table_name);
      else:
        return self.pg.default_table_name;

  def _sort_members (members):
    sorted_members = list(members);
    # define comparison function that sorts by qualifier as well
    def int_or_str(x):
      try: return int(x);
      except: return x;
    from past.builtins import cmp
    from functools import cmp_to_key
    sorted_members.sort(key=cmp_to_key(lambda a,b:
        cmp(list(map(int_or_str,a.name.split(':'))),list(map(int_or_str,b.name.split(':'))))));
    return sorted_members;
  _sort_members = staticmethod(_sort_members);

  def __init__ (self,label,members=[],name=None,
                     subgroups=None,
                     individual=False,individual_toggles=True,bookmark=True,
                     table_name="calibration.mep",table_in_ms=True,
                     **kw):
    """Creates a ParmGroup.
    'label' is the label of the group
    'members' is a list of MeqParm nodes
    'name'  is a descriptive name for the group. If not set, label is used.
    'subgroups' is a list of Subgroup objects, if parms are to be subgrouped
    'individual' is True if individual options are to be provided per parameter. Note that
                this is incompatible with subgroups
    'individual_toggles' is True if individual solvasbility toggles are to be provided in addition to
                subgroups.
    'bookmark' is True if bookmarks for the parameters are to be automatically provided
    'table_name' is the name of a MEP table
    'table_in_ms' is True if the table should reside inside the current MS (in which case the MS selector
                from the global Meow.Context is invoked to get the name of the MS)
    """
    self.label = label;
    self.name  = name or label;
    # sort group members by name
    self.nodes = sorted_nodes = self._sort_members(members);
    # now collect list of subgroups
    self.subgroups = subgroups or [];
    self.individual_toggles = individual_toggles;
    if subgroups and individual:
      raise TypeError("cannot combine subgroups and individual parameters");
    self.parm_subgroups = dict();
    for sg in self.subgroups:
      for parm in sg.nodes:
        self.parm_subgroups.setdefault(parm.name,[]).append(sg);
    self.nodes = sorted_nodes;
    # setup table names
    self.table_name = table_name;
    self.default_table_name = os.path.basename(table_name);
    if table_in_ms and Context.mssel and Context.mssel.msname:
      self.table_name = os.path.join(Context.mssel.msname,self.default_table_name);
    self._table_in_ms = table_in_ms;
    self._individual = individual;
    # put table name into the parms
    for node in self.nodes:
      node.initrec().table_name = self.table_name;
    # create bookmarks (if specified as a [W,H], it gives the number of parms to bookmark)
    if bookmark:
      if isinstance(bookmark,tuple) and len(bookmark) == 2:
        ncol,nrow = bookmark;
      else:
        ncol,nrow = 2,3;
      Bookmarks.make_node_folder("Parameters: %s"%label,
                        sorted_nodes,sorted=True,ncol=ncol,nrow=nrow);
    # add ourselves to global list of parmgroups
    global _all_parmgroups;
    _all_parmgroups.append(self);

  def add (self,*members):
    self.nodes += list(members);

  def make_controller (self,label=None,**kw):
    return self.Controller(self,label,**kw);

_all_solvejobs = [];

class SolveJob (object):
  def __init__ (self,label,name,*active_parmgroups):
    self.label = label;
    self.name  = name;
    self.tdloption_namespace = label.replace(":","_").replace(" ","_");
    self.active_parmgroups = active_parmgroups;
    # menu is only made on-demand
    self._jobmenu = None;
    # add to global list
    _all_solvejobs.append(self);

  def runtime_options (self):
    if self._jobmenu is None:
      opts = [ TDLOption('tile_size',"Tile size, in timeslots",[1,10,100],
                    doc="""Input data is normally sliced by time, and processed in chunks of the
                    indicated size. This will also be the effective parameter solution interval
                    (in time), unless you specify a different (smaller) value for the "Solution subinterval (time)" option below.""",
                    more=int,namespace=self),
               TDLOption('time_step',"Time stepping, in timeslots",[1,2,5,10],
                    doc="""Enter a step size N>1 to only process every Nth timeslot.""",
                    more=int,namespace=self)      ];
      # add options from parmgroups
      global _all_parmgroups;
      self.pg_controllers = [];
      other_opts = [];
      for pg in _all_parmgroups:
        solvable = pg in self.active_parmgroups;
        controller = pg.make_controller(self.label,solvable=solvable);
        if solvable:
          opts += controller.runtime_options();
        else:
          other_opts += controller.runtime_options();
        self.pg_controllers.append(controller);
      if other_opts:
        opts.append(TDLMenu("Simultaneously solve for other parameters",*other_opts));

      # add solver control
      self.solver_control = SolverControl(self.label);
      opts.append(TDLMenu("Solver options (for the brave)",*self.solver_control.runtime_options()));
      # add solve job
      opts.append(TDLJob(self._run_solve_job,self.name or "Run solution",job_id=self.label));
      # now make a runtime menu
      self._jobmenu = TDLMenu(self.name or "Solve for %s"%self.label,*opts);
    return [ self._jobmenu ];

  def run_solution (self,mqs,mssel=None,vdm=None,tiling=None,time_step=1,wait=False):
    """Helper function to put together TDL jobs.
    Starts a solution, setting the group to solvable""";
    mssel = mssel or Context.mssel;
    # make command lists for our parameters
    cmdlist = [];
    for pgc in self.pg_controllers:
      cmdlist += pgc.make_cmdlist();
    self.solver_control.update_state(mqs,cmdlist=cmdlist);
    # run the VisDataMux
    vdm = namify(vdm or Context.vdm or 'VisDataMux')
    return mqs.execute(vdm,mssel.create_io_request(tiling,time_step=time_step),wait=wait);

  def _run_solve_job (self,mqs,parent,wait=False,**kw):
    return self.run_solution(mqs,tiling=self.tile_size,time_step=self.time_step,wait=wait);

def num_solvejobs ():
  global _all_solvejobs;
  return len(_all_solvejobs);

def get_solvejob_options ():
  global _all_solvejobs;
  opts = [];
  for job in _all_solvejobs:
    opts += job.runtime_options();
  return opts;
