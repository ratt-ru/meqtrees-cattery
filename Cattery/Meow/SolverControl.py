# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from Timba.TDL import *


def namify (arg):
  if isinstance(arg,str):
    return arg;
  elif is_node(arg):
    return arg.name;
  else:
    raise TypeError("node name or node object expected, got '%s'"%str(type(arg)));



class SolverControl (object):
  def __init__ (self,label,name=None,node='solver'):
    self.label = label;
    self.name  = name or label;
    self.set_solver_node(node);
    # add TDLOptions
    class OptionsHolder (object):
      pass;
    self._opts = OptionsHolder();
    self._opts.tdloption_namespace = label+'.solver';
    self._option_list = [
      TDLOption('debug_level',"Solver debug level",[0,1,10],namespace=self._opts),
      TDLOption('colin_factor',"Collinearity factor",[0,1e-8,1e-6,1e-3,1e-1],default=0,more=float,namespace=self._opts),
      TDLOption('lm_factor',"Initial LM factor",[1,.1,.01,.001],default=3,more=float,namespace=self._opts),
      TDLOption('balanced_equations',"Assume balanced equations",False,namespace=self._opts),
      TDLOption('epsilon',"Convergence threshold",[.001,.0001,1e-5,1e-6],default=2,more=float,namespace=self._opts),
      TDLOption('num_iter',"Max iterations",[5,10,15,30,50,100,1000],default=2,more=int,namespace=self._opts),
      TDLOption('convergence_quota',"Subtiling convergence quota",[.9,1.],more=float,namespace=self._opts), \
      TDLOption('mt_solve',"Multithread solver if possible",True,namespace=self._opts) \
    ];

  def runtime_options (self):
    return self._option_list;

  def set_solver_node (self,node):
    self.node = namify(node);

  def make_state_record (self,solvables=None,cmdlist=None,**kw):
    # copy all options into solver defaults with the same name
    solver_defaults = record();
    for opt in self._option_list:
      solver_defaults[opt.symbol] = opt.value;
    # copy any extra options from functiohn invocation
    solver_defaults.update(kw);
    solver_defaults.epsilon_deriv      = solver_defaults["epsilon"];
    solver_defaults.save_funklets    = True;
    solver_defaults.last_update      = True;
    # add a solvable command
    # either 'solvables' is specified by name, or a cmdlist record/list of records is given
    if solvables:
      solvables = list(map(namify,solvables));
      solver_defaults.solvable         = record(command_by_list=(record(name=solvables,
                                                state=record(solvable=True)),
                                                record(state=record(solvable=False))))
    elif cmdlist:
      if not isinstance(cmdlist,(list,tuple)):
        cmdlist = (cmdlist,record(state=record(solvable=False)));
      else:
        cmdlist = list(cmdlist) + [ record(state=record(solvable=False)) ];
      solver_defaults.solvable = record(command_by_list=cmdlist);
    return solver_defaults;

  def update_state (self,mqs,solvables=None,cmdlist=None,wait=False,sync=True,**kw):
    state = self.make_state_record(solvables=solvables,cmdlist=cmdlist,**kw);
    mqs.setnodestate(self.node,state,wait=wait,sync=sync);

