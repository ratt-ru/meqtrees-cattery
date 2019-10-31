# -*- coding: utf-8 -*-
#
#% $Id$
#
#
# Copyright (C) 2002-2007
# The MeqTree Foundation &
# ASTRON (Netherlands Foundation for Research in Astronomy)
# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>,
# or write to the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from Timba.TDL import *
from . import Context
from . import Jones
from . import Bookmarks
from . import Utils

class _BaseTree (object):
  def __init__ (self,ns,array=None,observation=None):
    self.ns = ns;
    self.array = array or Context.array;
    self.observation = observation or Context.observation;
    if not self.array or not self.observation:
      raise ValueError("array or observation not specified in global Meow.Context, or in this function call");
    self._inputs = self._outputs = None;

  def set_inputs (self,inputs):
    if self._inputs is not None:
      if inputs is not None and inputs is not self._inputs:
        raise ValueError("this tree is already using a different set of inputs");
    else:
      self._inputs = inputs or self.array.spigots(flag_bit=1);

  def inputs (self):
    return self._inputs;

class ResidualTree (_BaseTree):
  def __init__ (self,ns,predict,array=None,observation=None):
    _BaseTree.__init__(self,ns,array,observation);
    self._predict = predict;

  def residuals (self,inputs=None):
    """Makes residual trees, by subtracting predict(p,q) from inputs(p,q).
    Returns residual node, which must be qualified with a (p,q) pair.
    """;
    # figure out our inputs
    self.set_inputs(inputs);
    outputs = self.ns.residual;
    if not outputs(*self.array.ifrs()[0]).initialized():
      for p,q in self.array.ifrs():
        outputs(p,q) << self._inputs(p,q) - self._predict(p,q);
    return outputs;

class SolveTree (ResidualTree):
  def __init__ (self,ns,predict,solve_ifrs=None,array=None,observation=None,
        weights=None,modulo=None,
        residuals=None,bookmarks=True):
    ResidualTree.__init__(self,ns,predict,array,observation);
    self._solve_ifrs = solve_ifrs or self.array.ifrs();
    self._make_residuals = residuals;
    self._weights = weights;
    self._modulo = modulo or 0;
    self.bookmarks = bookmarks;

  def outputs (self,inputs=None):
    """Makes solver tree, returns the set of reqseq output nodes.
    This is for compatibility with older scripts. The output nodes will be
    reqseqs returning either inputs or residuals, as chosen in the
    constructor.
    """;
    self.set_inputs(inputs);
    if self._make_residuals:
      outputs = self.residuals(inputs);
    else:
      outputs = self._inputs;
    return self.sequencers(inputs,outputs=outputs);

  def solver (self,inputs=None):
    """Makes solver tree given the inputs and the predicts, returns
    solver node.
    """;
    self.set_inputs(inputs);
    solver = self.ns.solver;
    if not solver.initialized():
      self._solver_name = solver.name;
      # create condeqs and request sequencers
      for p,q in self._solve_ifrs:
        children = [ self._inputs(p,q),self._predict(p,q) ];
        children[0].initrec().cache_num_active_parents = 1;
        if self._weights:
          children.append(self._weights(p,q));
        self.ns.ce(p,q) << Meq.Condeq(children=children,modulo=self._modulo);
      # create optimal poll order for condeqs, for efficient parallelization
      # (i.e. poll child 1:2, 3:4, 5:6, ..., 13:14 then the rest)
      # however, since _solve_ifrs may be a subset, we'll have to be a bit more clever
      polled = set();
      poll_order = [];
      for p,q in self._solve_ifrs:
        if p not in polled and q not in polled:
          poll_order.append(self.ns.ce(p,q).name);
          polled.add(p);
          polled.add(q);
          if len(poll_order) >= len(self.array.stations())/2:
            break;
      # add condeqs to solver
      solver << Meq.Solver(children=[self.ns.ce(p,q) for p,q in self._solve_ifrs],
                           child_poll_order=poll_order,flush_tables=True);
      # add solver bookmark
      if self.bookmarks:
        Bookmarks.Page("Solver").add(solver,viewer="Result Plotter");
    return solver;

  def sequencers (self,inputs=None,outputs=None):
    """Makes solver tree and adds sequencers to execute an
    'outputs' branch after the solve is complete.
    """;
    # figure out our inputs
    self.set_inputs(inputs);
    if outputs is None:
      raise TypeError("mandatory 'outputs' argument missing");

    reqseq = self._outputs = self.ns.reqseq;
    if not reqseq(*self.array.ifrs()[0]).initialized():
      solver = self.solver();
      # create condeqs and request sequencers
      for p,q in self.array.ifrs():
        reqseq(p,q) << Meq.ReqSeq(solver,outputs(p,q),result_index=1);
        outputs(p,q).set_options(cache_num_active_parents=1);

    return reqseq;

  def define_solve_job (self,jobname,jobid,solvables,tile_sizes=[1,10,100],vdm=None):
    """Defines a "solve job", with a runtime menu entry and its own set of options
    'jobname':    a descriptive name (will appear in the menu)
    'jobid':      unique identifier string
    'solvables':  list of solvables (Parm nodes or node names)
    'tile_sizes': list of suggested tilings
    'vdm':        VisDataMux node which runs things, Context.vdm is used if None
    """;
    # make sure solvables is a list of names
    def namify (arg):
      if isinstance(arg,str):
        return arg;
      elif is_node(arg):
        return arg.name;
      else:
        raise TypeError("'solvables' contains an object of illegal type '"+str(type(arg))+"'");
    solvables = [ namify(x) for x in solvables ];

    # init option vars
    tiling = 'tiling_'+jobid;
    solver_opts = '_solver_'+jobid;
    globals()[tiling] = 1;

    # figure out VDM name
    vdm = vdm or Context.vdm;
    if vdm is None:
      vdm = 'VisDataMux';
    elif is_node(vdm):
      vdm = vdm.name;
    elif not isinstance(vdm,str):
      raise TypeError("'vdm' argument must be a node or a node name");

    # define TDL job to run this solve job
    def run_solve_job (mqs,parent,**kw):
      # put together list of enabled solvables
      Utils.run_solve_job(mqs,solvables,
          solver_node=self.ns.solver.name,
          vdm_node=vdm,
          tiling=globals()[tiling],
          options=getattr(Utils,solver_opts));

    # define submenu
    TDLRuntimeMenu(jobname,
      TDLOption(tiling,"Tile size",tile_sizes,more=int),
      TDLMenu("Solver options",*Utils.solver_options(solver_opts)),
      TDLJob(run_solve_job,"Run solution")
    );

def define_inspector (nodeseries,qlist,*qualifier_lists,**kw):
  """Defines an inspector node for a node series. Returns a definition record, so it should\
  be <<'d into a node. E.g.:
      ns.my_inspector << define_inspector(mynodes,sources,stations);
  'nodeseries' should be an under-qualified node.
  All remaining arguments are treated as lists of qualifiers to be looped over. Note that these
  may also be objects with a "name" attribute (e.g. Meow.SkyComponent or Meow.Direction).
  Optional keywords:
    label (=None)     used to supply a base plot label; if not given, then nodeseries.name is used.
    freqavg(=True)    if True, frequency averging is performed on each node in the series. Set to False
                      if you know your nodes do not have a freq axis.
  """;
  label = kw.get('label',None) or getattr(nodeseries,'name','');
  freqavg = kw.get('freqavg',True);
  plot_labels = [];
  plot_children = [];
  # generate list of all qualifier combinations
  all_quals = [[getattr(q,'name',None) or str(q)] for q in qlist];
  for ql in qualifier_lists:
    all_quals = [ qq+[getattr(q,'name',None) or str(q)] for qq in all_quals for q in ql ];
  # filter list so that we only use nodes that are actually defined
  all_quals = [ qq for qq in all_quals if nodeseries(*qq).initialized() ];
  # create list of plot labels
  plot_labels = [ ":".join([label]+qq) for qq in all_quals ];
  if freqavg:
    children = [ nodeseries(*(qq+['freqavg'])) << Meq.Mean(nodeseries(*qq),reduction_axes=["freq"])
                  for qq in all_quals ];
  else:
    children = [ nodeseries(*qq) for qq in all_quals ];
  # now make a composer
  return Meq.Composer(plot_label=plot_labels,dims=[0],*children);

def inspector (outnode,nodes,bookmark=True):
  """Makes a generic inspector for a list of nodes.
  'outnode' is an output node to which the inspector is assigned.
  'nodes' is a list of nodes.
  If 'bookmark' is true, a single-page bookmark will automatically be added
  for the inspector (use string to give it a non-default name.)
  """
  if not nodes:
    raise ValueError("too few nodes specified in list");
  outnode << \
    Meq.Composer(
      plot_label=[ node.name for node in nodes ],dims=[0],
      mt_polling=True,
      *[ Meq.Mean(node,reduction_axes="freq") for node in nodes ]
    );
  if bookmark is True:
    bookmark = outnode.name;
  if bookmark:
    Bookmarks.Page(bookmark).add(outnode,viewer="Collections Plotter");
  return outnode;


def vis_inspector (outnode,visnodes,ifrs=None,array=None,bookmark=True,**kw):
  """Makes an inspector for visibility nodes.
  'outnode' is an output node to which the inspector is assigned.
  'visnodes' will be qualified with the ifrs pairs from the given array.
  If 'bookmark' is true, a single-page bookmark will automatically be added
  for the inspector (use string to give it a non-default name.)
  """
  if ifrs is None:
    array = array or Context.array;
    if not array:
      raise ValueError("array or ifrs not specified in global Meow.Context, or in this function call");
    ifrs = array.ifrs();
  outnode << \
    Meq.Composer(
      dims=[0],
      plot_label=[ "%s-%s"%(p,q) for p,q in ifrs ],
      mt_polling=True,
      *[ outnode(p,q) << Meq.Mean(visnodes(p,q),reduction_axes="freq")
         for p,q in ifrs ],**kw
    );
  if bookmark is True:
    bookmark = outnode.name.replace("_"," ");
  if bookmark:
    Bookmarks.Page(bookmark).add(outnode,viewer="Collections Plotter");
  return outnode;


def jones_inspector (outnode,jones,array=None,bookmark=True):
  """Makes an inspector for Jones matrices.
  'outnode' is an output node to which the inspector is assigned.
  'jones' will be qualified with the stations pairs from the given array.
  If 'bookmark' is true, a single-page bookmark will automatically be added
  for the inspector (use string to give it a non-default name.)
  """
  array = array or Context.array;
  if not array:
    raise ValueError("array not specified in global Meow.Context, or in this function call");
  outnode << \
    Meq.Composer(
      dims=[0],
      plot_label=[ str(p) for p in array.stations() ],
      mt_polling=True,
      *[ outnode(p) << Meq.Mean(jones(p),reduction_axes="freq")
         for p in array.stations() ]
    );
  if bookmark is True:
    bookmark = outnode.name;
  if bookmark:
    Bookmarks.Page(bookmark).add(outnode,viewer="Collections Plotter");
  return outnode;

def make_clipping_flagger (node,minval=None,maxval=None,flagmask=1):
  flaggers = [];
  if minval is not None:
    flaggers.append(node('clip_lt') << Meq.ZeroFlagger(Meq.Abs(node)-minval,oper='lt',flag_bit=flagmask));
  if maxval is not None:
    flaggers.append(node('clip_gt') << Meq.ZeroFlagger(Meq.Abs(node)-maxval,oper='gt',flag_bit=flagmask));
  return node('clipped') << Meq.MergeFlags(node,*flaggers);

def make_jones_norm_flagger (jp,minval=None,maxval=None,jnorm=None,jflag=None,freqmean=False,flagmask=1):
  """Makes nodes to compute the nporm of a Jones matrix, and a flagger based on clipping this to minval/maxval.
  jpnorm, if supplied, is the base node for the matrix norm -- else a default base node is made by using jp('norm'),""";
  flaggers = [];
  if jnorm is None:
    jnorm = jp('norm');
  jpsq = jp('sq') << Meq.MatrixMultiply(jp,jp('conj') << Meq.ConjTranspose(jp));
  jj = jp('2tr') << (jpsq(11) << Meq.Selector(jpsq,index=0)) + (jpsq(22) << Meq.Selector(jpsq,index=3));
  jj = jp('a2tr') << Meq.Abs(jj);
  jnorm << ( Meq.Mean(Meq.Sqrt(jj),reduction_axes="freq") if freqmean else Meq.Sqrt(jj) );
  if minval is not None:
    flaggers.append(jp('clip_lt') << Meq.ZeroFlagger(jnorm-minval,oper='lt',flag_bit=flagmask));
  if maxval is not None:
    flaggers.append(jp('clip_gt') << Meq.ZeroFlagger(jnorm-maxval,oper='gt',flag_bit=flagmask));
  if flaggers:
    return (jflag if jflag is not None else jp('clipped')) << Meq.MergeFlags(jp,*flaggers);
  else:
    return jp;
    
def make_jones_abs_flagger (jp,minval=None,maxval=None,jabs=None,jflag=None,diag=True,freqmean=False,flagmask=1):
  """Makes nodes to compute the absolute value of a Jones matrix, and a flagger based on clipping this to minval/maxval.
  jpabs, if supplied, is the base node for the matrix norm -- else a default base node is made by using jp('abs'),""";
  flaggers = [];
  if jabs is None:
    jabs = jp('abs');
  jabs << ( Meq.Mean(Meq.Abs(jp),reduction_axes="freq") if freqmean else Meq.Abs(jp) );
  if diag:
    jabs = jabs('sel') << Meq.Selector(jabs,index=[0,3],multi=True);
  if minval is not None:
    flaggers.append(jp('clip_lt') << Meq.ZeroFlagger(jabs-minval,oper='lt',flag_bit=flagmask));
  if maxval is not None:
    flaggers.append(jp('clip_gt') << Meq.ZeroFlagger(jabs-maxval,oper='gt',flag_bit=flagmask));
  if flaggers:
    return (jflag if jflag is not None else jp('clipped')) << Meq.MergeFlags(jp,merge_all=True,*flaggers);
  else:
    return jp;


def make_sinks (ns,outputs,array=None,
                post=None,
                vdm=None,
                spigots=True,
                output_col='DATA',**kw):
  array = array or Context.array;
  if not array:
    raise ValueError("array not specified in global Meow.Context, or in this function call");
  # make sinks
  sink = array.sinks(children=outputs,output_col=output_col,**kw);
  # make vdm
  if vdm is None:
    vdm = ns.VisDataMux;
  elif not is_node(vdm):
    raise TypeError("'vdm' argument should be a node");
  Context.vdm = vdm;

  # check 'post' argument, if it's a list, make a ReqMux called vdm:post
  # to process these guys
  if isinstance(post,(list,tuple)):
    if post:
      post = vdm('post') << Meq.ReqMux(*post);
    else:
      post = None;
  elif post and not is_node(post):
    raise TypeError("'post' argument should be a node or a list of nodes, %s given"%type(post));

  # figure out optimal poll order for children. We want to poll things 
  # so that per-antenna trees are parallelized as much as possible, so first let's poll
  # those IFRs that don't share antennas
  cpo = [];
  seen_stations = set();
  for p,q in array.ifrs():
    if p not in seen_stations and q not in seen_stations:
      cpo.append(sink(p,q,).name);
      seen_stations.add(p);
      seen_stations.add(q);
      if len(seen_stations) >= array.num_stations()-1:
        break;
  # now make the vdm
  vdm << Meq.VisDataMux(post=post,child_poll_order=cpo or None,*[sink(p,q) for p,q in array.ifrs()]);
  if spigots:
    if spigots is True:
      spigots = array.spigots();
    vdm.add_stepchildren(*[spigots(p,q) for p,q in array.ifrs()]);

  return vdm;

