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
from Timba.Meq import meq
from Timba.array import *
from Timba import mequtils
from Timba.dmi import hiid

import os
import random

import Meow

from Meow import Bookmarks
from Kittens import utils

from Calico import ParmTables
parmtab = parmtab_name = None;


def reload_parm_table (name):
  """This function is called whenever a new parmtable is selected.""";
  global parmtab;
  global parmtab_name;
  if name and ( parmtab is None or parmtab_name != name ):
    box = GUI.MessageBox("Reading parmtable","""<P>We're now reading the parmtable %s in order to
        determine what kind of parameters it contains. This may take a few seconds.</P>"""%name,GUI.Information,
        GUI.Button.NoButton,GUI.Button.NoButton);
    busy = GUI.BusyIndicator();
    box.show();
    try:
      parmtab = ParmTables.ParmTab(name);
      parmtab_name = name;
    finally:
      box.hide();
      busy = None;
    
# Parm options first
ptopt = TDLCompileOption('parmtable_name',"Select parmtable",TDLDirSelect("*.fmep"),mandatory=True);
ptopt.when_changed(reload_parm_table);

TDLRuntimeOption("tile_size","Tile size, in timeslots",[None,60,120,240],more=int);

# define some methods for forming up jones matrices from indivisdual parms
def jones_diag (xx,yy):
  return Meq.Matrix22(xx,0,0,yy);

def jones_full (*children):
  return Meq.Matrix22(*children);

# This list defines all the possible functions for building aggreaget objects such as Jones matrices 
# and complex numbers from individual parameters.
# For any given set of name qualifiers, the code below will attempt to find a function in the list to 
# map a set of parameters to compound objects.
# This function is expected to return a node definition.

aggregators = [
  (('r','i'),Meq.ToComplex),
  (('ampl','phase'),Meq.Polar),
  (('xx','xy','yx','yy'),Meq.Matrix22),
  (('XX','XY','YX','YY'),Meq.Matrix22),
  (('rr','rl','lr','ll'),Meq.Matrix22),
  (('RR','RL','LR','LL'),Meq.Matrix22),
  (('xx','yy'), lambda xx,yy:Meq.Matrix22(xx,0,0,yy)),
  (('rr','ll'), lambda rr,ll:Meq.Matrix22(rr,0,0,ll)), 
  (('XX','YY'), lambda xx,yy:Meq.Matrix22(xx,0,0,yy)),
  (('RR','LL'), lambda rr,ll:Meq.Matrix22(rr,0,0,ll)) 
];



def _define_forest (ns,**kw):
  if parmtab is None:
    raise RuntimeError("please select a valid parmtable in compile-time options");
  # make one MeqParm per each funklet
  parmdef = Meq.Parm(table_name=parmtable_name);
  parmnodes = dict([(name,(ns[name] << parmdef)) for name in parmtab.funklet_names()]);

  # reprocess flag is reset every time something changes in the parmnodes dict,
  # so we simply repeat this loop until things stop changing
  reprocess = True;
  while reprocess:
    reprocess = False;
    # see if we have recipe for making compound objects based on last funklet name qualifier
    nqlist = [ name.rsplit(':',1) for name in parmnodes.keys() if ':' in name ];
    # set of basenames and set of all qualifiers
    basenames = set([nq[0] for nq in nqlist]);
    quals = set([nq[1] for nq in nqlist]);
    for quallist,agg_func in aggregators:
      # look for an aggregator whose qualifiers are a subset of the given set
      if set(quallist) <= quals:
        # for each basename, attempt to find funklets for every qualifier in the list
        for basename in basenames:
          fqnames = [ basename+':'+qq for qq in quallist ];
          funks   = [ parmnodes.get(name,None) for name in fqnames ];
          # if every basename:qq combination gives a valid funklet, call the aggregator
          if None not in funks:
            aggnode = ns[basename] << agg_func(*funks);
            # remove aggregated funklets from the parmnodes dict
            for name in fqnames:
              del parmnodes[name];
            # add new node to dict
            parmnodes[basename] = aggnode;
            reprocess = True;

  # sort nodes by qualified name
  nodes = list(parmnodes.items());
  from functools import cmp_to_key
  nodes.sort(key=cmp_to_key(lambda a,b:ParmTables.cmp_qualified_names(a[0],b[0])));
  nodes = [ namenode[1] for namenode in nodes ];
  
  # make intermediate math node
  pall = ns.all_parms_raw << Meq.Composer(children=nodes,mt_polling=True,dims=[0]);
  pfm = ns.all_parms_fm << Meq.Mean(pall,reduction_axes="freq");
#  ns.all_minus_fm << pall - pfm;
#  ns.all_over_fm << pall/pfm;

  # make root node
  ns.root << Meq.ReqMux(
    pfm
#    ns.all_minus_fm,
#    ns.all_over_fm,
#    ns.fft_all_over_fm << Meq.FFTBrick(Meq.Selector(ns.all_over_fm,index=0,multi=True),axes_in=(hiid("time"),hiid("freq")),axes_out=(hiid("l"),hiid("m")))
  );
  Bookmarks.Page("Composite view").add(ns.all_parms_fm,viewer="Collections Plotter");
  # make bookmark folder
  Bookmarks.make_node_folder("All parameters",nodes,sorted=True);
  # close parmtable, else meqserver will be unable to get at it
  parmtab.close();

  # now make runtime options
  global active_axes;
  active_axes = [];
  reg_domain_opts = [];

  for iaxis in range(mequtils.max_axis):
    stats = parmtab.axis_stats(iaxis);
    if not stats.empty():
      active_axes.append((iaxis,stats.name,stats));
      reg_domain_opts.append( 
        TDLOption("num_cells_%s"%stats.name,"Number of grid points in %s"%stats.name,[len(stats.cells)],more=int)
      );

  TDLRuntimeMenu("View parameters over regularly gridded domain",
    *(reg_domain_opts + [TDLJob(_view_parameters_over_regularized_domain,"View over gridded domain")]));

  
def exec_cells (mqs,node,cells):
  if tile_size:
    time_grid = Timba.array.array(cells.grid.time);
    time_size = Timba.array.array(cells.cell_size.time);
    domid = 0;
    for i0 in range(0,len(time_grid),tile_size):
      i1 = min(len(time_grid),i0+tile_size);
      meq.add_cells_axis(cells,hiid('time'),grid=time_grid[i0:i1],cell_size=time_size[i0:i1]);
      req = meq.request(cells,rqid=meq.requestid(domain_id=domid));
      domid += 1;
      mqs.execute(node,req,wait=False);
  else:
    req = meq.request(cells);
    mqs.execute(node,req,wait=False);
  
  
def _tdl_job_View_parameters_over_intrinsic_domains (mqs,parent,**kw):
  # make domain encompassing entire parmtable
  exec_cells(mqs,'root',parmtab.subdomain_cells());

def _view_parameters_over_regularized_domain (mqs,parent,**kw):
  # make domain encompassing entire parmtable
  cells_per_axis = {};
  for iaxis,name,stats in active_axes:
    cells_per_axis["num_%s"%name] = globals()['num_cells_%s'%name];
  exec_cells(mqs,'root',parmtab.envelope_cells(**cells_per_axis));


