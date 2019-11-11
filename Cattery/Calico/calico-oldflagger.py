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
import os
import random

import Meow

import Meow.Bookmarks
import Meow.Utils
import Meow.StdTrees
import Meow.Context

def abs_clip (inputs,maxval=None,minval=None):
  # flag condition is abs(x) >= threshold
  flagged = inputs("absclip");
  for p,q in Meow.Context.array.ifrs():
    inp = inputs(p,q);
    flaggers = [];
    if maxval is not None:
      fc = flagged("fc1",p,q) << Meq.Abs(inp) - maxval;
      flaggers.append(flagged("zf1",p,q) << Meq.ZeroFlagger(fc,flag_bit=4,oper="GE"));
    if minval is not None:
      fc = flagged("fc2",p,q) << Meq.Abs(inp) - minval;
      flaggers.append(flagged("zf2",p,q) << Meq.ZeroFlagger(fc,flag_bit=4,oper="LE"));
    flagged(p,q) << Meq.MergeFlags(inp,*flaggers);
  return flagged;

def rms_clip (inputs,threshold_sigmas):
  # flag condition is abs(x)-mean(abs(x)) >= flag_threshold*rms
  flagged = inputs("rmsclip");
  for p,q in Meow.Context.array.ifrs():
    inp = inputs(p,q);
    a = flagged("abs",p,q) << Meq.Abs(inp);
    delta = flagged("delta",p,q) << Meq.Abs(a-Meq.Mean(a));
    stddev_a = flagged("stddev",p,q) << Meq.StdDev(delta);
    fc = flagged("fc",p,q) << delta - threshold_sigmas*stddev_a;
    flagged(p,q) << Meq.MergeFlags(inp,Meq.ZeroFlagger(fc,flag_bit=4,oper="GE"));
  return flagged;

# MS options first
mssel = Meow.Context.mssel = Meow.MSUtils.MSSelector(has_input=True,hanning=True,has_output=False,tile_sizes=[100,200,500],flags=True);
# MS compile-time options
TDLCompileOptions(*mssel.compile_options());
# MS run-time options
TDLRuntimeOptions(*mssel.runtime_options());
## also possible:

TDLCompileOption('flag_xx_yy',"Flag on XX/YY values only",True);
TDLCompileOption('flag_all_corrs',"Merge flags across all correlations",True);
TDLCompileOption('avg_freq',"Average input over frequency",True);
TDLCompileOption('flag_absmax',"Flag on absolute value >=",[None,1.,2.],more=float);
TDLCompileOption('flag_absmin',"Flag on absolute value <=",[None,1.,2.],more=float);
TDLCompileOption('flag_rms',"Flag on rms sigmas >= ",[None,3.,5.,10.],more=float);

def _define_forest(ns):
  ANTENNAS = mssel.get_antenna_set(list(range(1,15)));
  array = Meow.IfrArray(ns,ANTENNAS,mirror_uvw=False);
  observation = Meow.Observation(ns);
  Meow.Context.set(array,observation);

  outputs = spigots = array.spigots(corr=mssel.get_corr_index(),flag_bit=1);
  Meow.Bookmarks.make_node_folder("Input visibilities by baseline",
    [ spigots(p,q) for p,q in array.ifrs() ],sorted=True,ncol=2,nrow=2);

  # extract xx/yy if asked
  if flag_xx_yy:
    outputs = ns.xxyy;
    for p,q in array.ifrs():
      outputs(p,q) << Meq.Selector(spigots(p,q),index=[0,3],multi=True);
      
  # add freq averaging if needed
  if avg_freq:
    for p,q in array.ifrs():
      ns.freqavg(p,q) << Meq.Mean(outputs(p,q),reduction_axes=['freq']);
    outputs = ns.freqavg;
  
  # make an inspector for spigots, we'll add more to this list
  inspectors = [
    Meow.StdTrees.vis_inspector(ns.inspect('spigots'),spigots,bookmark=False)
  ];

  # flag on absolute value first
  if flag_absmax is not None or flag_absmin is not None:
    outputs = abs_clip(outputs,flag_absmax,flag_absmin);
    inspectors.append(Meow.StdTrees.vis_inspector(ns.inspect('abs'),outputs,bookmark=False));

  # then flag on rms
  if flag_rms is not None:
    outputs = rms_clip(outputs,flag_rms);
    inspectors.append(Meow.StdTrees.vis_inspector(ns.inspect('rms'),outputs,bookmark=False));
    
  # recreate 2x2 result (if only flagging on xx/yy)
  if flag_xx_yy:
    for p,q in array.ifrs():
      out = outputs(p,q);
      xx = ns.xx_out(p,q) << Meq.Selector(out,index=0);
      yy = ns.yy_out(p,q) << Meq.Selector(out,index=1);
      ns.make4corr(p,q) << Meq.Matrix22(xx,0,0,yy);
    outputs = ns.make4corr;
      
  # merge flags across correlations if asked
  if flag_all_corrs:
    for p,q in array.ifrs():
      ns.mergecorrflags(p,q) << Meq.MergeFlags(outputs(p,q));
    outputs = ns.mergecorrflags;
    
  # finally, merge flags with spigots (so that we can inspect output flags with original data)
  for p,q in array.ifrs():
    ns.output(p,q) << Meq.MergeFlags(spigots(p,q),outputs(p,q));
  outputs = ns.output;
  inspectors.append(Meow.StdTrees.vis_inspector(ns.inspect('output'),outputs,bookmark=False));
  
  # make sinks and vdm
  Meow.StdTrees.make_sinks(ns,outputs,
                           post=inspectors,output_col='');
  Meow.Bookmarks.make_node_folder("Output visibilities by baseline",
    [ outputs(p,q) for p,q in array.ifrs() ],sorted=True,ncol=2,nrow=2);

  # put all inspectors into bookmarks
  pg = Meow.Bookmarks.Page("Vis Inspectors",2,2);
  for node in inspectors:
    pg.add(node,viewer="Collections Plotter");

  # finally, setup imaging options
  imsel = mssel.imaging_selector(npix=512);
  TDLRuntimeMenu("Imaging options",*imsel.option_list());

def _tdl_job_run_flagging (mqs,parent,**kw):
  req = mssel.create_io_request();
  mqs.execute('VisDataMux',req,wait=False);



if __name__ == '__main__':
    Timba.TDL._dbg.set_verbose(5);
    ns = NodeScope();
    _define_forest(ns);


    ns.Resolve();
    pass

