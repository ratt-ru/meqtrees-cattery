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

from Timba.TDL import *
from Timba.Meq import meq
from Timba.array import *
import os
import random

import Meow

from Meow import Bookmarks,Context
import Meow.StdTrees

# MS options first
mssel = Context.mssel = Meow.MSUtils.MSSelector(has_input=True,read_flags=True,has_output=False,tile_sizes=[10,100,200]);
# MS compile-time options
TDLCompileOptions(*mssel.compile_options());
# MS run-time options
TDLRuntimeMenu("Data selection & flag handling",*mssel.runtime_options());

def _define_forest(ns):
  # setup contexts from MS
  mssel.setup_observation_context(ns);
  array = Meow.Context.array;
  observation = Meow.Context.observation;

  stas = array.stations();

  # make spigot nodes
  spigots = spigots0 = array.spigots(corr=mssel.get_corr_index());

  for p,q in array.ifrs():
    spigots('abs',p,q) << Meq.Abs(spigots(p,q));

  # ...and an inspector for them
  Meow.StdTrees.vis_inspector(ns.inspector('input'),spigots,vells_label=Context.correlations,
                              bookmark="Inspect input visibilities");
  Meow.StdTrees.vis_inspector(ns.inspector('ampl'),spigots('abs'),vells_label=Context.correlations,
                              bookmark="Inspect mean visibility amplitudes");
  Bookmarks.make_node_folder("Input visibilities by baseline",
    [ spigots(p,q) for p,q in array.ifrs() ],sorted=True,ncol=2,nrow=2);

  ns.inspectors << Meq.ReqMux(ns.inspector('input'),ns.inspector('ampl'));

  ns.VisDataMux << Meq.VisDataMux(post=ns.inspectors);

  # add imaging options
  imsel = mssel.imaging_selector(npix=512,arcmin=120);
  TDLRuntimeMenu("Make an image from this MS",*imsel.option_list());

def _tdl_job_View_MS (mqs,parent,**kw):
  req = mssel.create_io_request();
  mqs.execute('VisDataMux',req,wait=False);
