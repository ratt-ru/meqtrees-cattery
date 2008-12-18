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
"""This implements a Jones module for WSRT beams with pointing errors.
THis code is untested and probably won't compile.
""";

from Timba.TDL import *
from Meow import Context
from Meow import StdTrees
from Meow import ParmGroup

import random
import math
from math import sqrt,atan2

DEG = math.pi/180.;
ARCMIN = DEG/60;
ARCSEC = DEG/3600;

def WSRT_cos3_beamdiff (dE,lm,dlm,freqscale,*dum):
  """computes a derivative of a WSRT cos^3 beam for a given lm,dlm.
  lm is a tuple of l,m coordinates
  dl,dm are poitning offsets.
  """
  l,m = lm;
  r = ns.r << atan2(sqrt((1-m**2)*l*l+m**2),sqrt(1-m**2)*sqrt(1-l**2));
  ns.cosr << Meq.Cos(ns.rclip<<Meq.Min(ns.rscale<<freqscale*ns.r,cutoff);
  denom = dE('denom') << 3*Meq.Sqr(ns.cosr)*(-Meq.Sin(ns.rclip))*ns.freqscale/ns.r;
  dEdlm = dE('dlm') << dlm/denom;
  dE << Meq.Norm(dEdlm);
  return dE;

def compute_jones (Jones,sources,stations=None,label="beam",pointing_offsets=None,inspectors=[],**kw):
  """Computes beam gain for a list of sources.
  The output node, will be qualified with either a source only, or a source/station pair
  """;
  stations = stations or Context.array.stations();
  ns = Jones.Subscope();
  beamscale = ns.beamscale << Meq.Parm(wsrt_beam_size_factor,tags="beam solvable");
  freqscale = ns.freqscale << (beamscale*1e-9)*Meq.Freq();
  # add solution for beamscale
  global pg_beam;
  pg_beam = ParmGroup.ParmGroup(label,[beamscale],table_name="%s.fmep"%label,bookmark=True);
  ParmGroup.SolveJob("cal_"+label+"_scale","Calibrate beam scale",pg_beam);
  
  # are pointing errors configured?
  if pointing_offsets:
    # create nodes to compute actual pointing per source, per antenna
    for p in Context.array.stations():
      for src in sources:
        lm = ns.lm(src.direction,p) << src.direction.lm() + pointing_offsets(p);
        beam_model(Jones(src,p),lm,freqscale,p);
  # no pointing errors
  else:
    # if not per-station, use same beam for every source
    if beam_model._not_per_station:
      for src in sources:
        lmnst = src.direction.lmn_static();
        if lmnst:
          beam_model(Jones(src),lmnst[0:2],freqscale);
        else:
          beam_model(Jones(src),src.direction.lm(),freqscale);
      inspectors.append(Jones('inspector') << StdTrees.define_inspector(Jones,sources));
      return SourceBeam(Jones);
    else:
      for src in sources:
        for p in Context.array.stations():
          beam_model(Jones(src,p),src.direction.lm(),freqscale,p);
  return Jones;

_model_option = TDLCompileOption('beam_model',"Beam model",
  [WSRT_cos3_beam]
);

_wsrt_option_menu = TDLCompileMenu('WSRT beam model options',
  TDLOption('wsrt_beam_size_factor',"Beam size factor (1/GHz)",[64.],more=float),
  TDLOption('wsrt_beam_cutoff',"Cosine argument cutoff (degrees)",[100.],more=float)
);

def _show_option_menus (model):
  _wsrt_option_menu.show(model==WSRT_cos3_beam);

_model_option.when_changed(_show_option_menus);

def runtime_options ():
  return [];