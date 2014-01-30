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
from Meow import Context

import random
import math
import ErrorGens

DEG = math.pi/180.;
ARCMIN = DEG/60;
ARCSEC = DEG/3600;

def compute_jones (Jones,stations=None,**kw):
  stations = stations or Context.array.stations();
  
  # get error generator function for gain and phase
  gaingen = _gain_errgen.node_maker();
  phasegen = _phase_errgen.node_maker();
  
  for p in stations:
    Jj = Jones(p);
    ns = Jj.Subscope();
    # generate nodes
    gaingen(ns.xa,station=p);
    gaingen(ns.ya,station=p);
    phasegen(ns.xp,station=p);
    phasegen(ns.yp,station=p);
    # put them together into matrix
    Jj << Meq.Matrix22(
      ns.x << Meq.Polar(ns.xa,ns.xp),0,0,ns.y << Meq.Polar(ns.ya,ns.yp)
    );
  return Jones;

_gain_errgen = ErrorGens.Selector("gain",1,[.5,1.5]);
_phase_errgen = ErrorGens.Selector("phase",0,60,unit=("deg",DEG));

TDLCompileOptions(*_gain_errgen.options());
TDLCompileOptions(*_phase_errgen.options());
