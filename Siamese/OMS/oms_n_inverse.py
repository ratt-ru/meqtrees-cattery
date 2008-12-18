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

import math

DEG = math.pi/180.;
ARCMIN = DEG/60;
ARCSEC = DEG/3600;
C = 299792458;  # speed of light

def compute_jones (Jones,sources,stations=None,**kw):
  """Computes conjugate N for a list of sources.
  """;
  stations = stations or Context.array.stations();
  ns = Jones.Subscope();
  for p in stations:
    ns.w(p) << Meq.Selector(Context.array.uvw(Context.observation.phase_centre)(p),index=2);
  
  ns.wl << (2*math.pi*Meq.Freq())/C;
  
  for src in sources:
    n1 = ns.n1(src) << ns.wl*(src.direction.n()-1);
    for p in stations:
      Jones(src,p) << Meq.Polar(1.,n1*ns.w(p));
  
  return Jones;

