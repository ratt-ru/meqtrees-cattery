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

from Timba.TDL import *
import math
from Meow import Jones
from Meow import Context

H = 300000;           # height of ionospheric layer, in meters
Lightspeed = 3e+8;

## forget about this option for now: were trying to do a more realistic coordinate
## system, so no point in rotating ionospheres
# TDLCompileOption("iono_rotate","Rotate ionosphere with sky",True);

def compute_piercings (ns,source_list,stations=None):
  """Creates nodes to compute the "piercing points" of each
  source in source_list, for each antenna in array.""";
  stations = stations or Context.array.stations();
  xyz = Context.array.xyz();
  
  for p in stations:
    # project our stations into an ionospheric coordinate system xyz1
    # for now, we just project everything to z=0, and take station 1 as origin
    # xyz_proj gives positions projected to equatorial plane (z=0)
    xyzp = ns.xyz_proj(p) << Meq.Paster(xyz(p),0,index=2);
    # xyz1 gives positions relative to station1
    xyz1 = ns.xyz1(p) << xyzp - ns.xyz_proj(stations[0]);
    # and extract the xy1 component
    ns.xy1(p) << Meq.Selector(xyz1,index=[0,1],multi=True);

    for src in source_list:
      # now get the source azimuth and ZA as seen from station p
      az = src.direction.az(xyz(p));
      za = math.pi/2 - src.direction.el(xyz(p));
      # pR: distance to pierce point (from point directly overhead of station)
      # is height*tg ZA
      r = ns.pR(src,p) << H*Meq.Tan(za);
      # dx/dy: relative xy of pierce point. This is given by azimuth (0 is north.)
      dxy =  ns.dxy(src,p) << Meq.Composer( ns.dx(src,p) << r*Meq.Sin(az),
                                                                   ns.dy(src,p) << r*Meq.Cos(az));
      # pxy are the absolute piercing coordinatesthen
      ns.pxy(src.name,p) << ns.xy1(p) + dxy;
  return ns.pxy;
  
def compute_za_cosines (ns,source_list,stations=None):
  """Creates nodes to compute the zenith angle cosine of each
  source in source_list, for each antenna in array.""";
  stations = stations or Context.array.stations();
  xyz = Context.array.xyz();
  za_cos = ns.za_cos;
  for src in source_list:
    for p in stations:
      za_cos(src,p) << Meq.Cos(math.pi/2 - src.direction.el(xyz(p)));
      
  return za_cos;

def compute_zeta_jones_from_tecs (zeta,tecs,source_list,stations):
  """Creates the Z Jones for ionospheric phase, given TECs (per source, 
  per station).""";
  stations = stations or Context.array.stations();
  print [src.name for src in source_list];
  print stations;
  for src in source_list:
    for p in stations:
      zeta(src,p) << Meq.Polar(1,-25*Lightspeed*tecs(src,p)/Meq.Freq());
  return zeta;

