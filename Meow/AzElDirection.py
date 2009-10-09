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
from Timba.Meq import meq
from Direction import Direction,lm_to_radec
from Parameterization import Parameterization
import Context
import math

class AzElDirection (Direction):
  """A AzElDirection represents a direction on the sky, specified
  in az,el coordinates relative to some xyz position. Note that xyz should be a 3-vector node.
  If None, then Context.array.xyz0() is used.
  """;
  def __init__(self,ns,name,az,el,xyz=None,
               quals=[],kwquals={}):
    Parameterization.__init__(self,ns,name,
                              quals=quals,kwquals=kwquals);
    self._xyz = xyz or Context.array.xyz0();
    self._add_parm('az',az,tags="direction");
    self._add_parm('el',el,tags="direction");
    self.static = None;

  def azel (self,xyz=None):
    """Returns az-el 2-vector node for this direction.
    BUG: xyz will be ignored if passed in. Should really implement azel->azel conversion
    (at different coordinates) in meqtrees.""";
    azel = self.ns.azel;
    if not azel.initialized():
      azel << Meq.Composer(self._parm('az'),self._parm('el'));
    return azel;
     
  def radec (self):
    """Returns ra-dec two-pack for this direction.""";
    radec = self.ns.radec;
    if not radec.initialized():
      radec << Meq.AzElRaDec(azel=self.azel(),xyz=self._xyz);
    return radec;
 
