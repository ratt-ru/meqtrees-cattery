# -*- coding: utf-8 -*-
# file ../Meow/Position.py

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
from Parameterization import *


class Position (Parameterization):
  """A Position represents a location in 1,2,3,4 dimensional space.
  Its dimensionality depends on the coordinates (x,y,z,t) that are
  specified at creation. These may be numbers or nodes.
  """;
  def __init__(self, ns, name,
               x=None, y=None, z=None, t=None,
               quals=[],kwquals={}):
    Parameterization.__init__(self,ns,name,
                              quals=quals,kwquals=kwquals)
    self.dims = [0,0,0,0]       # i.e. [x,y,z,t]
    if not y==None:
      self.dims[0] = 1
      self._add_parm('xpos', x, tags=['position','xpos'])
    if not y==None:
      self.dims[1] = 1
      self._add_parm('ypos', y, tags=['position','ypos'])
    if not z==None:
      self.dims[2] = 1
      self._add_parm('zpos', z, tags=['position','zpos'])
    if not t==None:
      self.dims[3] = 1
      self._add_parm('tpos', t, tags=['position','tpos'])


  def xyz (self):
    """Returns the xyt 3-vector node for this position."""
    if self.dims[0]*self.dims[1]*self.dims[2]==0:
      raise ValueError,'Position does not have enough dimensions' 
    xyz = self.ns.xyz
    if not xyz.initialized():
      x = self._parm('xpos')
      y = self._parm('ypos')
      z = self._parm('zpos')
      # print str(x),str(y),str(z)
      xyz << Meq.Composer(x,y,z)
    return xyz


#================================================================
    
if __name__ == '__main__':
    ns = NodeScope()
    p1 = Position(ns,'p1', x=1, y=2, z=(ns << Meq.Parm(6.5)))
    xyz = p1.xyz()
    print str(xyz)
