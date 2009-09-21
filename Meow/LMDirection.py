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

class LMDirection (Direction):
  """A LMDirection represents a direction on the sky, specified
  in l,m coordinates relative to some other direction dir0 
  (the phase center of the global observation in Meow.Context, by default).
  """;
  def __init__(self,ns,name,l,m,dir0=None,
               static=True,quals=[],kwquals={}):
    Parameterization.__init__(self,ns,name,
                              quals=quals,kwquals=kwquals);
    self._dir0 = Context.get_dir0(dir0);
    self._add_parm('l',l,tags="direction");
    self._add_parm('m',m,tags="direction");
    self.static = static;
    if static and self._is_constant('l') and self._is_constant('m'):
      n = math.sqrt(1-l*l-m*m);
      self._add_parm('n',n,tags="direction");
      self.static = l,m,n;
      self.static_radec = {};
    else:
      self.static = None;
      
  def radec (self):
    """Returns ra-dec two-pack for this direction.""";
    radec = self.ns.radec;
    if not radec.initialized():
      radec << Meq.LMRaDec(radec_0=self._dir0.radec(),lm=self.lm(self._dir0));
    return radec;
 
  def _lm (self,dir0=None):
    """Helper function: creates L,M nodes as needed, returns them as an (l,m) 
    tuple. dir0 is a direction relative to which lm is computed, at the
    moment it is not used.
    """;
    return (self._parm("l"),self._parm("m"));

  def lmn (self,dir0=None):
    """Returns LMN three-pack for this component.
    dir0 is a direction relative to which lm is computed.
    BUG here: reference direction is ignored
    """;
    lmn = self.ns.lmn;
    if not lmn.initialized():
      if self.static:
        lmn << Meq.Constant(value=self.static);
      else:
        l,m = self._lm();
        n = self.ns.n << Meq.Sqrt(1-Meq.Sqr(l)-Meq.Sqr(m));
        lmn << Meq.Composer(l,m,n);
    return lmn;
    
  def lmn_static (self,dir0=None):
    """Returns static LMN tuple, given a reference direction dir0, or using the global phase 
    center if not supplied. 
    Both this direction and the reference direction must be static, otherwise None is returned.
    BUG here: reference direction is ignored
    """;
    return self.static;
    
  def radec_static (self,dir0=None):
    """Returns static RA-Dec tuple, given a reference direction dir0, or using the global phase 
    center if not supplied. 
    Both this direction and the reference direction must be static, otherwise None is returned.
    BUG here: reference direction is ignored
    """;
    dir0 = Context.get_dir0(dir0);
    if not self.static or not dir0.radec_static():
      return None;
    ra0,dec0 = dir0.radec_static();
    l,m,n = self.static;
    # see if we have already computed this lmn
    radec = self.static_radec.get((ra0,dec0),None);
    if radec is None:
      radec = self.static_radec[(ra0,dec0)] = lm_to_radec(l,m,ra0,dec0);
    return radec;
