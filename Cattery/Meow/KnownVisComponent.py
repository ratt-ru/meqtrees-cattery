# -*- coding: utf-8 -*-
#
#% $Id: SkyComponent.py 6583 2008-11-24 09:27:58Z oms $
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
from .SkyComponent import *
from .Direction import *
from . import Context

class KnownVisComponent (SkyComponent):
  """A KnownVisComponent is simply a SkyComponent whose visibilities are
  provided by the caller, via the 'visnodes' argument (which should be qualified
  by a station pair.)
  """;
  def __init__(self,ns,name,visnodes,direction=None):
    SkyComponent.__init__(self,ns,name,direction or Context.observation.phase_center);
    self.visnodes = visnodes;

  def visibilities (self,array=None,observation=None,nodes=None,**kw):
    return self.visnodes;

