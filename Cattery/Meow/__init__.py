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
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

import sys

_meow_path = __path__[0];

# various Meow options that can be enabled/disabled by imported modules
DiscoverMaximumW = False;

# do nothing if we're imported on the kernel-side, as reloading
# the octopython stuff confuses the hell out of the kernel
if "Timba.meqkernel" in sys.modules:
  pass;
  
else:
  import __main__
  _verbosity = getattr(__main__,'_meow_verbosity',1);

  def dprint (msg,level=1):
    if level <= _verbosity:
      print(msg);

  from Cattery.Meow.CorruptComponent import CorruptComponent
  from Cattery.Meow.Position import Position
  from Cattery.Meow.Direction import Direction
  from Cattery.Meow.LMDirection import LMDirection
  from Cattery.Meow.AzElDirection import AzElDirection
  from Cattery.Meow.LMApproxDirection import LMApproxDirection
  from Cattery.Meow.IfrArray import IfrArray
  from Cattery.Meow.Observation import Observation
  from Cattery.Meow.Parameterization import Parameterization
  from Cattery.Meow.Patch import Patch
  from Cattery.Meow.SkyComponent import SkyComponent
  from Cattery.Meow.PointSource import PointSource
  from Cattery.Meow.KnownVisComponent import KnownVisComponent
  from Cattery.Meow.GaussianSource import GaussianSource
  from Cattery.Meow.SixpackComponent import SixpackComponent
  from Cattery.Meow.FITSImageComponent import FITSImageComponent
  from Cattery.Meow.Parm import Parm
  from Cattery.Meow import StdTrees
  from Cattery.Meow import Bookmarks
  from Cattery.Meow import MSUtils
  from Cattery.Meow import Utils

  __all__ = [
              CorruptComponent,
              Position,
              Direction,
              LMDirection,
              LMApproxDirection,
              FITSImageComponent,
              GaussianSource,
              IfrArray,
              Observation,
              Parameterization,
              Patch,
              PointSource,
              SixpackComponent,
              SkyComponent,
              Bookmarks,
              MSUtils,
              Utils,
              _meow_path,
              StdTrees,
  ];
