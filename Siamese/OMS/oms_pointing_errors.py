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

import ErrorGens

DEG = math.pi/180.;
ARCMIN = DEG/60;
ARCSEC = DEG/3600;

def compute_pointings(nodes,stations=None,**kw):
  """Computes pointing errors for a list of stations.""";
  # if the error generator is set to "no error", return None
  if not _pe_errgen.has_errors():
    return None;
  node_maker = _pe_errgen.node_maker();
  ns = nodes.Subscope();
  # create nodes to compute pointing errors per antenna
  for p in stations or Context.array.stations():
    node_maker(ns.l(p),station=p);
    node_maker(ns.m(p),station=p);
    nodes(p) << Meq.Composer(ns.l(p),ns.m(p));
  return nodes

# make an error generator for pointings
_pe_errgen = ErrorGens.Selector("pointing",0,10,label="pe",unit=("arcsec",ARCSEC));

TDLCompileOptions(*_pe_errgen.options());

