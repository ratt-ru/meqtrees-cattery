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

from Timba.TDL import *
from Meow import Context

import math

from . import ErrorGens

DEG = math.pi/180.;
ARCMIN = DEG/60;
ARCSEC = DEG/3600;

def compute_pointings(nodes,stations=None,**kw):
  """Computes pointing errors for a list of stations.""";
  # if the error generator is set to "no error", return None
  if not _pe_errgen_l.has_errors() and not _pe_errgen_m.has_errors():
    return None;
  node_maker1 = _pe_errgen_l.node_maker();
  node_maker2 = _pe_errgen_m.node_maker();
  ns = nodes.Subscope();
  stations = stations or Context.array.stations();
  if station_subset != STATIONS_ALL:
    subset = set(re.split("[\s,]+",station_subset));
  else:
    subset = set(stations);
  # create nodes to compute pointing errors per antenna
  for p in stations:
    if p in subset:
      node_maker1(ns.l(p),station=p,axis='l');
      node_maker2(ns.m(p),station=p,axis='m');
      nodes(p) << Meq.Composer(ns.l(p),ns.m(p));
    else:
      nodes(p) << Meq.Composer(0,0);
  return nodes

# make an error generator for pointings
_pe_errgen_l = ErrorGens.Selector("pointing in l",0,10,label="pe_l",unit=("arcsec",ARCSEC));
_pe_errgen_m = ErrorGens.Selector("pointing in m",0,10,label="pe_m",unit=("arcsec",ARCSEC));

STATIONS_ALL = "all";

TDLCompileOption("station_subset","Stations with pointing errors",[STATIONS_ALL],more=str,doc="""<P>
  Specify a list of (space- or comma-separated) station names, or "all".</P>""");

TDLCompileOptions(*_pe_errgen_l.options());
TDLCompileOptions(*_pe_errgen_m.options());

