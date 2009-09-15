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

"""This module implements a simple dipole projection matrix for an NS-EW pair of
dipoles that is flat on the ground. The matrix looks like
<TABLE>
<TR><TD>cos(az)</TD><TD>-sin(az)*sin(el)</TD></TR>
<TR><TD>sin(az)</TD><TD>cos(az)*sin(el)</TD></TR>
</TABLE>

Author: O. Smirnov &lt;<tt>smirnov@astron.nl</tt>&gt;""";

__default_label__ = "L";
__default_name__  = "dipole projection matrix";

from Timba.TDL import *
from Meow import Context
from Meow import StdTrees

def proj_matrix (src,xyz):
  az,el = src.direction.az(xyz),src.direction.el(xyz);
  azc = az('cos') << Meq.Cos(az);
  azs = az('sin') << Meq.Sin(az);
  els = el('sin') << Meq.Sin(el);
  return Meq.Matrix22(azc,-azs*els,azs,azc*els);

def compute_jones (Jones,sources,stations=None,inspectors=[],label='L',**kw):
  """Creates the dipole projection matrix.""";
  stations = stations or Context.array.stations;
  ns = Jones.Subscope();
  insp = Jones.scope.inspector(label);
  if per_source:
    if per_station:
      for src in sources:
	for p in stations:
	  Jones(src,p) << proj_matrix(src,Context.array.xyz(p));
      insp << StdTrees.define_inspector(Jones,sources,stations,label=label);
    else:
      for src in sources:
	Jones(src) << proj_matrix(src,Context.array.xyz0());
	for p in stations:
	  Jones(src,p) << Meq.Identity(Jones(src));
      insp << StdTrees.define_inspector(Jones,sources,label=label);
  else: # not per source
    if per_station:
      for p in stations:
	Jones(p) << proj_matrix(sources[0],Context.array.xyz(p));
	for src in sources:
	  Jones(src,p) << Meq.Identity(Jones(p));
      insp << StdTrees.define_inspector(Jones,stations,label=label);
    else:
      Jones << proj_matrix(sources[0],Context.array.xyz0());
      for src in sources:
	for p in stations:
	  Jones(src,p) << Meq.Identity(Jones);
      insp << StdTrees.define_inspector(Jones,label=label);

  # add inspectors
  StdTrees.inspector(Jones.scope.inspector(label,'AzEl') ,[src.direction.azel() for src in sources],bookmark=False);
  inspectors += [ insp,Jones.scope.inspector(label,'AzEl') ];

  return Jones;

TDLCompileOption('per_station',"Per-station projection (for large arrays)",True);
TDLCompileOption('per_source',"Per-source projection (for large fields)",True);


