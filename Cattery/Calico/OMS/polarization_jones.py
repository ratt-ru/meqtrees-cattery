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
"""This implements Jones modules related to polarization.
  DecoupledLeakage
  CoupledLeakage
This has not been actively used, since for polarization it is easier to use
a solvable_jones.FullRealImag for the bandpass, and solve for all coefficients.
""";

from Timba.TDL import *
import Meow
from Meow import Context
from Meow import Jones
from Meow import ParmGroup
from Meow import StdTrees

class DecoupledLeakage (object):
  def __init__ (self):
    self.options = [];

  def runtime_options (self):
    return self.options;

  def compute_jones (self,nodes,stations=None,tags=None,label='',inspectors=[],**kw):
    stations = stations or Context.array.stations();
    rot_def = Meow.Parm(0);
    nr = Jones.decoupled_rotation_matrix(nodes("R"),rot_def,tags=tags,series=stations);
    ne = Jones.decoupled_ellipticity_matrix(nodes("E"),rot_def,tags=tags,series=stations);
    
    for p in stations:
      nodes(p) << Meq.MatrixMultiply(nr(p),ne(p));

    # make parmgroups for phases and gains
    self.pg_rot = ParmGroup.ParmGroup(label+"_leakage",
                    nodes.search(tags="solvable"),
                    table_name="%s_leakage.mep"%label,bookmark=True);

    # make solvejobs
    ParmGroup.SolveJob("cal_"+label+"_leakage","Calibrate %s (leakage)"%label,self.pg_rot);
    
    # make inspector for parameters
    StdTrees.inspector(nodes('inspector'),self.pg_rot.nodes,bookmark=False);
    inspectors.append(nodes('inspector'));

    return nodes;

class CoupledLeakage (object):
  def __init__ (self):
    self.options = [];

  def runtime_options (self):
    return self.options;

  def compute_jones (self,nodes,stations=None,tags=None,label='',inspectors=[],**kw):
    stations = stations or Context.array.stations();
    rot_def = Meow.Parm(0);
    nr = Jones.rotation_matrix(nodes("R"),rot_def,tags=tags,series=stations);
    ne = Jones.ellipticity_matrix(nodes("E"),rot_def,tags=tags,series=stations);
    
    for p in stations:
      nodes(p) << Meq.MatrixMultiply(nr(p),ne(p));

    # make parmgroups for phases and gains
    self.pg_rot = ParmGroup.ParmGroup(label+"_leakage",
                    nodes.search(tags="solvable"),
                    table_name="%s_leakage.mep"%label,bookmark=True);

    # make solvejobs
    ParmGroup.SolveJob("cal_"+label+"_leakage","Calibrate %s (leakage)"%label,self.pg_rot);
    
    # make inspector for parameters
    StdTrees.inspector(nodes('inspector'),self.pg_rot.nodes,bookmark=False);
    inspectors.append(nodes('inspector'));

    return nodes;
