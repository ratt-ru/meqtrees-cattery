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
"""This is a Pointing module implementing solvable pointing errors.
""";

from Timba.TDL import *
from Meow import Context
from Meow import ParmGroup

import math

DEG = math.pi/180.;
ARCMIN = DEG/60;
ARCSEC = DEG/3600;

def compute_pointings (nodes,stations=None,label="pnt",return_parms=None,**kw):
  """Computes pointing errors for a list of stations.""";
  # if the error generator is set to "no error", return None
  ns = nodes.Subscope();
  # create parms for pointing errors per antenna
  pps = [];
  for p in stations or Context.array.stations():
    dl = ns.dl(p) << Meq.Parm(0,tags="pointing solvable");
    dm = ns.dm(p) << Meq.Parm(0,tags="pointing solvable");
    nodes(p) << Meq.Composer(dl,dm);
    pps += [dl,dm];
  if return_parms is None:
    # add parmgroup
    global pg_pointing;
    pg_pointing = ParmGroup.ParmGroup(label,pps,table_name="%s.fmep"%label,bookmark=True);
    ParmGroup.SolveJob("cal_"+label,"Calibrate pointing errors",pg_pointing);
  else:
    return_parms += pps;

  return nodes
