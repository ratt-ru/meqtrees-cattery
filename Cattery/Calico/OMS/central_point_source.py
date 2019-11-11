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

"""This is a Sky Model module.
This implements a sky model with a single central point source.
""";

from Timba.TDL import *
import Meow

TDLCompileOption("spectral_index","Spectral index",[None,0.],more=float,
    doc="""Spectral index of source, use None for a flat spectrum""");
TDLCompileOption("ref_freq","SpI reference frequency, MHz",1421,more=float);

options = [];

def runtime_options ():
  return options;

def source_list (ns,name="cps"):
  # figure out spectral index parameters
  if spectral_index is not None:
    spi_def = Meow.Parm(spectral_index);
    freq0_def  = ref_freq;
  else:
    spi_def = freq0_def = None;
  i_def = Meow.Parm(1);
  quv_def = Meow.Parm(0);

  srcdir = Meow.LMDirection(ns,name,0,0);
  src = Meow.PointSource(ns,name,srcdir,I=i_def,Q=quv_def,U=quv_def,V=quv_def,spi=spi_def,freq0=freq0_def);

  ## define a parmgroup for source parameters
  ## now make a solvejobs for the source
  #pg_src = ParmGroup("source",
              #src.coherency().search(tags="solvable"),
              #table_name="sources.mep",
              #individual=True,
              #bookmark=True);
  ## now make a solvejobs for the source
  #options.append(pg_src.make_solvejob_menu("Calibrate source fluxes"));

  return [ src ];

