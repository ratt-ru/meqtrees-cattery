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

"""<P>This is the old-style WSRT cos^3 beam model. This has been superseded by
Siamese.OMS.analytic_beams and Calico.OMS.wsrt_beams.</P> 

<P>Author: O. Smirnov.</P>""";
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from Timba.TDL import *
from Meow import Context

import random
import math

from . import ErrorGens


DEG = math.pi/180.;
ARCMIN = DEG/60;
ARCSEC = DEG/3600;

def WSRT_cos3_beam (E,lm,*dum):
  """computes a gaussian beam for the given direction
  'E' is output node
  'lm' is direction (2-vector node)
  """
  ns = E.Subscope();
  # According to the WSRT Observer's Guide
  # (http://www.astron.nl/wsrt/wsrtGuide/node6.html#SECTION00067000000000000000)
  # the full power WSRT beam can be modelled roughly as
  #   cos^6(c x freq x r)
  # where freq is in GHz (hence the 1e-9 factor below), r is distance
  # from pointing center (in same units as expected by the cos() function, i.e.
  # in degrees if cos() is taken as a function of degrees), and "the constant c=68 is,
  # to first order, wavelength independent at GHz frequencies (declining to  c=66 at
  # 325 MHz and c = 63 at 4995 MHz)."
  # Hence the E-Jones factor is simply cos^3.
  #
  ns.lmsq << Meq.Sqr(lm);
  ns.lsq  << Meq.Selector(ns.lmsq,index=0);
  ns.msq  << Meq.Selector(ns.lmsq,index=1);
  ns.Ex << Meq.Pow(Meq.Cos(Meq.Sqrt(ns.lsq*(1+beam_ellipticity)+ns.msq*(1-beam_ellipticity))*(wsrt_beam_size_factor*1e-9)*Meq.Freq()),3);
  ns.Ey << Meq.Pow(Meq.Cos(Meq.Sqrt(ns.lsq*(1-beam_ellipticity)+ns.msq*(1+beam_ellipticity))*(wsrt_beam_size_factor*1e-9)*Meq.Freq()),3);
  E << Meq.Matrix22(ns.Ex,0,0,ns.Ey);
  return E;

# this beam model is not per-station
WSRT_cos3_beam._not_per_station = True;

def compute_jones (Jones,sources,stations=None,pointing_offsets=None,**kw):
  """Computes beam gain for a list of sources.
  The output node, will be qualified with either a source only, or a source/station pair
  """;
  stations = stations or Context.array.stations();
  ns = Jones.Subscope();
  # are pointing errors configured?
  if pointing_offsets:
    # create nodes to compute actual pointing per source, per antenna
    for p in stations:
      for src in sources:
        lm = ns.lm(src.direction,p) << src.direction.lm() + pointing_offsets(p);
        beam_model(Jones(src,p),lm,p);
  # no pointing errors
  else:
    # if not per-station, use same beam for every source
    if beam_model._not_per_station:
      for src in sources:
        beam_model(Jones(src),src.direction.lm());
        for p in Context.array.stations():
          Jones(src,p) << Meq.Identity(Jones(src));
    else:
      for src in sources:
        for p in Context.array.stations():
          beam_model(Jones(src,p),src.direction.lm(),p);
  return Jones;

_model_option = TDLCompileOption('beam_model',"Beam model",
  [WSRT_cos3_beam]
);

_wsrt_option_menu = TDLCompileMenu('WSRT beam model options',
  TDLOption('wsrt_beam_size_factor',"Beam size factor",[68.],more=float),
  TDLOption('beam_ellipticity',"Beam ellipticity",[.01],more=float)
);

def _show_option_menus (model):
  _wsrt_option_menu.show(model==WSRT_cos3_beam);

_model_option.when_changed(_show_option_menus);
