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
"""This module implements an analytic VLA beam model given by 
  Uson & Cotton (2008) "Beam squint and Stokes V with off-axis feeds"
  Bibliographic code 2008A&A...486..647U.

Authors: Fred Dulwich and Shannon Jaeger"""

from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from Timba.TDL import *
from Meow import Context

import random
import math


DEG = math.pi/180.;
ARCMIN = DEG/60;
ARCSEC = DEG/3600;

def VLA_beam_squint (E,lm,p):
  """
  Computes a beam for the given direction.
  'E' is output node.
  'lm' is direction (2-vector node).
  'p' is the current station.

  This routine computes the VLA beam using the model given in
  Uson & Cotton (2008) "Beam squint and Stokes V with off-axis feeds"
  Bibliographic code 2008A&A...486..647U.
  """
  ns = E.Subscope();

  # Beam coefficients c1 to c6 from Uson & Cotton (2008).
  cn = [-0.56249985, 0.21093573, -0.03954289, 0.00443319, -0.00031761, 0.00001109];

  # If not using parallactic angle, compute the (constant)
  # sine and cosine of the feed angle for the dx and dy terms.
  if not vla_pa:
    ca = math.cos(vla_feed_angle * DEG);
    sa = math.sin(vla_feed_angle * DEG);
  else:
    # Get the parallactic angle (depends on antenna number).
    radec = Context.observation.radec0();
    xyz = Context.array.xyz(p);
    ns.pa << Meq.ParAngle(radec=radec, xyz=xyz);

    # Compute the sine and cosine of (feedAngle - parAngle) as MeqNodes.
    ca = ns.cos_angle << Meq.Cos(vla_feed_angle * DEG - ns.pa);
    sa = ns.sin_angle << Meq.Sin(vla_feed_angle * DEG - ns.pa);
  # endif

  # Find the beam squint (depends on wavelength).
  ns.squint << ((vla_squint / 3600.) * DEG * 2.997925e+8) / Meq.Freq();

  # Define k and construct the dx and dy terms, taking squint into account.
  k = 1.496e-9 * (25.0 / d_ant);
  ns.lm0 << ns.squint * Meq.Composer(-sa, ca);

  # The E_x term.
  ns.lmsq_x << Meq.Sqr((lm - ns.lm0) / DEG);
  ns.lsq_x  << Meq.Selector(ns.lmsq_x, index=0);
  ns.msq_x  << Meq.Selector(ns.lmsq_x, index=1);
  ns.u_x << Meq.Sqr(k * Meq.Freq()) * (ns.lsq_x + ns.msq_x);
  ns.Ex << 2 * Meq.Add(0.5, *[c*Meq.Pow(ns.u_x, i+1) for i,c in enumerate(cn)])

  # The E_y term.
  ns.lmsq_y << Meq.Sqr((lm + ns.lm0) / DEG);
  ns.lsq_y  << Meq.Selector(ns.lmsq_y, index=0);
  ns.msq_y  << Meq.Selector(ns.lmsq_y, index=1);
  ns.u_y << Meq.Sqr(k * Meq.Freq()) * (ns.lsq_y + ns.msq_y);
  ns.Ey << 2 * Meq.Add(0.5, *[c*Meq.Pow(ns.u_y, i+1) for i,c in enumerate(cn)])

  # Return the beam.
  # Commented out lines produce the difference image (view this the hard way via the tree).
  #  ns.diff << ns.Ex - ns.Ey;
  #  E << Meq.Matrix22(ns.Ey + ns.diff, 0, 0, -ns.Ex-ns.diff);
  E << Meq.Matrix22(ns.Ex, 0, 0, ns.Ey);
  return E;


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
      for src in sources:
        for p in stations:
          beam_model(Jones(src,p),src.direction.lm(),p);
  return Jones;

_model_option = TDLCompileOption('beam_model',"Beam model",
  [VLA_beam_squint]
);

_vla_option_menu = TDLCompileMenu('(E)VLA beam model options',
  TDLOption('vla_squint', "(E)VLA squint factor, arcsec/m", [237.56], more=float,
    doc="This is wavelength dependent:<br>squint&nbsp;=&nbsp;+/-&nbsp;squintFactor&nbsp;*&nbsp;lambda&nbsp;(arcsec/m)."),
  TDLOption('vla_feed_angle', "(E)VLA feed position angle, degrees", [0], more=float,
    doc="The orientation of the feed on the feed circle."),
  TDLOption('d_ant',"(E)VLA effective antenna diameter, metres", [24.5], more=float,
    doc="The effective antenna diameter (24.5 metres at 1.4 GHz)."),
  TDLOption('vla_pa', "Apply parallactic angle", True,
    doc="If set, the parallactic angle is used for the beam squint calculation.")
);

def _show_option_menus (model):
  _vla_option_menu.show(model==VLA_beam_squint);

_model_option.when_changed(_show_option_menus);
