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

"""<P>This module implements voltage beam patterns that are read in from FITS files and interpolated
via a PyNode. This implementation allows for a lot of flexibility in how the beam  is
interpolated.</P>

<P align="right">Author: O. Smirnov &lt;<tt>smirnov@astron.nl</tt>&gt;</P>""";

from Timba.TDL import *
import Meow
from Meow import Context
from Meow import StdTrees


import os
import os.path

TDLCompileOption("filename_pattern","Filename pattern",["beam_$(xy)_$(reim).fits"],more=str,doc="""<P>
  Pattern for beam FITS filenames. The real and imaginary parts of each element of the beam Jones matrix should be stored in a separate FITS file. A number of variables will be substituted in the beam pattern, these may be introduced as
  "$var" or "$(var)". Use "$$" for a "$" character. The following variables are recognized:</P>
  <UL>
  <LI><B>$xy</B> or <B>$corr</B>: the correlation index, in lowercase 
  (xx/xy/yx/yy or rr/rl/lr/ll)</LI>
  <LI><B>$XY</B> or <B>$CORR</B>: the correlation index, in uppercase 
  (XX/XY/YX/YY or RR/RL/LR/LL)</LI>
  <LI><B>$reim</B> or <B>$ReIm</B> or <B>$REIM</B>: "re", "Re" or "RE" for real, "im", "Im" or "IM" for imaginary</LI>
  <LI><B>$realimag</B> or <B>$RealImag</B> or <B>$REALIMAG</B>: "real", "Real" or "REAL" for real, "imag", "Imag" or "IMAG" for imaginary</LI>
  <UL>""");
TDLCompileOption("missing_is_null","Use null for missing beam elements",True,doc="""<P>
  If True, then 0 will be used when a beam file is missing. Useful if you only have beams for
  the diagonal elements (XX/YY), or only real beams.
  </P>""");
TDLCompileOption("spline_order","Spline order for interpolation",[1,2,3,4,5],default=3);
TDLCompileOption("normalize_gains","Normalize max beam gain to 1",False);
TDLCompileOption("sky_rotation","Include sky rotation",True,doc="""<P>
  If True, then the beam will rotate on the sky with parallactic angle. Use for e.g. alt-az mounts.)
  </P>""");

CORRS = Context.correlations;
REIM = "re","im";
REALIMAG = dict(re="real",im="imag");

import Utils

def make_beam_filename (filename_pattern,corr,reim):
  """Makes beam filename for the given correlation and real/imaginary component (one of "re" or "im")"""
  return Utils.substitute_pattern(filename_pattern,
    corr=corr.lower(),xy=corr.lower(),CORR=corr.upper(),XY=corr.upper(),
    reim=reim.lower(),REIM=reim.upper(),ReIm=reim.title(),
    realimag=REALIMAG[reim].lower(),REALIMAG=REALIMAG[reim].upper(),
    RealImag=REALIMAG[reim].title());

def make_beam_node (beam,pattern,*children):
  """Makes beam interpolator node for the given filename pattern.""";
  filename_real = [];
  filename_imag = [];
  for corr in CORRS:
    # make FITS images or nulls for real and imaginary part
    filename_real.append(make_beam_filename(pattern,corr,'re'));
    filename_imag.append(make_beam_filename(pattern,corr,'im'));
  # if we end up with only one unique filename, make a scalar-mode interpolator
  if len(set(filename_real)) == 1 and len(set(filename_imag)) == 1:
    filename_real = filename_real[0];
    filename_imag = filename_imag[0];
  # now make interpolator node
  import InterpolatedBeams
  beam << Meq.PyNode(class_name="FITSBeamInterpolatorNode",module_name=InterpolatedBeams.__file__,
                     filename_real=filename_real,filename_imag=filename_imag,normalize=normalize_gains,
                     missing_is_null=missing_is_null,spline_order=spline_order,verbose=0,
                     children=children);

def compute_jones (Jones,sources,stations=None,pointing_offsets=None,inspectors=[],label='E',**kw):
  stations = stations or Context.array.stations;
  ns = Jones.Subscope();

  # declare an inspector node for the Jones matrix -- will be defined below
  insp = Jones.scope.inspector(label);
  inspectors += [ insp ];

  # if sky rotation is in effect, make rotation matrices
  if sky_rotation:
    for p in stations:
      pa = ns.pa(p) << Meq.ParAngle(radec=Context.observation.radec0(),xyz=Context.array.xyz(p));
      ca = ns.cos_pa(p) << Meq.Cos(pa);
      sa = ns.sin_pa(p) << Meq.Sin(pa);
      ns.parot(p) << Meq.Matrix22(ca,-sa,sa,ca);

  # loop over sources
  for src in sources:
    # If sky rotation and/or pointing offsets are in effect, we have a per-station beam.
    # Otherwise the beam is the same for all stations.
    if sky_rotation or pointing_offsets:
      for p in stations:
        lm = src.direction.lm();
        # apply rotation to put sources into the antenna frame
        if sky_rotation:
          lm = ns.lmrot(src,p) << Meq.MatrixMultiply(ns.parot(p),src.direction.lm());
        # apply offset (so pointing offsets are interpreted in the azel frame, if rotating)
        if pointing_offsets:
          lm = ns.lmoff(src,p) << lm + pointing_offsets(p);
        # now make the beam node
        make_beam_node(Jones(src,p),filename_pattern,lm);
    else:
      make_beam_node(Jones(src),filename_pattern,src.direction.lm());
      for p in stations:
        Jones(src,p) << Meq.Identity(Jones(src));

  # define an inspector
  if sky_rotation or pointing_offsets:
    # Jones inspector is per-source, per-station
    insp << StdTrees.define_inspector(Jones,sources,stations,label=label);
  else:
    # Jones inspector is per-source
    insp << StdTrees.define_inspector(Jones,sources,label=label);

  return Jones;
