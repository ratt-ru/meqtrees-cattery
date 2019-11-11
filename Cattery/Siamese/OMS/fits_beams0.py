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

"""<P>This module implements voltage beam patterns that are read in from FITS files. This implementation uses the Resampler/Compounder combo, which scales poorly with tile size.
It is recommended that you use small tile sizes (8~16) with this module.</P>

<P align="right">Authors: T. Willis &lt;<tt>Tony.Willis@nrc-cnrc.gc.ca</tt>&gt; & O. Smirnov &lt;<tt>smirnov@astron.nl</tt>&gt;</P>""";
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division
from Timba.TDL import *
import Meow
from Meow import Context
from Meow import StdTrees


import os
import os.path

TDLCompileOption("filename_pattern","Filename pattern",["beam_%(xy)s_%(reim)s.fits"],more=str,doc="""<P>
  Pattern for beam FITS filenames. The real and imaginary parts of each element of the beam Jones matrix should be stored in a separate FITS file. The following tokens in the pattern will be
  substituted:</P>
  <UL>
  <LI><B>%(xy)s</B> or <B>%(corr)s</B>: the correlation index, in lowercase 
  (xx/xy/yx/yy or rr/rl/lr/ll)</LI>
  <LI><B>%(XY)s</B> or <B>%(CORR)s</B>: the correlation index, in uppercase 
  (XX/XY/YX/YY or RR/RL/LR/LL)</LI>
  <LI><B>%(reim)s</B> or <B>%(ReIm)s</B> or <B>%(REIM)s</B>: "re", "Re" or "RE" for real, "im", "Im" or "IM" for imaginary</LI>
  <LI><B>%(realimag)s</B> or <B>%(RealImag)s</B> or <B>%(REALIMAG)s</B>: "real", "Real" or "REAL" for real, "imag", "Imag" or "IMAG" for imaginary</LI>
  <UL>""");
TDLCompileOption("missing_is_null","Use null for missing beam elements",True,doc="""<P>
  If True, then 0 will be used when a beam file is missing. Useful if you only have beams for
  the diagonal elements (XX/YY), or only real beams.
  </P>""");
TDLCompileOption("norm_beams","Normalize beams",True,doc="""<P>
  If True, then the beam gains will be normalized to unity.)
  </P>""");
TDLCompileOption("sky_rotation","Include sky rotation",True,doc="""<P>
  If True, then the beam will rotate on the sky with parallactic angle. Use for e.g. alt-az mounts.)
  </P>""");

CORRS = Context.correlations;
REIM = "re","im";
REALIMAG = dict(re="real",im="imag");

from . import Utils

def make_beam_filename (corr,reim):
  """Makes beam filename for the given correlation and real/imaginary component (one of "re" or "im")"""
  # make tokens dictionary for substitution into the filename pattern
  tokens = dict(
    corr=corr.lower(),xy=corr.lower(),CORR=corr.upper(),XY=corr.upper(),
    reim=reim.lower(),REIM=reim.upper(),ReIm=reim.title(),
    realimag=REALIMAG[reim].lower(),REALIMAG=REALIMAG[reim].upper(),
    RealImag=REALIMAG[reim].title(),
  );
  filename = Utils.substitute_pattern(filename_pattern,**tokens);
  # apply % operation, for comparison with old %(var)s syntax.
  return filename%tokens;

def make_beam_nodes (beam):
  """Makes appropriate nodes to read in beams and put together a Jones matrix.
  beam is the base node.""";
  for corr in CORRS:
    # make FITS images or nulls for real and imaginary part
    for reim in REIM:
      filename = make_beam_filename(corr,reim);
      if os.path.exists(filename):
        beam(corr,reim) << Meq.FITSImage(filename=filename,cutoff=1.0,mode=2);
      elif missing_is_null:
        print("Can't find %s, using null"%filename);
        beam(corr,reim) << 0;
      else:
        raise ValueError("beam file %s does not exist"%filename);
    # put it together into a complex voltage gain
    beam(corr) << Meq.ToComplex(*[beam(corr,ri) for ri in REIM]);
    if norm_beams:
      beam(corr)("pol") << beam(corr,"re") * beam(corr,"re") + beam(corr,"im") * beam(corr,"im")
  if norm_beams:
    beam("sqrt_x") << Meq.Sqrt(beam("XX")("pol") + beam("XY")("pol"))
    beam("norm_x") << Meq.Max(beam("sqrt_x"),reduction_axes=["L", "M"])
    beam("sqrt_y") << Meq.Sqrt(beam("YX")("pol") + beam("YY")("pol"))
    beam("norm_y") << Meq.Max(beam("sqrt_y"),reduction_axes=["L", "M"])
    for corr in CORRS:
      if corr == "XX" or corr == "XY":
        beam(corr)("norm") << beam(corr) / beam("norm_x")
      if corr == "YX" or corr == "YY":
        beam(corr)("norm") << beam(corr) / beam("norm_y")
    beam << Meq.Matrix22(*[beam(corr)("norm") for corr in CORRS]);
  else:
  # now put it together into one beam
    beam << Meq.Matrix22(*[beam(corr) for corr in CORRS]);

def compute_jones (Jones,sources,stations=None,pointing_offsets=None,inspectors=[],label='E',**kw):
  stations = stations or Context.array.stations;
  ns = Jones.Subscope();

  # make the gridded beam pattern node
  beam = ns.beam;
  make_beam_nodes(beam);
  # declare an inspector node for the Jons matrix -- will be defined below
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
        # now make the resampler/compounder combo. dep_mask=0xFF is a little cache-defeating magic
        # (becaise I'm too lazy to fiure out the proper symdep to put in)
        res = Jones("res",src,p) << Meq.Resampler(ns.beam,dep_mask=0xFF);
        Jones(src,p) << Meq.Compounder(lm,res,common_axes=[hiid('l'),hiid('m')]);
    else:
      # compute E for a given source. Use Resampler (though this is not really adequate!)
      res = Jones("res",src) << Meq.Resampler(ns.beam,dep_mask=0xFF);
      Jones(src) << Meq.Compounder(src.direction.lm(),res,common_axes=[hiid('l'),hiid('m')]);
      for p in stations:
        Jones(src,p) << Meq.Identity(Jones(src))

  # add a bookmark for the beam pattern
  Meow.Bookmarks.Page("%s-Jones beam pattern"%label).add(beam,viewer="Result Plotter");

  # define an inspector
  if sky_rotation or pointing_offsets:
    # Jones inspector is per-source, per-station
    insp << StdTrees.define_inspector(Jones,sources,stations,label=label);
  else:
    # Jones inspector is per-source
    insp << StdTrees.define_inspector(Jones,sources,label=label);


  return Jones;
