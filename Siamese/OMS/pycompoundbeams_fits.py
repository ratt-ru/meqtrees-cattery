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
import random
import math
DEG = math.pi/180;

TDLCompileOption("filename_pattern","Filename pattern",["beam_$(elem)_$(xy)_$(reim).fits"],more=str,doc="""<P>
  Pattern for beam FITS filenames. The real and imaginary parts of each element of the beam Jones matrix should be stored in a separate FITS file. A number of variables will be substituted in the beam pattern, these may be introduced as
  "$var" or "$(var)". Use "$$" for a "$" character. The following variables are recognized:</P>
  <UL>
  <LI><B>$elem</B>: element number (0-based).</LI>
  <LI><B>$ele1</B>: element number (1-based).</LI>
  <LI><B>$xy</B>: x or y, for the x or y component of the beam.</LI>
  <LI><B>$XY</B>: X or Y, likewise.</LI>
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
TDLCompileOption("ampl_interpolation","Use amplitude interpolation for beams",False,doc="""<P>
Check the box if you want beam interpolation done with amplitude. The default is to just do real and imaginary voltages separately.</P>""");
TDLCompileOption("sky_rotation","Include sky rotation",True,doc="""<P>
  If True, then the beam will rotate on the sky with parallactic angle. Use for e.g. alt-az mounts.)
  </P>""");
TDLCompileOption("num_elements","Number of beam elements (for FPA/AA compound beams)",[None],more=int,doc="""<P>
  Set to None for a single-pixel feed, so that a single beam pattern is used. Otherwise, set to the number
  of elements in the FPA or AA. Each element's beam will be different.
  """);
TDLCompileOption("weight_filename_x","Initial weights for X beam",TDLFileSelect("*.bw"),more=int,doc="""<P>
  Beam weights will be read from this file -- a whitespace-separated list of N complex numbers (as e.g. "1.5+0.1j") is expected</P>"""); 
TDLCompileOption("weight_filename_y","Initial weights for Y beam",TDLFileSelect("*.bw"),more=int,doc="""<P>
  Beam weights will be read from this file -- a whitespace-separated list of N complex numbers (as e.g. "1.5+0.1j") is expected</P>"""); 
TDLCompileOption("min_ampl_var","Minimum amplitude variation, dB",[0,0.2],more=float);
TDLCompileOption("max_ampl_var","Maximum amplitude variation, dB",[0,0.5],more=float);
TDLCompileOption("min_phase_var","Minimum phase variation, deg",[0,5],more=float);
TDLCompileOption("max_phase_var","Maximum phase variation, deg",[0,10],more=float);
TDLCompileOption("min_period_var","Minimum variation period, hours",2,more=float);
TDLCompileOption("max_period_var","Maximum variation period, hours",24,more=float);
TDLCompileOption("do_correct","Correct for gain of first source",False,doc="""<P>If true, then the gain towards the first source in
      the model will be divided out. This simulates (to first order) correction by selfcal.</P>""");

CORRS = Context.correlations;
REIM = "re","im";
REALIMAG = dict(re="real",im="imag");

from Siamese.OMS import Utils

def make_beam_filename (filename_pattern,elem,xy,reim):
  """Makes beam filename for the given x/y, element (0-based), and real/imaginary component (one of "re" or "im")"""
  return Utils.substitute_pattern(filename_pattern,
    elem="%d"%elem,ele1="%d"%(elem+1),
    xy=xy.lower(),XY=xy.upper(),
    reim=reim.lower(),REIM=reim.upper(),ReIm=reim.title(),
    realimag=REALIMAG[reim].lower(),REALIMAG=REALIMAG[reim].upper(),
    RealImag=REALIMAG[reim].title());

def make_beam_node (beam,pattern,*children):
  """Makes beam interpolator node for the given filename pattern.""";
  filename_real = [];
  filename_imag = [];
  for corr in "x","y":
    for elem in range(num_elements):
      # make FITS images or nulls for real and imaginary part
      filename_real.append(make_beam_filename(pattern,elem,corr,'re'));
      filename_imag.append(make_beam_filename(pattern,elem,corr,'im'));
  # if we end up with only one unique filename, make a scalar-mode interpolator
  # now make interpolator node
  import CompoundInterpolatedBeams
  beam << Meq.PyNode(class_name="FITSCompoundBeamInterpolatorNode",module_name=CompoundInterpolatedBeams.__file__,
                     filename_real=filename_real,filename_imag=filename_imag,normalize=normalize_gains,
                     missing_is_null=missing_is_null,spline_order=spline_order,verbose=0, 
                     ampl_interpolation=ampl_interpolation,
                     children=children);

def compute_jones (Jones,sources,stations=None,pointing_offsets=None,inspectors=[],label='E',**kw):
  stations = stations or Context.array.stations;
  ns = Jones.Subscope();
  JE = Jones("elem");

  # loop over sources to create per-element beamgains
  for src in sources:
    # If sky rotation and/or pointing offsets are in effect, we have a per-station beam.
    # Otherwise the beam is the same for all stations.
    if sky_rotation or pointing_offsets:
      for p in stations:
        lm = src.direction.lm();
        # apply rotation to put sources into the antenna frame
        if sky_rotation:
          xyz = Context.array.xyz(p);
          pa_rot = Context.observation.phase_centre.pa_rot(xyz);
          lm = ns.lmrot(src,p) <<  Meq.MatrixMultiply(pa_rot,src.direction.lm());
        # apply offset (so pointing offsets are interpreted in the azel frame, if rotating)
        if pointing_offsets:
          lm = ns.lmoff(src,p) << lm + pointing_offsets(p);
        # now make the beam node
        make_beam_node(JE(src,p),filename_pattern,lm);
    else:
      make_beam_node(JE(src),filename_pattern,src.direction.lm());
  
  # now create per-station weights
  wx = map(complex,open(weight_filename_x).read().split());
  wy = map(complex,open(weight_filename_y).read().split());
  
  if len(wx) != num_elements or len(wy) != num_elements:
    raise ValueError,"weights files contain %d (X) and %d (Y) complex weights, %d expected"%(len(wx),len(wy),num_elements);
  
  # w0:x and w0:y are the nominal weight vectors
  ns.w0('x') << Meq.Constant(value=wx);
  ns.w0('y') << Meq.Constant(value=wy);
 
  # create perturbed weights
  a0 = 10**(min_ampl_var/20)-1;
  a1 = 10**(max_ampl_var/20)-1;
  for p in stations:
    for xy in 'x','y':
      werr = ns.werr(xy,p);
      # amplitude and phase period and offset
      for ap in 'ampl','phase':
        p0 = werr("period",ap) << Meq.Constant(value=[random.uniform(min_period_var*3600,max_period_var*3600)/(2*math.pi) for i in range(num_elements)]);
        t0 = werr("offset",ap) << Meq.Constant(value=[random.uniform(0,2*math.pi) for i in range(num_elements)]);
        werr("sin",ap) << Meq.Sin((Meq.Time()/p0)+t0);
      # amplitude excursion
      e0 = werr("maxampl") << Meq.Constant(value=[random.uniform(a0,a1) for i in range(num_elements)]);
      ep = werr("maxphase") << Meq.Constant(value=[random.uniform(min_phase_var*DEG,max_phase_var*DEG) for i in range(num_elements)]);
      # weight errors
      werr << Meq.Polar(1+e0*werr("sin","ampl"),ep*werr("sin","phase"));
      ns.weight(xy,p) << ns.w0(xy)*werr;
 
  # put these together into Jones matrices
  for src in sources:
    jesrc = JE(src);
    for p in stations:
      if sky_rotation or pointing_offsets:
        jesrc = JE(src,p);
      # jesrc returns a (2,N) matrix of element gains towards this source. Multiply this by the weights to get the "x" and "y" beams
      ex = Jones(src,p,"x") << Meq.MatrixMultiply(jesrc,ns.weight("x",p));
      ey = Jones(src,p,"y") << Meq.MatrixMultiply(jesrc,ns.weight("y",p));
      if do_correct:
        J0 = Jones(src,p,"ref") << Meq.Composer(ex,ey,dims=[2,2]);
        if src is sources[0]:
          Jones(src,p) << 1;
          Jones(src,p,"inv") << Meq.MatrixInvert22(J0);
        else:
          Jones(src,p) << Meq.MatrixMultiply(Jones(sources[0],p,"inv"),J0);
      else:
        Jones(src,p) << Meq.Composer(ex,ey,dims=[2,2]);
  
  # declare an inspector node for the Jones matrix -- will be defined below
  insp = Jones.scope.inspector(label);

  # define inspectors
  insp << StdTrees.define_inspector(Jones,sources,stations,label=label);
  insp1 = insp("werr") << StdTrees.define_inspector(ns.werr,("x","y"),stations,label="%s weight drifts"%label);
  insp2 = insp("weights") << StdTrees.define_inspector(ns.weight,("x","y"),stations,label="%s weights"%label);

  inspectors += [ insp,insp1,insp2 ];

  return Jones;
