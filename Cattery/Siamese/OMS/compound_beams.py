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

"""<P>This module implements voltage beam patterns that are read from several 
FITS files and interpolated via a PyNode. Multiple beam patterns are read
in and added together.</P>

<P align="right">Author: O. Smirnov &lt;<tt>smirnov@astron.nl</tt>&gt;</P>""";
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division
from Timba.TDL import *
import Meow
from Meow import Context
from Meow import StdTrees


import os
import os.path
import random
import math
import pickle
DEG = math.pi/180;

MAXPATT = 5;
TDLCompileOption("num_patterns","Number of beam patterns to combine",list(range(1,MAXPATT+1)));
for n in range(1,MAXPATT+1):
  TDLCompileOption("filename_pattern_%d"%n,"Filename pattern %d"%n,["beam_$(xy)_$(reim).fits"],more=str,doc="""<P>
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
TDLCompileOption("ampl_interpolation","Use amplitude interpolation for beams",False,doc="""<P>
Check the box if you want beam interpolation done with amplitude. The default is to just do real and imaginary voltages separately.</P>""");
TDLCompileOption("l_beam_offset","Offset beam pattern in L (deg)",[0.0], more=float,
doc="""<P>By default,the beam reference position (as given by the FITS header) is placed at l=m=0 on the sky, i.e. 
at the phase centre. You can use this option to offset the beam pattern.</P>"""),
TDLCompileOption("m_beam_offset","Offset beam pattern in M (deg)",[0.0], more=float,
doc="""<P>By default,the beam reference position (as given by the FITS header) is placed at l=m=0 on the sky, i.e. 
at the phase centre. You can use this option to offset the beam pattern.</P>"""),
TDLCompileOption("sky_rotation","Include sky rotation",True,doc="""<P>
  If True, then the beam will rotate on the sky with parallactic angle. Use for e.g. alt-az mounts.)
  </P>""");
  

CORRS = Context.correlations;
REIM = "re","im";
REALIMAG = dict(re="real",im="imag");

from Siamese.OMS import Utils
from Siamese.OMS.pybeams_fits import make_beam_filename

def make_beam_node (beam,pattern,l_offset,m_offset,*children):
  """Makes beam interpolator node for the given filename pattern.""";
  filename_real = [];
  filename_imag = [];
  for corr in "x","y":
    for elem in range(1,num_patterns):
      # make FITS images or nulls for real and imaginary part
      filename_real.append(make_beam_filename(pattern,elem,corr,'re'));
      filename_imag.append(make_beam_filename(pattern,elem,corr,'im'));
  # if we end up with only one unique filename, make a scalar-mode interpolator
  # now make interpolator node
  from . import CompoundInterpolatedBeams
  beam << Meq.PyNode(class_name="FITSCompoundBeamInterpolatorNode",module_name=CompoundInterpolatedBeams.__file__,
                     filename_real=filename_real,filename_imag=filename_imag,
                     l_0=l_offset,m_0=m_offset,
                     normalize=normalize_gains,
                     missing_is_null=missing_is_null,spline_order=spline_order,verbose=0,
                     ampl_interpolation=ampl_interpolation,
                     children=children);

def make_norm (J,Jnorm):
  """Returns the "norm" of a Jones matrix, as tr(|AA^H|)/2.""";
  J('sq') << Meq.MatrixMultiply(J,J("conj")<<Meq.ConjTranspose(J));
  ja = J('abssq') << Meq.Abs(J('sq'));
  jxx = J('abssq','xx') << Meq.Selector(ja,index=0);
  jyy = J('abssq','yy') << Meq.Selector(ja,index=3);
  Jnorm << Meq.Sqrt((jxx+jyy)/2);

def compute_jones (Jones,sources,stations=None,pointing_offsets=None,inspectors=[],label='E',**kw):
  stations = stations or Context.array.stations;
  ns = Jones.Subscope();
  JE = Jones("elem");

  per_station = sky_rotation or pointing_offsets;

  # read offsets file
  if read_offsets_file:
    offsets = list(map(float,open(offsets_file).read().split()));
    if len(offsets) < (beam_number+1)*2:
      raise ValueError("beam number %d not found in offsets file"%beam_number);
    l_offset,m_offset = offsets[beam_number*2:(beam_number+1)*2];
    if invert_l:
      l_offset = -l_offset;
  else:
    l_offset,m_offset = 0,0;

  # loop over sources to create per-element beamgains
  for src in sources:
    # If sky rotation and/or pointing offsets are in effect, we have a per-station beam.
    # Otherwise the beam is the same for all stations.
    if per_station:
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
        make_beam_node(JE(src,p),filename_pattern,l_offset,m_offset,lm);
    else:
      make_beam_node(JE(src),filename_pattern,l_offset,m_offset,src.direction.lm());


  # now load weights
  wx = pickle.load(open(weight_filename_x, "rb"));
  wy = pickle.load(open(weight_filename_y, "rb"));
  if wx.shape[1] != num_elements or wy.shape[1] != num_elements:
    raise ValueError("""weights files contain weights for %d (X) and %d (Y) complex
                      elements, %d expected"""%(wx.shape[1],wy.shape[1],num_elements));
  if beam_number > wx.shape[0] or beam_number > wy.shape[0]:
    raise ValueError("beam number %d not found in weights files"%beam_number);

  # w0:x and w0:y are the nominal weight vectors
  ns.w0('x') << Meq.Constant(value=wx[beam_number,:]);
  ns.w0('y') << Meq.Constant(value=wy[beam_number,:]);

  if sim_element_errors:
    # create perturbed weights
    a0 = 10**(min_ampl_var/20)-1;
    a1 = 10**(max_ampl_var/20)-1;
    for p in stations:
      for xy in 'x','y':
        werr = ns.werr(xy,p);
        # amplitude and phase period and offset
        for ap in 'ampl','phase':
          p0 = werr("period",ap) << Meq.Constant(value=[random.uniform(min_period_var*3600,max_period_var*3600)/(2*math.pi) for i in range(num_elements)]);
          if start_phased_up:
            werr("sin",ap) << Meq.Sin(Meq.Time()/p0);
          else:
            t0 = werr("offset",ap) << Meq.Constant(value=[random.uniform(0,2*math.pi) for i in range(num_elements)]);
            werr("sin",ap) << Meq.Sin((Meq.Time()/p0)+t0);
          
        # amplitude excursion
        e0 = werr("maxampl") << Meq.Constant(value=[random.uniform(a0,a1) for i in range(num_elements)]);
        ep = werr("maxphase") << Meq.Constant(value=[random.uniform(min_phase_var*DEG,max_phase_var*DEG) for i in range(num_elements)]);
        # weight errors
        werr << Meq.Polar(1+e0*werr("sin","ampl"),ep*werr("sin","phase"));
        ns.weight(xy,p) << ns.w0(xy)*werr;

  # compute matrix norms based on nominal matrices
  if per_station:
    quallist = [ [src,p] for src in sources for p in stations ];
  else:
    quallist = [ [src] for src in sources ];
  for qq in quallist:
    ex0 = Jones(*(qq+["x0"])) << Meq.MatrixMultiply(JE(*qq),ns.w0("x"));
    ey0 = Jones(*(qq+["y0"])) << Meq.MatrixMultiply(JE(*qq),ns.w0("y"));
    J0 = Jones(*(qq+["nominal"])) << Meq.Composer(ex0,ey0,dims=[2,2]);
    if do_normalize or (qq[0] is sources[0] and not do_correct):
      make_norm(J0,Jones(*(qq+["norm"])));

  if do_normalize:
    # norm towards source is average per-station norm
    if per_station:
      for src in sources:
        Jones(src,"norm") << Meq.Add(*[Jones(src,p,"norm") for p in stations])/len(stations);
  elif not do_correct:
    if per_station:
      for src in sources:
        Jones(sources[0],"norm") << Meq.Add(*[Jones(sources[0],p,"norm") for p in stations])/len(stations);
    for src in sources[1:]:
      Jones(src,"norm") << Meq.Identity(Jones(sources[0],"norm"));
  else:
    for src in sources:
      Jones(src,"norm") << 1;

  # put these together into Jones matrices
  for src in sources:
    # jesrc/JJ will eventually point to the unqualified element beam node
    # of the unqualified jones node. Depending on whether we're per-station or
    # not, this is qualified with src,p or just src.
    jesrc = JE(src);
    jnorm = Jones(src,"norm");
    JJ = Jones(src);

    if sim_element_errors:
      for p in stations:
        if per_station:
          jesrc = JE(src,p);
          JJ = JJ(p);
        # jesrc returns a (2,N) matrix of element gains towards this source.
        # Multiply this by the weights to get the "x" and "y" beams
        ex = Jones(src,p,"x") << Meq.MatrixMultiply(jesrc,ns.weight("x",p));
        ey = Jones(src,p,"y") << Meq.MatrixMultiply(jesrc,ns.weight("y",p));
        if do_correct:
          J0 = Jones(src,p,"ref") << Meq.Composer(ex,ey,dims=[2,2])/jnorm;
          if src is sources[0]:
            Jones(src,p) << Meq.Constant(value=[1,0,0,1],dims=[2,2]);
            Jones(src,p,"inv") << Meq.MatrixInvert22(J0);
          else:
            Jones(src,p) << Meq.MatrixMultiply(Jones(sources[0],p,"inv"),J0);
        else:
          Jones(src,p) << Meq.Composer(ex,ey,dims=[2,2])/jnorm;
    # no element errors, use nominal beam
    else:
      if per_station:
        for p in stations:
          if do_correct:
            J0 = Jones(src,p,"ref") << Jones(src,p,"nominal")/jnorm;
            if src is sources[0]:
              Jones(src,p) << Meq.Constant(value=[1,0,0,1],dims=[2,2]);
              Jones(src,p,"inv") << Meq.MatrixInvert22(J0);
            else:
              Jones(src,p) << Meq.MatrixMultiply(Jones(sources[0],p,"inv"),J0);
          else:
            Jones(src,p) << Jones(src,p,"nominal")/jnorm;
      else:
        if do_correct:
          J0 = Jones(src,"ref") << Jones(src,"nominal")/jnorm;
          if src is sources[0]:
            Jones(src) << Meq.Constant(value=[1,0,0,1],dims=[2,2]);
            Jones(src,"inv") << Meq.MatrixInvert22(J0);
          else:
            Jones(src) << Meq.MatrixMultiply(Jones(sources[0],"inv"),J0);
        else:
          Jones(src) << Jones(src,"nominal")/jnorm;
        for p in stations:
          Jones(src,p) << Meq.Identity(Jones(src));


  # declare an inspector node for the Jones matrix -- will be defined below
  insp = Jones.scope.inspector(label);

  # define inspectors
  insp << StdTrees.define_inspector(Jones,sources,stations,label=label);
  inspectors += [ insp ];
  if sim_element_errors:
    insp1 = insp("werr") << StdTrees.define_inspector(ns.werr,("x","y"),stations,label="%s weight drifts"%label);
    insp2 = insp("weights") << StdTrees.define_inspector(ns.weight,("x","y"),stations,label="%s weights"%label);

    inspectors += [ insp1,insp2 ];

  return Jones;
