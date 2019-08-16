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
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division
from Timba.TDL import *
from Timba import pynode
from Timba.Meq import meq
from Timba import mequtils

import Meow
from Meow import Context
from Meow import StdTrees

import os
import os.path
import math
import numpy
import random
import Kittens.utils
_verbosity = Kittens.utils.verbosity(name="vb");
dprint = _verbosity.dprint;
dprintf = _verbosity.dprintf;


DEG = math.pi/180;

from . import EMSSVoltageBeam
from . import InterpolatedVoltageBeam
from .InterpolatedVoltageBeam import unite_shapes,unite_multiple_shapes

SYM_SEPARATE = None;
SYM_X = "X";
SYM_Y = "Y";
COORD_LM = "lm";
COORD_THETAPHI = "thetaphi";


TDLCompileOption("filename_pattern","Filename pattern",["beam_$(hv).pat"],more=str,doc="""<P>
  Pattern for beam filenames. Each beam file contains two elements of the E-Jones matrix. A number of variables will be
  substituted in the filename pattern, these may be introduced as
  "$var" or "$(var)". Use "$$" for a "$" character. The following variables are recognized:</P>
  <UL>
  <LI><B>$xy</B> or <B>$XY</B>: x or y, or X or Y, for the first or second row of the E-Jones matrix</LI>
  <LI><B>$hv</B> or <B>$HV</B>: v or h, or V or H, for the first or second row of the E-Jones matrix</LI>
  <LI><B>$label</B>: pattern label. Multiple-resolution patterns may be specified, they will all be added together.</LI>
  <LI><B>$freq</B>: frequency label. Multiple per-frequency patterns may be specified, they will be interpolated in.</LI>
  <UL>""");
TDLCompileOption("pattern_labels","Beam pattern labels",[None],more=str);
TDLCompileOption("freq_labels","Beam frequency labels",[None],more=str);
TDLCompileOption("beam_symmetry","Use 90-degree symmetry",{
  None:"no, load separate X/Y patterns",
  SYM_X:"yes, XX XY (v) given",
  SYM_Y:"yes, YX YY (h) given"},doc=
  """<P>If the beam has 90-degree rotational symmetry, then its XX/XY and YX/YY components can be derived
  from one another via a 90 degree rotation. If this is the case, then you may use just one set of pattern files
  (either one, as specified by this setting), and let the interpolator apply rotation to derive the other set.
  The beam pattern filename above should't contain a substitutable $(hv) or $(xy) element then.
  </P>
  """);
TDLCompileOption("normalization_factor","Normalize beam patterns by",[1],more=float,
  doc="""<P>Beam patterns will be divided by the specified normalization factor.</P>""");
TDLCompileOption("rotate_xy","Rotate polarization vectors from theta-phi to xy frame",True,
  doc="""<P>Set this if the beam amplitudes are specified in the <i>&theta;&phi;</i> frame, and need to be rotated back to the
  <I>xy</I> frame.</P>""");
TDLCompileOption("spline_order","Spline order for interpolation",[1,2,3,4,5],default=3);
TDLCompileOption("hier_interpol","Use hierarchical interpolation (lm, then frequency)",True);

coord_opt = TDLCompileOption("interpol_coord","Coordinate interpolation",
  {COORD_LM:"l,m (fast)",COORD_THETAPHI:"theta,phi (all-sky accurate)"},
  doc="""<P>This determines the coordinate conversion strategy for beam interpolation. Using l,m coordinates is faster,
  but less accurate in the far sidelobes (and cannot handle backlobes, i.e. sources over 90&deg; at all). Use theta, phi 
  coordinates for a more accurate but slower interpolation in the far sidelobes, or if you have sources more that 90&deg; 
  away from the pointing centre.</P>""");

l0opt = TDLCompileOption("l_beam_offset","Offset beam pattern in L (deg)",[0.0], more=float,
doc="""<P>By default,the beam reference position (<i>&theta;=&phi;=0</i>) is placed at <I>l=m=</I>0 on the sky, i.e.
at the phase centre. You can use this option to offset the beam pattern.</P>""");

m0opt = TDLCompileOption("m_beam_offset","Offset beam pattern in M (deg)",[0.0], more=float,
doc="""<P>By default,the beam reference position (<i>&theta;=&phi;=0</i>) is placed at <I>l=m=</I>0 on the sky, i.e.
at the phase centre. You can use this option to offset the beam pattern.</P>""");

coord_opt.when_changed(lambda val:(l0opt.show(val is COORD_LM),m0opt.show(val is COORD_LM)));

TDLCompileOption("sky_rotation","Include sky rotation",True,doc="""<P>
  If True, then the beam will rotate on the sky with parallactic angle. Use for alt-az mounts.</P>""");
TDLCompileOption("horizon_masking","Include horizon masking",False,doc="""<P>
  If True, then horizon masking will be included. This may slow things down.</P>""");
TDLCompileOption("randomize_rotation","Include random rotational offsets",False,doc="""<P>
  If True, then each station's beam will be rotated by an arbitrary amount. This can help randomize the sidelobes.</P>""");
TDLCompileOption("verbose_level","Debugging message level",[None,1,2,3],more=int);

def _cells_grid (obj,axis):
  """helper function to get a grid out of the cells object. Returns None if none is found"""
  if hasattr(obj,'cells') and hasattr(obj.cells.grid,axis):
    return obj.cells.grid[axis];
  else:
    return None;


class EMSSPolarBeamInterpolatorNode (pynode.PyNode):
  def __init__ (self,*args):
    pynode.PyNode.__init__(self,*args);
    # Maintain a global dict of VoltageBeam objects per each filename set, so that we reuse them
    # We also reset this dict each time a node is created (otherwise the objects end up being reused
    # even after the tree has been rebuilt.)
    global _voltage_beams;
    _voltage_beams = {};
    self._vbs = None;

  def update_state (self,mystate):
    """Standard function to update our state""";
    mystate('filename',[]);
    mystate('spline_order',3);
    mystate('hier_interpol',True);
    mystate('interpol_lm',True);
    mystate('l_0',0.0);
    mystate('m_0',0.0);
    mystate('verbose',0);
    mystate('missing_is_null',False);
    mystate('beam_symmetry',None);
    mystate('normalization_factor',1);
    mystate('rotate_xy',True);
    # other init
    mequtils.add_axis('l');
    mequtils.add_axis('m');
    self._freqaxis = mequtils.get_axis_number("freq");
    _verbosity.set_verbose(self.verbose);

  def init_voltage_beams (self):
    """initializes VoltageBeams supplied to this node"""
    # maintain a global dict of VoltageBeam objects per each filename set, so that we reuse them
    if self._vbs is not None:
      return self._vbs;
    # we get VoltageBeam objects from global dict, or init new ones if not already defined
    global _voltage_beams;
    self._vbs = [];
    # setup rotation and y-arg combinations according to symmetry modes
    # only V/X beam supplied. H/Y must be rotated by 90 degrees and columns swapped.
    if self.beam_symmetry == SYM_X:
      rotations = 0,90;
      yargflips = False,True;
    # only H/Y beam supplied. V/X rotated/swapped
    elif self.beam_symmetry == SYM_Y:
      rotations = 90,0;
      yargflips = True,False;
    # else both beams supplied
    else:
      rotations = 0,0;
      yargflips = False,False;
    # loop over files
    for per_label in self.filename:
      vbmat = [];
      for xy,rotate,yargflip,files_per_xy in zip(('x','y'),rotations,yargflips,per_label):
        for xy1,yarg0 in zip(('x','y'),(False,True)):
          yarg = yarg0^yargflip;
          vbkey = tuple(files_per_xy),yarg,rotate;
          vb = _voltage_beams.get(vbkey)
          if vb is None:
            dprint(1,"Loading beams for %s%s from"%(xy,xy1),files_per_xy);
            vb = _voltage_beams[vbkey] = EMSSVoltageBeam.EMSSVoltageBeamPS(files_per_xy,y=yarg,rotate=rotate,
                        spline_order=self.spline_order,hier_interpol=self.hier_interpol,
                        rotate_xy=self.rotate_xy,normalization_factor=self.normalization_factor,
                        proj_theta=False,verbose=self.verbose);
          vbmat.append(vb);
      self._vbs.append(vbmat);
    return self._vbs;

  def interpolate_per_source (self,lm_list,dl,dm,grid,vbs,thetaphi=False,rotate=None,masklist=None):
    # Loops over all source coordinates in the lm tensor, interpolates beams for them, and returns the resulting vellsets.
    # lm is a list of [[l0,m0],[l1,m1],...] source coordinates
    # This is the "fail-safe" method, as it interpolates per-source, and therefore allows for the source lm arrays to have different
    # time/frequency dependencies per source.
    vellsets = [];
    for isrc,(l,m) in enumerate(lm_list):
      # apply pointing offsets, if any
      if dl is not None:
        # unite shapes just in case, since l/m and dl/dm may have time/freq axes
        l,dl,m,dm = unite_multiple_shapes(l,dl,m,dm);
      # if mask is set, make sure it has the same shape
      mask = masklist[isrc] if masklist is not None else None;
      if mask is not None:
        l,m,mask = unite_multiple_shapes(l,m,mask);
      grid['l'],grid['m'] = l,m;
      # loop over all 2x2 matrices (we may have several, they all need to be added)
      E = [None]*4;
      for vbmat in vbs:
        for i,vb in enumerate(vbmat):
          beam = vb.interpolate(freqaxis=self._freqaxis,lm=self.interpol_lm,thetaphi=thetaphi,rotate=rotate,mask=mask,**grid);
          if E[i] is None:
            E[i] = meq.complex_vells(beam.shape);
          E[i][...] += beam[...];
      # make vellsets
      for ej in E:
        vellsets.append(meq.vellset(ej));
    return vellsets;
    
  def transform_coordinates (self,l,m):
    return l,m;

  def interpolate_batch (self,lm_list,dl,dm,grid,vbs,thetaphi=False,rotate=None,masklist=None):
    # A faster version of interpolate_per_source(), which assumes that all lm's (as well as the masks, if given) 
    # have the same shape, and stacks them into a single array for a single interpolation call.
    # If there's a shape mismatch, it'll fall back to interpolate_per_source
    # 'maskarr', if given, should be an list of per-source mask arrays 
    nsrc = len(lm_list);
    maskcube = None;
    for isrc,(l,m) in enumerate(lm_list):
      mask = masklist[isrc] if masklist is not None else None;
      # apply pointing offsets, if any
      if dl is not None:
        # unite shapes just in case, since l/m and dl/dm may have time/freq axes
        l,dl,m,dm = unite_multiple_shapes(l,dl,m,dm);
        l,m = l-dl,m-dm;
      # unite l,m shapes just in case, and transform
      l,m,mask = unite_multiple_shapes(l,m,mask);
      if not isrc:
        lm_shape = l.shape;
        cubeshape = [nsrc]+list(lm_shape);
        lcube = numpy.zeros(cubeshape,float);
        mcube = numpy.zeros(cubeshape,float);
        if masklist is not None:
          maskcube = numpy.zeros(cubeshape,bool);
      else:
        if l.shape != lm_shape:
          dprint(1,"l/m shapes unequal at source %d, falling back to per-source interpolation"%isrc);
          return self.interpolate_per_source(lm_list,dl,dm,grid,vbs,thetaphi=thetaphi,rotate=rotate,masklist=masklist);
      lcube[isrc,...] = l;
      mcube[isrc,...] = m;
      if mask is not None:
        maskcube[isrc,...] = mask;
#        if mask.any():
#          dprint(2,"source %d has %d slots masked"%(isrc,mask.sum()));
    # if 'rotate' is specified, it needs to be promoted to the same cube shape, and ravelled
    if rotate is not None:
      lcube,mcube,maskcube,rotate = unite_multiple_shapes(lcube,mcube,maskcube,rotate.reshape([1]+list(rotate.shape)));
      cubeshape = list(lcube.shape);
      rotate = rotate.ravel();
    # ok, we've stacked things into lm cubes, interpolate
    grid['l'],grid['m'] = lcube.ravel(),mcube.ravel();
    # loop over all 2x2 matrices (we may have several, they all need to be added)
    E = [None]*4;
    for vbmat in vbs:
      for i,vb in enumerate(vbmat):
        beam = vb.interpolate(freqaxis=self._freqaxis,extra_axes=1,thetaphi=thetaphi,rotate=rotate,**grid);
        if E[i] is None:
          E[i] = beam;
        else:
          E[i] += beam;
    # The l/m cubes have a shape of [nsrcs,lm_shape].
    # These are raveled for interpolation, so the resulting Es have a shape of [nsrcs*num_lm_points,num_freq]
    # Reshape them properly. Note that there's an extra "source" axis at the front, so the frequency axis
    # is off by 1.
    if len(cubeshape) <= self._freqaxis+1:
      cubeshape = list(cubeshape) + [1]*(self._freqaxis - len(cubeshape) + 2);
    cubeshape[self._freqaxis+1] = len(grid['freq']);
    E = [ ej.reshape(cubeshape) for ej in E ];
    # apply mask cube, if we had one
    if maskcube is not None:
      # the E planes now have the same shape as the maskcube, but with an extra frequency axis. Promote the maskcube
      # accordingly
      me = unite_multiple_shapes(maskcube,*E);
      maskcube = me[0];
      E = me[1:];
      dprint(2,"maskcube sum",maskcube.sum());
      for ej in E:
        ej[maskcube] = 0;
    # now tease the E's apart plane by plane
    # make vellsets
    vellsets = [];
    for isrc in range(nsrc):
      for ej in E:
        dprint(2,"source %d has %d null gains"%(isrc,(ej[isrc,...]==0).sum()));
        ejplane = ej[isrc,...];
        value = meq.complex_vells(ejplane.shape,ejplane);
        vellsets.append(meq.vellset(value));
    return vellsets;

  def get_result (self,request,lm,dlm=None,azel=None,pa=None):
    # get list of VoltageBeams
    vbs = self.init_voltage_beams();
    # now, figure out the lm and time/freq grid
    # lm may be a 2/3-vector or an Nx2/3 tensor
    dims = getattr(lm,'dims',[len(lm.vellsets)]);
    if len(dims) == 2 and dims[1] in (2,3):
      lm_list = [ (lm.vellsets[i].value,lm.vellsets[i+1].value) for i in range(0,len(lm.vellsets),dims[1]) ];
      tensor = True;
    elif len(dims) == 1 and dims[0] in (2,3):
      lm_list = [ (lm.vellsets[0].value,lm.vellsets[1].value) ];
      tensor = False;
    else:
      raise TypeError("expecting a 2/3-vector or an Nx2/3 matrix for child 0 (lm)");
    nsrc = len(lm_list);
    dl = dm = masklist = rotate = None;
    # pointing offsets are optional
    if dlm is not None and len(dlm.vellsets) == 2:
      dl,dm = dlm.vellsets[0].value,dlm.vellsets[1].value;
    # az/el is optional. If specified, then horizon masking is enabled
    if azel is not None and len(azel.vellsets) > 1:
      if len(azel.vellsets) != nsrc*2:
        raise TypeError("expecting a Nx2 matrix for child 2 (azel)");
      masklist = [ azel.vellsets[i].value<0 for i in range(1,nsrc*2,2) ];
    # PA is optional. If specified, then this is the rotation angle for interpolation
    if pa is not None:
      if len(pa.vellsets) != 1:
        raise TypeError("expecting a single value for child 3 (pa)");
      rotate = pa.vellsets[0].value;
    # setup grid dict that will be passed to VoltageBeam.interpolate
    grid = dict();
    for axis in 'time','freq':
      values = _cells_grid(lm,axis);
      if values is None:
        values = _cells_grid(request,axis);
      if values is not None:
        if numpy.isscalar(values) or values.ndim == 0:
          values = [float(values)];
        grid[axis] = values;
    # accumulate per-source EJones tensor
    vellsets = self.interpolate_batch(lm_list,dl,dm,grid,vbs,rotate=rotate,masklist=masklist);
    # create result object
    cells = request.cells if vbs[0][0].hasFrequencyAxis() else getattr(lm,'cells',request.cells);
    result = meq.result(vellsets[0],cells=cells);
    result.vellsets[1:] = vellsets[1:];
    result.dims = (nsrc,2,2) if tensor else (2,2);
    return result;


class EMSSPolarBeamRaDecInterpolatorNode (EMSSPolarBeamInterpolatorNode):

  @staticmethod
  def lm_to_radec (l,m,ra0,dec0):
    """Returns ra,dec corresponding to l,m w.r.t. direction ra0,dec0""";
    # see formula at http://en.wikipedia.org/wiki/Orthographic_projection_(cartography)
    rho = numpy.sqrt(l**2+m**2);
    null = (rho==0);
    if sum(null) == rho.size:
      return ra0,dec0;
    cc = numpy.arcsin(rho);
    ra = ra0 + numpy.arctan2( l*numpy.sin(cc),rho*numpy.cos(dec0)*numpy.cos(cc)-m*numpy.sin(dec0)*numpy.sin(cc) );
    dec = numpy.arcsin( numpy.cos(cc)*numpy.sin(dec0) + m*numpy.sin(cc)*numpy.cos(dec0)/rho );
    ra[null]  = ra0[null];
    dec[null] = dec0[null];
    return ra,dec;
  
  def get_result (self,request,radec,radec0,dlm=None,azel=None,pa=None):
    # get list of VoltageBeams
    vbs = self.init_voltage_beams();
    # now, figure out the radec and time/freq grid
    # radec may be a 2/3-vector or an Nx2/3 tensor
    dims = getattr(radec,'dims',[len(radec.vellsets)]);
    if len(dims) == 2 and dims[1] == 2:
      radec_list = [ (radec.vellsets[i].value,radec.vellsets[i+1].value) for i in range(0,len(radec.vellsets),dims[1]) ];
      tensor = True;
    elif len(dims) == 1 and dims[0] == 2:
      radec_list = [ (radec.vellsets[0].value,radec.vellsets[1].value) ];
      tensor = False;
    else:
      raise TypeError("expecting a 2-vector or an Nx2 matrix for child 0 (radec)");
    nsrc = len(radec_list);
    # radec0 is a 2-vector
    if len(radec0.vellsets) != 2:
      raise TypeError("expecting a 2-vector for child 1 (radec0)");
    ra0,dec0 = radec0.vellsets[0].value,radec0.vellsets[1].value;
    # get offsets and rotations
    masklist = rotate = None;
    # pointing offsets are optional -- apply them in here
    if dlm is not None and len(dlm.vellsets) == 2:
      dl,dm = dlm.vellsets[0].value,dlm.vellsets[1].value;
      ra0,dl,dec0,dm = unite_multiple_shapes(ra0,dl,dec0,dm);
      ra0,dec0 = self.lm_to_radec(dl,dm,ra0,dec0); 
    # az/el is optional. If specified, then horizon masking is enabled
    if azel is not None and len(azel.vellsets) > 1:
      if len(azel.vellsets) != nsrc*2:
        raise TypeError("expecting a Nx2 matrix for child 3 (azel)");
      masklist = [ azel.vellsets[i].value<0 for i in range(1,nsrc*2,2) ];
    # PA is optional. If specified, then this is the rotation angle for interpolation
    if pa is not None:
      if len(pa.vellsets) != 1:
        raise TypeError("expecting a single value for child 4 (pa)");
      rotate = pa.vellsets[0].value;
    thetaphi_list = [];
    # precompute these -- may be functions of time if time-variable pointing is in effect
    sind0,cosd0 = numpy.sin(dec0),numpy.cos(dec0);
    # now go over all coordinates
    for ra,dec in radec_list:
      # This is based on the formula in Tigger/Coordinates.py, except we flip the sign of phi (by swapping ra0 and ra around),
      # since rotation is meant to go N through W here
      dra = numpy.subtract(*unite_shapes(ra,ra0));
      cosra,sinra = numpy.cos(dra),numpy.sin(dra);
      sind0x,sind,cosd0x,cosd,cosra,sinra = unite_multiple_shapes(sind0,numpy.sin(dec),cosd0,numpy.cos(dec),cosra,sinra);
      theta = numpy.arccos(sind0x*sind + cosd0x*cosd*cosra);
      phi = numpy.arctan2(-cosd*sinra,-cosd*sind0x*cosra+sind*cosd0x);
      thetaphi_list.append((theta,phi));
    dprint(2,"theta/phi",thetaphi_list);
    # setup grid dict that will be passed to VoltageBeam.interpolate
    grid = dict();
    for axis in 'time','freq':
      values = _cells_grid(radec,axis);
      if values is None:
        values = _cells_grid(request,axis);
      if values is not None:
        if numpy.isscalar(values) or values.ndim == 0:
          values = [float(values)];
        grid[axis] = values;
    # accumulate per-source EJones tensor
    vellsets = self.interpolate_batch(thetaphi_list,None,None,grid,vbs,thetaphi=True,rotate=rotate,masklist=masklist);
    # create result object
    cells = request.cells;
    result = meq.result(vellsets[0],cells=cells);
    result.vellsets[1:] = vellsets[1:];
    result.dims = (len(radec_list),2,2) if tensor else (2,2);
    return result;


from Siamese.OMS import Utils
XY2HV = dict(x='v',y='h');

def make_beam_filename (filename_pattern,xy,label,freq):
  """Makes beam filename for the given correlation and real/imaginary component (one of "re" or "im")"""
  hv = XY2HV[xy.lower()];
  return Utils.substitute_pattern(filename_pattern,
    xy=xy.lower(),XY=xy.upper(),
    hv=hv,HV=hv.upper(),
    label=label,freq=freq);

def make_beam_node (beam,pattern,lm,radec0=None,dlm=None,azel=None,pa=None):
  """Makes beam interpolator node for the given filename pattern.""";
  filenames = [];
  labels = re.split('[,\s]',pattern_labels) if pattern_labels else [""];
  freqs = re.split('[,\s]',freq_labels) if freq_labels else [""];
  for label in labels:
    per_label = [];
    for xy in "xy":
      per_xy = [];
      for freq in freqs:
        filename = make_beam_filename(pattern,xy,label,freq);
        if not os.path.exists(filename):
          raise RuntimeError("Can't find beam pattern file %s"%filename);
        per_xy.append(filename);
      per_label.append(per_xy);
    filenames.append(per_label);
  # form up list of children
  children = [ lm ] if radec0 is None else [ lm,radec0 ];
  children += [ dlm or 0,azel or 0];
  if pa is not None:
    children += [ pa ];
  # now make interpolator node
  from . import InterpolatedBeams
  beam << Meq.PyNode(class_name="EMSSPolarBeamInterpolatorNode" if radec0 is None else "EMSSPolarBeamRaDecInterpolatorNode",
                     module_name=__file__,
                     filename=filenames,
                     missing_is_null=False,spline_order=spline_order,verbose=verbose_level or 0,
                     l_beam_offset=l_beam_offset*DEG,m_beam_offset=m_beam_offset*DEG,
                     beam_symmetry=beam_symmetry,
                     normalization_factor=normalization_factor,rotate_xy=rotate_xy,
                     children=children);


def compute_jones (Jones,sources,stations=None,pointing_offsets=None,inspectors=[],label='E',**kw):
  stations = stations or Context.array.stations;
  ns = Jones.Subscope();

  # declare an inspector node for the Jones matrix -- will be defined below
  insp = Jones.scope.inspector(label);
  inspectors += [ insp ];

  radec0 = Context.observation.phase_centre.radec() if interpol_coord is COORD_THETAPHI else None;

  # loop over sources
  for src in sources:
    xy = src.direction.radec() if interpol_coord is COORD_THETAPHI else src.direction.lm();
    # If sky rotation and/or horizon masking and/or pointing offsets are in effect, we have a per-station beam.
    # Otherwise the beam is the same for all stations.
    if sky_rotation or pointing_offsets or horizon_masking:
      for p in stations:
        azel = src.direction.azel(Context.array.xyz(p)) if horizon_masking else None;
        pa = Context.observation.phase_centre.pa(Context.array.xyz(p)) if sky_rotation else None;
        # apply offset in the node (so pointing offsets are interpreted in the azel frame, if rotating)
        make_beam_node(Jones(src,p),filename_pattern,xy,radec0=radec0,
                       dlm=pointing_offsets and pointing_offsets(p),azel=azel,pa=pa);
    else:
      make_beam_node(Jones(src),filename_pattern,xy,radec0=radec0);
      for p in stations:
        Jones(src,p) << Meq.Identity(Jones(src));

  # define an inspector
  if sky_rotation or pointing_offsets or horizon_masking:
    # Jones inspector is per-source, per-station
    insp << StdTrees.define_inspector(Jones,sources,stations,label=label);
  else:
    # Jones inspector is per-source
    insp << StdTrees.define_inspector(Jones,sources,label=label);

  return Jones;


def compute_jones_tensor (Jones,sources,stations,lmn=None,pointing_offsets=None,inspectors=[],label='E',**kw):
  stations = stations or Context.array.stations;
  ns = Jones.Subscope();

  radec0 = Context.observation.phase_centre.radec() if interpol_coord is COORD_THETAPHI else None;
  
  # if using lm-interpolation, make lmn tensor
  if radec0 is None:
    # see if sources have a "beam_lm" or "_lm_ncp" attribute
    lmsrc = [ src.get_attr("beam_lm",None) or src.get_attr("_lm_ncp",None) for src in sources ];
    # if all source have the attribute, create lmn tensor node (and override the lmn argument)
    if all([lm is not None for lm in lmsrc]):
      lmn = ns.lmT << Meq.Constant(lmsrc);
    # failing that (and if lmn was not supplied), make tensor here 
    if not lmn:
      lmn = ns.lmT << Meq.Composer(dims=[0],*[ src.direction.lm() for src in sources ]);
    xy = lmn;
  # else make a radec tensor. First try to make it constant, failing that, use a composer
  else:
    radec = [ src.direction.radec_static() for src in sources ];
    if all(radec):
      xy = ns.radecT << Meq.Constant(radec);
    else:
      xy = ns.radecT << Meq.Composer(dims=[0],*[ src.direction.radec() for src in sources ]);
      
  # if horizon masking is in effect, make an azel tensor
  if horizon_masking:
    azel = ns.azelT << Meq.Composer(dims=[0],*[ src.direction.azel() for src in sources ]);
  else:
    azel = None;

  # create station tensor
  # If sky rotation and/or pointing offsets are in effect, we have a per-station beam.
  # Otherwise the beam is the same for all stations.
  if sky_rotation or pointing_offsets or randomize_rotation or horizon_masking:
    for p in stations:
      angle = None;
      if randomize_rotation:
        angle = ns.rot0(p) << random.uniform(-math.pi,math.pi);
      if sky_rotation:
        pa = Context.observation.phase_centre.pa(Context.array.xyz(p));
        angle = (ns.rot1(p) << angle + pa) if angle is not None else pa;
      # apply offset in the node (so pointing offsets are interpreted in the azel frame, if rotating)
      make_beam_node(Jones(p),filename_pattern,xy,radec0=radec0,dlm=pointing_offsets and pointing_offsets(p),azel=azel,pa=angle);
    ns.inspector1 << Meq.Composer(dims=[len(stations)*len(sources),2,2],*[Jones(p) for p in stations]);
  else:
    make_beam_node(Jones,filename_pattern,xy,radec0=radec0);
    ns.inspector1 << Meq.Composer(Jones,dims=[len(sources),2,2]);
    Jones = lambda p,J=Jones:J;
    
  # add inspector
  ns.inspector << Meq.Mean(ns.inspector1,reduction_axes=["freq"],
      plot_label=["%s:%s"%(src.name,p) for p in stations for src in sources]);
  inspectors.append(ns.inspector);
  
  return Jones;
    
