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
"""<P>This implements primary beams expressed as analytic functions. For now, only the WSRT cos^3 beam 
model is provided, but readers are encouraged to contribute models of their own.</P>

<P>See also Calico.OMS.wsrt_beams for a WSRT beam model with solvable pointing and scale.</P> 

<P>Author: O. Smirnov</P>""";
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division
from Timba.TDL import *
from .AnalyticBeams import ClippedSincBeam
from Meow.Direction import radec_to_lmn
import Meow
from Meow import Context
from Meow import StdTrees
from Meow import ParmGroup

import math
from math import sqrt,atan2

DEG = math.pi/180.;
ARCMIN = DEG/60.;
ARCSEC = DEG/3600.;

class WSRT_cos3_beam (object):
  label = "wsrt_cos3";
  name  = "WSRT cos^3 beam model";
  
  def __init__ (self):
    self._options = [
      TDLOption('bf',"Beam factor (1/GHz)",[65.],more=float,namespace=self),
      TDLOption('ell',"Ellipticity",[0,0.01],more=float,namespace=self,
        doc="""<P>This makes the X and Y beam patterns elliptical (with major and minor axes of 1&#177;<I>e</I>) rather
        than circular. Positive <I>e</I> corresponds to the X dipole beam being elongated in the East-West direction,
        and the Y beam being elongated North-South.</P>"""),
      TDLOption('clip',"Clip beam once gain goes below",[.1],more=float,namespace=self,
        doc="""This makes the beam flat outside the circle corresponding to the specified gain value."""),
      TDLOption('newstar_clip',"Use NEWSTAR-compatible beam clipping (very naughty)",
        False,namespace=self,
        doc="""NEWSTAR's primary beam model does not implement clipping properly: it only makes the
        beam flat where the model gain goes below the threshold. Further away, where the gain (which
        is cos^3(r) essentially) rises above the threshold again, the beam is no longer clipped, and so 
        the gain goes up again. This is a physically incorrect model, but you may need to enable this
        if working with NEWSTAR-derived models, since fluxes of far-off sources will have been estimated with
        this incorrect beam."""),
      TDLOption('dish_sizes',"Dish size(s), m.",[25],more=str,namespace=self)
    ];
    self._prepared = False;
    
  def option_list (self):
    return self._options;
  
  def prepare (self,ns):
    if not self._prepared:
      self.ns = ns;
      try:
        self._sizes = list(map(float,str(self.dish_sizes).split()));
      except:
        raise RuntimeError("illegal dish sizes specification");
      self._per_station = len(self._sizes)>1;
      if self.ell != 0:
        self._ellnode = ns.ell << Meq.Constant(value=Timba.array.array([self.ell,-self.ell]));
      else:
        self._ellnode = None;
      self._prepared = True;

  def compute (self,E,lm,pointing=None,p=0):
    """computes a cos^3 beam for the given direction, using NEWSTAR's
    cos^^(fq*B*scale*r) model (thus giving a cos^6 power beam).
    r=sqrt(l^2+m^2), which is not entirely accurate (in terms of angular distance), but close enough 
    to be good for NEWSTAR, so good enough for us.
    'E' is output node
    'lm' is direction (2-vector node, or l,m tuple)
    'p' is station number (only if is_per_station() returned True);
    """
    # work out beam arguments
    clip = -self.clip if self.newstar_clip else self.clip;
    size = self.bf*1e-9*self._sizes[min(p,len(self._sizes)-1)]/25.;
    ns = E.Subscope();
    # pointing
    if pointing is None:
      pointing = self.ns.dlm_0 << Meq.Constant([0,0]);
    # different invocations for elliptical and non-elliptical beams
    if self._ellnode is not None:
      E << Meq.WSRTCos3Beam(size,lm,pointing,self._ellnode,clip=clip);
    else:
      E << Meq.WSRTCos3Beam(size,lm,pointing,clip=clip);
    return E;
    
  def compute_tensor (self,E,lm,pointing=None,p=0):
    """computes a cos^3 beam for the given lm direction tensor
    'E' is output node
    'lm' is direction (node returning Nx2 tensor)
    'p' is station number
    """
    # fortunately, exactly the same code works for the tensor case
    return self.compute(E,lm,pointing,p);

INV_C = 1.0 / 299792458

class circular_aperture_beam (WSRT_cos3_beam):
  label = "circular_aperture";
  name  = "Circular aperture beam model";
  
  def __init__ (self):
    self._options = [
      TDLOption('bf',"Beam factor",[1.0],more=float,namespace=self),
      TDLOption('dish_sizes',"Dish size(s)",[25],more=str,namespace=self)
    ];
    self._prepared = False;
    self.ell = 0;
    
  def compute (self,E,lm,pointing=None,p=0):
    """Circular Apertures. With equal dishes, should generate Airy disks"""
    ns = E.Subscope();
    if pointing == None:
      pointing = Meq.Constant([0,0])
    if isinstance(lm,(list,tuple)):
      lm = ns.lm1 << Meq.Constant(lm)
    r = ns.r << Meq.Norm(lm-pointing);
    size = self._sizes[min(p,len(self._sizes)-1)];
    # The Max() is to avoid div by zero cases.
    ns.theta_r << Meq.Max(INV_C * r * size * math.pi * self.bf * Meq.Freq(), 1e-9);
    ns.E_theta << 2 * Meq.Bessel(ns.theta_r, order=1) / ns.theta_r;
    E << Meq.Matrix22(ns.E_theta, 0, 0, ns.E_theta);
    return E;

  compute_tensor = None;
#  def compute_tensor (self,E,lm,pointing=None,p=0):
#    return self.compute(E,lm,pointing,p);

class clipped_sinc_beam (WSRT_cos3_beam):
  label = "clipped_sinc";
  name  = "Clipped sinc function beam model";
  
  def __init__ (self):
    self._options = [
      TDLOption('bf',"Beam factor",[1.],more=float,namespace=self),
      TDLOption('dish_sizes',"Dish size(s)",[25],more=str,namespace=self)
    ];
    self._prepared = False;
    self.ell = 0;
    
  def compute (self,E,lm,pointing=None,p=0):
    ns = E.Subscope();
    size = self._sizes[min(p,len(self._sizes)-1)];
    # beamscale is determined by beam factor, plus frequency (divide by 1.4e+9 so that
    # beamscale is exactly 265.667 at 1.4 GHz), plus dish size
    scale = 265.667*self.bf*(1/1.4e+9)*(size/25);
    children = [lm,pointing] if pointing else [lm];
    E << Meq.PyNode(scale=scale,
      class_name="ClippedSincBeam",
      module_name=ClippedSincBeam.__file__,
      *children);
    return E;

  def compute_tensor (self,E,lm,pointing=None,p=0):
    return self.compute(E,lm,pointing,p);


def compute_jones (Jones,sources,stations=None,
                    label="beam",pointing_offsets=None,inspectors=[],
                    solvable_sources=set(),**kw):
  """Computes beam gain for a list of sources.
  The output node, will be qualified with either a source only, or a source/station pair
  """;
  stations = stations or Context.array.stations();
  ns = Jones.Subscope();
  
  # figure out beam model
  for beam_model in MODELS:
    if globals().get('use_%s'%beam_model.label,None):
      break;
  else:
    raise RuntimeError("no beam model selected");
  beam_model.prepare(ns);
  
  # this dict will hold LM tuples (or nodes) for each source.
  lmsrc = {};
  # see if sources have a "beam_lm" attribute, use that for beam offsets
  for src in sources:
    lm = src.get_attr("beam_lm",None) or src.get_attr("_lm_ncp",None);
    if lm:
      src.set_attr(label+'r',math.sqrt(lm[0]**2+lm[1]**2)/math.pi*(180*60));
      lmsrc[src.name] = ns.lm(src) << Meq.Constant(value=Timba.array.array(lm));
    # else try to use static lm coordinates
    else:
      # else try to use static lm coordinates
      lmnst = src.direction.lmn_static();
      if lmnst:
        lm = lmnst[0:2];
        src.set_attr(label+'r',math.sqrt(lm[0]**2+lm[1]**2)/math.pi*(180*60));
        lmsrc[src.name] = ns.lm(src) << Meq.Constant(value=Timba.array.array(lm));
      # else use lmn node
      else:
        lmsrc[src.name] = src.direction.lm();
  
  for src in sources:
    for ip,p in enumerate(stations):
      beam_model.compute(Jones(src,p),lmsrc[src.name],pointing_offsets and pointing_offsets(p),ip);
          
  return Jones;



def compute_jones_tensor (Jones,sources,stations=None,lmn=None,
                          label="beam",pointing_offsets=None,inspectors=[],
                          **kw):
  """Computes beam gain tensor for a list of sources.
  The output node, will be qualified with either a source only, or a source/station pair
  """;
  stations = stations or Context.array.stations();
  ns = Jones.Subscope();
  
  # figure out beam model
  for beam_model in MODELS:
    if globals().get('use_%s'%beam_model.label,None):
      break;
  else:
    raise RuntimeError("no beam model selected");
  
  if getattr(beam_model,'compute_tensor',None) is None:
    return None;
    
  beam_model.prepare(ns);
  
  # see if sources have a "beam_lm" or "_lm_ncp" attribute
  lmsrc = [ src.get_attr("beam_lm",None) or src.get_attr("_lm_ncp",None) for src in sources ];
  
  # if all source have the attribute, create lmn tensor node (and override the lmn argument)
  if all([lm is not None for lm in lmsrc]):
    lmn = ns.lmnT << Meq.Constant(lmsrc);
    
  # if lmn tensor is not set for us, create a composer
  if lmn is None:
    lmn = ns.lmnT << Meq.Composer(dims=[0],*[ src.direction.lm() for src in sources ]);

  # create station tensors
  for ip,p in enumerate(stations):
    beam_model.compute_tensor(Jones(p),lmn,pointing_offsets and pointing_offsets(p),ip);
          
  return Jones;


# available beam models
MODELS = [
  WSRT_cos3_beam(),circular_aperture_beam(),clipped_sinc_beam()
];

_model_option = TDLCompileMenu("Beam model",
  exclusive='model_type',
  *[ TDLMenu(model.name,toggle="use_%s"%model.label,*model.option_list())
      for model in MODELS ]);

