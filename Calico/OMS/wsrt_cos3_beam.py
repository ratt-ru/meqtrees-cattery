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
"""<P>This module implements a WSRT beam E-Jones, with optional solvable parameters.</P>
<P align="right">Author: O. Smirnov &lt;<tt>smirnov@astron.nl</tt>&gt;</P>""";

__default_label__ = "E";
__default_name__  = "WSRT cos^3 beam";

from Timba.TDL import *
from Meow.Direction import radec_to_lmn
import Meow
from Meow import Context
from Meow import StdTrees
from Meow import ParmGroup

import math
from math import sqrt,atan2

import solvable_pointing_errors

DEG = math.pi/180.;
ARCMIN = DEG/60;
ARCSEC = DEG/3600;

PER_ARRAY = "entire array"
PER_STATION = "per station";

class WSRTCos3Beam (object):
  def __init__ (self,label=__default_label__,solvable=False):
    self.tdloption_namespace = label+"_wsrt_beam";
    self.label = label;
    self._options = [
      TDLOption('bf',"Beam scale factor (1/GHz)",[65.],more=float,namespace=self),
      TDLOption('ellipticity',"Beam ellipticity",[0.],more=float,namespace=self,
              doc="""<P>This makes the X and Y beam patterns elliptical (with major and minor axes of 1&#177;<I>e</I>) rather
              than circular. Positive <I>e</I> corresponds to the X dipole beam being elongated in the East-West direction,
              and the Y beam being elongated North-South.</P>"""),
      TDLOption('beam_clip',"Clip beam once gain goes below",[.1],more=float,namespace=self,
        doc="""This makes the beam flat outside the circle corresponding to the specified gain value."""),
      TDLOption('newstar_mode',
        "Use NEWSTAR-compatible beam clipping (very naughty)",False,namespace=self,
        doc="""NEWSTAR's primary beam model does not implement clipping properly: it only makes the
        beam flat where the model gain goes below the threshold. Further away, where the gain (which
        is cos^3(r) essentially) rises above the threshold again, the beam is no longer clipped, and so
        the gain goes up again. This is a physically incorrect model, but you may need to enable this
        if working with NEWSTAR-derived models, since fluxes of far-off sources will have been estimated with this incorrect beam.""")
    ];
    if solvable:
      self._options += [
        TDLOption("solve_pointings","Solvable pointing offsets",False,namespace=self),
        TDLOption("solve_scale","Solvable beam scale",[None,PER_ARRAY,PER_STATION],namespace=self),
        TDLOption("solve_ell","Solvable ellipticity",[None,PER_ARRAY,PER_STATION],namespace=self),
      ];
    else:
      self.solve_pointings = self.solve_scale = self.solve_ell = None;
    self.beamshape = None;
    
  def compile_options (self):
    return self._options;
    
  # initializes parm nodes and other solvables common to both the tensor and matrix 
  # implementations
  def init_parameters (self,ns,sources,stations,inspectors=[]):
    if self.beamshape is not None:
      return;
    # create solvables if enabled
    parms = [];
    parmgroups = [];
    parmdef_scale = Meq.Parm(self.bf*1e-9,tags="beam solvable");
    parmdef_ell = Meq.Parm(self.ellipticity,tags="beam solvable");
    self.per_station = False;

    # solvable beam scale
    if self.solve_scale is PER_STATION:
      self.beamshape = ns.beamshape;
      for p in stations:
        parms.append(beamshape(p) << parmdef_scale);
      inspectors.append(ns.inspector("scale") << 
            StdTrees.define_inspector(ns.beamshape,stations,label=self.label));
      self.per_station = True;
    elif self.solve_scale is PER_ARRAY:
      parms.append(ns.beamshape << parmdef_scale);
      self.beamshape = lambda p:ns.beamshape;
    else:
      ns.beamshape ** (self.bf*1e-9);
      self.beamshape = lambda p:ns.beamshape;
    
    # solvable ellipticities
    if self.solve_ell is PER_STATION:
      self.ell = ns.ell;
      for p in stations:
        ell_xy = ns.ell_xy(p) << parmdef_ell;
        parms.append(ell_xy);
        self.ell(p) << Meq.Composer(ell_xy,-ell_xy);
      inspectors.append(ns.inspector("ellipticity") << 
            StdTrees.define_inspector(ns.ell_xy,stations,label=self.label));
      self.per_station = True;
    elif self.solve_ell is PER_ARRAY:
      ell_xy = ns.ell_xy << parmdef_ell;
      parms.append(ell_xy);
      ns.ell << Meq.Composer(ell_xy,-ell_xy);
      self.ell = lambda p:ns.ell;
    elif self.ellipticity != 0:
      ns.ell ** Meq.Constant([self.ellipticity,-self.ellipticity]);
      self.ell = lambda p:ns.ell;
    else:
      self.ell = lambda p:None;
      
    # make parm group, if any solvables have been created
    if parms:
      parmgroups.append(ParmGroup.Subgroup("beam shape",list(parms)));
      
    # solvable pointings
    if self.solve_pointings:
      self.dlm = ns.dlm;
      parmdef_0 = Meq.Parm(0,tags="beam solvable");
      pparms = [];
      for p in stations:
        dl = ns.dl(p) << parmdef_0;
        dm = ns.dm(p) << parmdef_0;
        ns.dlm(p) << Meq.Composer(dl,dm);
        pparms += [dl,dm];
      parmgroups.append(ParmGroup.Subgroup("pointing offsets",pparms));
      parms += pparms;
      inspectors.append(ns.inspector('dlm') << 
              StdTrees.define_inspector(ns.dlm,stations,label=self.label));
      self.per_station = True;
    else:
      ns.dlm_null ** Meq.Constant([0,0]);
      self.dlm = lambda p:ns.dlm_null;

    # add solve jobs
    self.solvable = bool(parms);
    if self.solvable:
      self.pg_beam = ParmGroup.ParmGroup(self.label,parms,subgroups=parmgroups,table_name="%s.fmep"%self.label,bookmark=True);
      ParmGroup.SolveJob("cal_"+self.label,"Calibrate beam parameters",self.pg_beam);

  def make_beam_nodes (self,E,bf,lm,ell=None,dlm=None):
    """computes a cos^3 beam for the given direction, using NEWSTAR's
    cos^^(fq*B*r) model (thus giving a cos^6 power beam).
    r=sqrt(l^2+m^2), which is not entirely accurate (in terms of angular distance), but close enough
    to be good for NEWSTAR, so good enough for us.
    'E' is output node
    'lm' is direction (2-vector node, or l,m tuple)
    """
    clip = -self.beam_clip if self.newstar_mode else self.beam_clip;
    if ell is None:
      if dlm is None:
        E << Meq.WSRTCos3Beam(bf,lm,clip=clip);
      else:
        E << Meq.WSRTCos3Beam(bf,lm,dlm,clip=clip);
    else:
      E << Meq.WSRTCos3Beam(bf,lm,dlm,ell,clip=clip);
    return E;

  def compute_jones (self,Jones,sources,stations=None,inspectors=[],**kw):
    """Computes beam gain for a list of sources.
    The output node, will be qualified with either a source only, or a source/station pair
    """;
    stations = stations or Context.array.stations();
    ns = Jones.Subscope();

    # init solvables etc.
    self.init_parameters(ns,sources,stations,inspectors);

    # this dict will hold LM tuples (or nodes) for each source.
    lmsrc = {};
    # see if sources have a "beam_lm" attribute, use that for beam offsets
    for src in sources:
      lm = src.get_attr("beam_lm",None) or src.get_attr("_lm_ncp",None);
      if lm:
        src.set_attr(self.label+'r',math.sqrt(lm[0]**2+lm[1]**2)/math.pi*(180*60));
        lmsrc[src.name] = ns.lm(src) << Meq.Constant(value=Timba.array.array(lm));
      # else try to use static lm coordinates
      else:
        # else try to use static lm coordinates
        lmnst = src.direction.lmn_static();
        if lmnst:
          lm = lmnst[0:2];
          src.set_attr(self.label+'r',math.sqrt(lm[0]**2+lm[1]**2)/math.pi*(180*60));
          lmsrc[src.name] = ns.lm(src) << Meq.Constant(value=Timba.array.array(lm));
        # else use lmn node
        else:
          lmsrc[src.name] = src.direction.lm();

    if self.per_station:
      for src in sources:
        for p in stations:
          self.make_beam_nodes(Jones(src,p),self.beamshape(p),lmsrc[src.name],self.ell(p),self.dlm(p));
    else:
      p0 = stations[0];
      for src in sources:
        self.make_beam_nodes(Jones(src,p0),self.beamshape(p0),lmsrc[src.name],self.ell(p0),self.dlm(p0));
        for p in stations[1:]:
          Jones(src,p) << Meq.Identity(Jones(src,p0));

    # make inspectors
    inspectors.append(ns.inspector << StdTrees.define_inspector(
              Jones,sources,stations,label=self.label));

    return Jones;


  def compute_jones_tensor (self,Jones,sources,stations=None,lmn=None,inspectors=[],**kw):
    """Computes beam gain tensor for a list of sources.
    The output node will be qualified with a station.""";

    # if beam is solvable, but only for a subset of sources, then  
    stations = stations or Context.array.stations();
    ns = Jones.Subscope();
    
    # init solvables etc.
    self.init_parameters(ns,sources,stations,inspectors);
    
    # see if sources have a "beam_lm" or "_lm_ncp" attribute
    lmsrc = [ src.get_attr("beam_lm",None) or src.get_attr("_lm_ncp",None) for src in sources ];
    
    # if all source have the attribute, create lmn tensor node (and override the lmn argument)
    if all([lm is not None for lm in lmsrc]):
      lmn = ns.lmnT << Meq.Constant(lmsrc);
      
    # if lmn tensor is not set for us, create a composer
    if lmn is None:
      lmn = ns.lmnT << Meq.Composer(dims=[0],*[ src.direction.lm() for src in sources ]);

    # create station tensors
    if self.per_station:
      for p in stations:
        self.make_beam_nodes(Jones(p),self.beamshape(p),lmn,self.ell(p),self.dlm(p));
      retval = Jones;
    else:
      p0 = stations[0];
      self.make_beam_nodes(Jones,self.beamshape(p0),lmn,self.ell(p0),self.dlm(p0));
      retval = lambda p,J=Jones:J;
            
    return retval;

