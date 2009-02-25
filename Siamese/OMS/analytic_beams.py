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

from Timba.TDL import *
from Meow.Direction import radec_to_lmn
import Meow
from Meow import Context
from Meow import StdTrees
from Meow import ParmGroup

import math
from math import sqrt,atan2
import sets

DEG = math.pi/180.;
ARCMIN = DEG/60;
ARCSEC = DEG/3600;

class WSRT_cos3_beam (object):
  label = "wsrt_cos3";
  name  = "WSRT cos^3 beam model";
  
  def __init__ (self):
    self._options = [
      TDLOption('bf',"Beam factor (1/GHz)",[65.],more=float,namespace=self),
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
      TDLOption('dish_sizes',"Dish size(s)",[25],more=str,namespace=self)
    ];
    self._prepared = False;
    
  def option_list (self):
    return self._options;
  
  def prepare (self):
    if self._prepared:
      return;
    try:
      self._sizes = map(float,str(self.dish_sizes).split());
    except:
      raise RuntimeError,"illegal dish sizes specification";
    self._per_station = len(self._sizes)>1;
    self._prepared = True;
  
  def is_per_station (self):
    self.prepare();
    return self._per_station;

  def compute (self,E,lm,p=0):
    """computes a cos^3 beam for the given direction, using NEWSTAR's
    cos^^(fq*B*scale*r) model (thus giving a cos^6 power beam).
    r=sqrt(l^2+m^2), which is not entirely accurate (in terms of angular distance), but close enough 
    to be good for NEWSTAR, so good enough for us.
    'E' is output node
    'lm' is direction (2-vector node, or l,m tuple)
    'p' is station number (only if is_per_station() returned True);
    """
    ns = E.Subscope();
    if isinstance(lm,(list,tuple)):
      l,m = lm;
      r = ns.r << sqrt(l*l+m*m);
    else:
      r = ns.r << Meq.Norm(lm);
    clip = self.clip;
    if self.newstar_clip:
      clip = -clip;
    size = self._sizes[min(p,len(self._sizes)-1)];
    E << Meq.WSRTCos3Beam(self.bf*(size/25.)*1e-9,r,clip=clip);
    return E;

class SkyJonesNotPerStation (object):
  """This class implements a Sky Jones contract, given a set of nodes that are
  per-source but (possibly) not per-station""";
  def __init__ (self,jones,p0):
    self.jones = jones;
    # p0 is the first station index. It is used to check if a beam is per-station
    self.p0 = p0;
  def __call__ (self,src,p=None):
    beam = self.jones(src);
    # if station is specified, return per-station beam (if available),
    # else the common source beam
    if p is not None:
      if beam(p).initialized():
        return beam(p);
      else:
        return beam;
    # else return potentially qualifiable beam
    else:
      if beam(self.p0).initialized():
        return beam;
      else:
        return lambda p:beam;
  def search (self,*args,**kw):
    return self.jones.search(*args,**kw);

def compute_jones (Jones,sources,stations=None,
                    label="beam",pointing_offsets=None,inspectors=[],
                    solvable_sources=sets.Set(),**kw):
  """Computes beam gain for a list of sources.
  The output node, will be qualified with either a source only, or a source/station pair
  """;
  stations = stations or Context.array.stations();
  ns = Jones.Subscope();
  
  # figure out beam model
  for beam_model in MODELS:
    if globals().get('use_%s'%model.label,None):
      break;
  else:
    raise RuntimeError,"no beam model selected";
  beam_model.prepare();
  
  # this dict will hold LM tuples (or nodes) for each source.
  lmsrc = {};
  # see if sources have a "beam_lm" attribute, use that for beam offsets
  for src in sources:
    lm = src.get_attr("beam_lm",None);
    if lm:
      l,m = lmsrc[src.name] = lm;
      src.set_attr(label+'r',math.sqrt(l**2+m**2)/math.pi*(180*60));
    # else try to use static lm coordinates
    else:
      # else try to use static lm coordinates
      lmnst = src.direction.lmn_static();
      if lmnst:
        l,m = lmsrc[src.name] = lmnst[0:2];
        src.set_attr(label+'r',math.sqrt(l**2+m**2)/math.pi*(180*60));
      # else use lmn node
      else:
        lmsrc[src.name] = src.direction.lm();
    
  # set of all sources, will later become set of sources for which
  # a beam hasn't been computed
  all_sources = sets.Set([src.name for src in sources]);
  
  # if pointing errors are enabled, compute beams for sources with a PE
  if pointing_offsets:
    # is a subset specified?
    if pe_subset != "all":
      pe_sources = sets.Set(pe_subset.split(" "));
      # make sure pe_sources does not contain unknown sources
      pe_sources &= all_sources;
      # remove PE_enabled sources from global list
      all_sources -= pe_sources;
    # else all sources have a PE
    else:
      pe_sources = all_sources;
      all_sources = [];
    # create nodes to compute actual pointing per source, per antenna
    for name in pe_sources:
      solvable_sources.add(name);
      lm = lmsrc[name];
      # if LM is a static constant still, make constant node for it
      if isinstance(lm,(list,tuple)):
        lm = ns.lm(name) << Meq.Constant(value=Timba.array.array(lm));
      for ip,p in enumerate(Context.array.stations()):
        # make offset lm
        lm1 = lm(p) << lm + pointing_offsets(p);
        # make beam model
        beam_model.compute(Jones(name,p),lm1,ip);
  
  # now, all_sources is the set of sources for which we haven't computed the beam yet
  if all_sources:
    # if not per-station, use same beam for every source
    if not beam_model.is_per_station():
      for name in all_sources:
        beam_model.compute(Jones(name),lmsrc[name]);
      inspectors.append(Jones('inspector') << StdTrees.define_inspector(Jones,list(all_sources)));
      # in this case we return a wrapper around the base Jones node, to make sure
      # that the station index is applied or not applied depending on whether a per-station
      # beam is available
      return SkyJonesNotPerStation(Jones,p0=Context.array.stations()[0]);
    # else make per-source, per-station beam
    else:
      for name in all_sources:
        for ip,p in enumerate(Context.array.stations()):
          beam_model.compute(Jones(name,p),lmsrc[name],ip);
          
  return Jones;

# this will be set to True on a per-source basis, or if beam scale is solvable
solvable = False;

# available beam models
MODELS = [
  WSRT_cos3_beam()
];


_model_option = TDLCompileMenu("Beam model",
  exclusive='model_type',
  *[ TDLMenu(model.name,toggle="use_%s"%model.label,*model.option_list())
      for model in MODELS ]);

TDLCompileOption('pe_subset',"Apply pointing errors (if any) to a subset of sources",
        ["all"],more=str,doc="""Selects a subset of sources to which pointing errors 
        (if any) are applied. Enter soure names separated by space""");


