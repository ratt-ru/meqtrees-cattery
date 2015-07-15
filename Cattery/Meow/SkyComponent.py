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

from Timba.TDL import *
from Timba.Meq import meq
from Parameterization import *
from Direction import *
import Context

class SkyComponent (Parameterization):
  """A SkyComponent represents an abstract sky model element.
  SkyComponents have an name and an associated direction.
  """;
  def __init__(self,ns,name,direction):
    Parameterization.__init__(self,ns,name);
    if isinstance(direction,Direction):
      self.direction = direction;
    else:
      if not isinstance(direction,(list,tuple)) or len(direction) != 2:
        raise TypeError,"direction: Direction object or (ra,dec) tuple expected";
      ra,dec = direction;
      self.direction = Direction(ns,name,ra,dec);
    # If the source uses station decomposition (i.e. if the sqrt_visibilities() method has been called),
    # this will be set to True.
    self.using_station_decomposition = False;
    # if source should include a time/bandwidth smearing correction, this will be true
    self.smearing = False;
    # user-defined attributes. May be used for anything.
    self.attrs = {};

  def enable_smearing (self,smearing=True):
    self.smearing = smearing;

  def is_polarized(self):
    return True

  def is_smeared (self):
    return self.smearing;

  def set_attr (self,attr,value):
    self.attrs[attr] = value;

  def get_attr (self,attr,default=None):
    return self.attrs.get(attr,default);

  def radec (self):
    """Returns ra-dec two-pack for this component's direction""";
    return self.direction.radec();

  def lmn (self,dir0=None):
    return self.direction.lmn(dir0);

  def get_solvables (self):
    return Parameterization.get_solvables(self) + self.direction.get_solvables();

  def is_station_decomposable (self):
    """Returns true if source can be decomposed into per-station contributions
    (i.e. if the sqrt_visibilities method below is implemented)""";
    return False;

  def sqrt_visibilities (self,array=None,observation=None,nodes=None):
    """If component can be decomposed into a per-station contribution,
    this creates the per-station nodes. If not, None can be returned.
    'array' is an IfrArray object, or None if the global context is to be used.
    'observation' is an Observation object, or None if the global context is
      to be used.
    If 'nodes' is None, creates sqrt-visibility nodes as ns.sqrtvis(name,...,p),
    where '...' are any extra qualifiers from observation.radec0().
    Otherwise 'nodes' is supposed to refer to an unqualified node, and
    sqrt-visibility nodes are created as nodes(p).
    Returns the actual unqualified sqrt-visibility node that was created, i.e.
    either 'nodes' itself, or the automatically named nodes. Else returns None
    (if decomposition is not supported).""";
    return None;

  def coherency (self,array=None,observation=None,nodes=None,**kw):
    """Returns nodes computing coherencies of component.
    'array' is an IfrArray object, or None if the global context is to be used.
    'observation' is an Observation object, or None if the global context is
      to be used.
    'nodes' is a base node (which will be qualified with p,q). If None,
    a default must be used.
    Returns something that must be qualified with (p,q) to get a coherency node.
    """;
    raise TypeError,type(self).__name__+".coherency() not defined";

  def smear_factor (self,array=None,dir0=None):
    """Returns smearing factor associated with this source and direction.
    Returns something that must be qualified with p,q, or can be None if there's no smearing.
    By default, uses the Direction implementation."""
    return self.direction.smear_factor(array,dir0);

  def visibilities  (self,array=None,observation=None,nodes=None,smear=False,**kw):
    """Creates nodes computing visibilities of component.
    'array' is an IfrArray object, or None if the global context is to be used.
    'observation' is an Observation object, or None if the global context is
      to be used.
    'smear' is True if smearing is to be applied.
    'nodes' is a base node (which will be qualified with p,q). If None,
    a default may be used.
    Returns something that must be qualified with (p,q) to get a visibility node.
    """;
    observation = Context.get_observation(observation);
    cohnodes = self.coherency(array,observation);
    # return coherency directly, if we're at phase centre
    if self.direction is observation.phase_centre:
      return cohnodes;
    # else apply, phase shifts and (optionally) smearing
    visnodes = nodes or self.ns.vis;
    array = Context.get_array(array);
    if not visnodes(*array.ifrs()[0]).initialized():
      # apply smear factor if asked to (and if it's implemented)
      smear = (smear or self.is_smeared()) and \
              self.smear_factor(array,observation.phase_centre);
      if smear:
        cohsm = visnodes('smear');
        for p,q in array.ifrs():
          cohsm(p,q) << smear(p,q)*cohnodes(p,q);
        cohnodes = cohsm;
      # phase shift
      self.direction.make_phase_shift(visnodes,cohnodes,array,
                                      Context.get_observation(observation).phase_center);
    return visnodes;

  def corrupt (self,jones,per_station=True,label=None):
    from Meow.CorruptComponent import CorruptComponent
    if per_station:
      return CorruptComponent(self.ns0,self,station_jones=jones,label=label);
    else:
      return CorruptComponent(self.ns0,self,jones=jones,label=label);
