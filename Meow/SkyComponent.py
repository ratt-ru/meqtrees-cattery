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
    
  def enable_smearing (self, smearing=True):
    self.smearing = smearing;
    
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
    
  def make_visibilities (self,nodes,array,observation,**kw):
    """Abstract method.
    Creates nodes computing nominal visibilities of this component 
    Actual nodes are then created as nodes(name,sta1,sta2) for all array.ifrs().
    Returns partially qualified visibility node which must be qualified 
    with an (sta1,sta2) pair.
    """;
    raise TypeError,type(self).__name__+".make_visibilities() not defined";
  
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
    
  def visibilities  (self,array=None,observation=None,nodes=None,**kw):
    """Creates nodes computing visibilities of component.
    'array' is an IfrArray object, or None if the global context is to be used.
    'observation' is an Observation object, or None if the global context is 
      to be used.
    If 'nodes' is None, creates visibility nodes as ns.visibility(name,...,p,q), 
    where '...' are any extra qualifiers from observation.radec0().
    Otherwise 'nodes' is supposed to refer to an unqualified node, and 
    visibility nodes are created as nodes(p,q).
    Returns the actual unqualified visibility node that was created, i.e. 
    either 'nodes' itself, or the automatically named nodes
    """;
    observation = Context.get_observation(observation);
    array = Context.get_array(array);
    if nodes is None:
      nodes = self.ns.visibility.qadd(observation.radec0());
    if not nodes(*(array.ifrs()[0])).initialized():
      self.make_visibilities(nodes,array,observation,**kw);
    return nodes;
    
  def corrupt (self,jones,per_station=True,label=None):
    from Meow.CorruptComponent import CorruptComponent
    if per_station:
      return CorruptComponent(self.ns0,self,station_jones=jones,label=label);
    else:
      return CorruptComponent(self.ns0,self,jones=jones,label=label);
      
