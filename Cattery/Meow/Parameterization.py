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
from Timba.array import *
import Meow

def create_polc(c00=0.0,deg_f=0,deg_t=0):
  """helper function to create a t/f polc with the given c00 coefficient,
  and with given order in t/f""";
  polc = meq.polc(zeros((deg_t+1, deg_f+1))*0.0);
  polc.coeff[0,0] = c00;
  return polc;

# type of polc object  
POLC_TYPE = type(meq.polc(0));

def resolve_parameter (name,node,value,tags=[],solvable=True,solvables=None):
  """Helper function, resolves 'value' to a parameter-like node.
  'name' is a parameter name (only used for error messages)
  'node' is an uninitialized node stub
  'value' specifies the parameter.
      (a) If value is numeric, creates a Meq.Constant() and binds it to 'node'.
      (b) If value is a node, returns it as-is.
      (c) If value is a Meow.Parm specification, creates a Meq.Parm(), and binds it to 'node'. If
          'solvable' is True, create Meq.Parm will have a 'solvable' tag.
  'solvables' may be None, or a list. If it is a list, then any solvable parameters 
      will be appended to this list. In case (c) this is the created Parm (if solvable=True). In
      case (b) a search will be performed on the node to extract any solvable parameters.
  """
  # make sure tags is a list, and add name
  if isinstance(tags,str):
    tags = tags.split(" ");
  else:
    tags = list(tags);
  tags.append(name);
  if isinstance(value,(int,float,complex)):
    return node << Meq.Constant(value=value,tags=tags);
  elif is_node(value):
    if solvables is not None:
      solvables += value.search(tags="solvable");
    return value;
  elif isinstance(value,Meow.Parm):
    if solvable:
      tags.append("solvable");
    node << value.make(tags);
    if solvables is not None and solvable:
      solvables.append(node);
    return node;
  else:
    raise TypeError,"Parameter '"+name+"' can only be defined as a constant, a node, or a Meow.Parm";
  

class Parameterization (object):
  """Parameterization is a base class for objects with parameters.
  It provides services for managing parms, etc.
  It also provides a "qualified scope" object -- available
  as self.ns -- that can be used just like a proper node scope
  to create all nodes associated with this object, automatically
  assigning them the given qualifiers (including the name, if not None).
  The global node scope is available as self.ns0.
  If needed, the set of qualifiers (including name) can be accessed via 
  quals()/kwquals().
  """;
  
  def __init__(self,ns,name,quals=[],kwquals={}):
    if not isinstance(ns,NodeScope):
      raise TypeError,"expected NodeScope, got "+ \
	getattr(ns,'__class__',type(ns)).__name__;
    self.ns0    = ns;
    self.name   = name;
    self._parmnodes = {};
    self._parmdefs = {};
    self._quals     = list(quals);
    self._kwquals   = kwquals;
    if name is not None:
      self._quals.append(name);
    if self._quals or self._kwquals:
      self.ns = ns.QualScope(*self._quals,**self._kwquals);
    else:
      self.ns = ns;
    # this will be a list of Parms that are used for any parameters. 
    # if empty, then all parameters are constants
    self._solvables = [];
    
  def _add_parm (self,name,value,tags=[],solvable=True):
    """Adds an entry for parameter named 'name'. No nodes are created yet;
    they will only be created when self._parm(name) is called later. This
    is useful because a Meow component can pre-define all necessary 
    parameters at construction time with _add_parm(), but then make nodes
    only for the ones that are actually in use with _parm().
    
    'value' can be:
    * a numeric constant, in which case a Meq.Constant is created
    * an actual node, in which case that node is used as is
    * a Meow.Parm, in which case a Meq.Parm is defined, with the given tags
      added on.
    If solvable=True, a "solvable" tag will be added on. This marks 
    potentially solvable parms. This is true by default since everything
    can be solved for; if you want to make a non-solvable parm, use False.
    """;
    if not isinstance(value,(int,float,complex,Meow.Parm)) and \
       not is_node(value):
      raise TypeError,"argument must be a constant, a node, or a Meow.Parm";
    self._parmdefs[name] = (value,tags,solvable);
    
  def get_value (self,name,default=None):
    """Gets default value for named parm, or None if it is a node or a funklet.""";
    if name not in self._parmdefs:
      if default is None:
        raise KeyError,"unknown Meow parm '"+name+"'";
      else:
        return default;
    value,tags,solvable = self._parmdefs[name];
    # now figure out the value
    if isinstance(value,Meow.Parm):
      value = value.value;
    if isinstance(value,(int,float,complex)):
      return value;
    return None;
    
  def _parm (self,name,value=None,tags=[],nodename=None,solvable=True):
    """Returns node representing parameter 'name'.
    If 'nodename' is None, node is named after parm, else
    another name may be given.
    If value is not supplied, then the parameter should have been previously
    defined via _add_parm(). Otherwise, you can define and create a parm 
    on-the-fly by suppying a value and tags as to _add_parm.
    
    """;
    if name in self._parmnodes:
      if value is not None:
        raise TypeError,"duplicate values for Meow parm '"+name+"'";
      return self._parmnodes[name];
    # else has to be defined first (possibly)
    if value is None:
      value,tags,solvable = self._parmdefs.get(name,(None,None,None));
      if value is None:
        raise TypeError,"Meow parm '"+name+"' not previously defined";
    else:
      self._add_parm(name,value,tags,solvable=solvable);
    # now define the node
    nodename = nodename or name;
    parmnode = self._parmnodes[name] = resolve_parameter(name,self.ns[nodename],
                                              value,tags,solvable,self._solvables);
    return parmnode;
    
  def get_solvables (self):
    """Returns list of solvable parameters for this object""";
    return self._solvables;

  def _get_constant (self,name):
    """Returns constant corresponding to given parm, or None if parm is non-constant""";
    value,tags,solvable = self._parmdefs.get(name,(None,None,None));
    if value is None:
      raise TypeError,"Meow parm '"+name+"' not previously defined";
    if isinstance(value,(int,float,complex)):
      return value;
    return None;
  
  def _is_constant (self,name):
    """Returns True if given parm is a constant""";
    value,tags,solvable = self._parmdefs.get(name,(None,None,None));
    return isinstance(value,(int,float,complex));
