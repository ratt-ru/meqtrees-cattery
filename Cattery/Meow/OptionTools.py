# -*- coding: utf-8 -*-
#
#% $Id: Utils.py 6330 2008-09-05 17:39:55Z oms $ 
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
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from Timba.TDL import *
from Timba.Meq import meq
import Meow

class ListOptionParser (object):
  def __init__ (self,minval=None,maxval=None,name=""):
    self.minval,self.maxval = minval,maxval;
    self.name = name or "index";
    
  def set_min (self,minval):
    self.minval = minval;
  def set_max (self,maxval):
    self.maxval = maxval;
    
  def _validate_index (self,index,minval,maxval):
    if self.minval is None:
      minval = min(minval,index);
    if self.maxval is None:
      maxval = max(maxval,index);
    if index < minval or index > maxval:
      raise ValueError("illegal %s '%d'"%(self.name,index));
    return index,minval,maxval;
    
  def apply (self,value,objlist=None,names=None):
    """Parses string list of selection values, returns list of selected objects from objlist.
    If neither 'names' nor 'objlist' is given, just parses and returns numbers.
    If 'names' is given, it should be a separate list of object names.
    'names' can be set to True to use str() to derive am object's name.
    'names' can be set to False to ignore name specifications completely""";
    value = value.strip();
    if not value:
      return objlist;
    minval = self.minval or 0;
    maxval = len(objlist)-1 if objlist else self.maxval;
    # names
    if names:
      if names is True:
        names = [ str(obj) for obj in (objlist or [])];
      elif len(names) != len(objlist):
        raise ValueError("number of names must match number of objects");
      name_to_index = dict([(name,i) for i,name in enumerate(names)]);
    # init set of selected indices. If spec starts with an exclusion ("-smth"),
    # start with full set, else with empty set
    if value.startswith('-'):
      subset = set(range(maxval or 0));
    else:
      subset = set();
    # now parse specs one  by one
    for spec in re.split("[\s,]+",value):
      if not spec:
        continue;
      # process "-" at start
      if spec.startswith("-"):
        oper = subset.difference_update;
        spec = spec[1:];
      else:
        oper = subset.update;
      # standard terms
      if spec.upper() == "ALL":
        oper(list(range((maxval or 0)+1)));
        continue;
      # single number
      match = re.match("^\d+$",spec);
      if match:
        index,minval,maxval = self._validate_index(int(spec),minval,maxval);
        oper([index]);
        continue;
      # [number]:[number]
      match = re.match("^(\d+)?:(\d+)?$",spec);
      if not match:
        if names and spec in name_to_index: 
          oper([name_to_index[spec]]);
        elif names is not False:
          raise ValueError("illegal %s '%s'"%(self.name,spec));
        continue;
      if match.group(1):
        index1 = int(match.group(1));
      elif self.minval is not None:
        index1 = self.minval;
      else:
        raise ValueError("illegal %s '%s'"%(self.name,spec));
      if match.group(2):
        index2 = int(match.group(2));
      elif self.maxval is not None:
        index2 = self.maxval;
      else:
        index2 = index1;
      index1,minval,maxval = self._validate_index(index1,minval,maxval);
      index2,minval,maxval = self._validate_index(index2,minval,maxval);
      # add to subset
      oper(list(range(index1,index2+1)));
    # return subset
    subset = sorted(subset);
    if objlist:
      return [ objlist[i] for i in subset ];
    else:
      return subset;

  def validator (self,value):
    self.apply(value,names=False);
    return True;
