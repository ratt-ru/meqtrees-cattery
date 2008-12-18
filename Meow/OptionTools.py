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


from Timba.TDL import *
from Timba.Meq import meq
import Meow

class ListOptionParser (object):
  def __init__ (self,minval=None,maxval=None,name="",allow_names=False):
    self.minval,self.maxval = minval,maxval;
    self.name = name or "index";
    self.allow_names = allow_names;
    
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
      raise ValueError,"illegal %s '%d'"%(self.name,index);
    return index,minval,maxval;
    
  def parse_list (self,value):
    """Parses string list of values, returns list of numbers""";
    if not value:
      return None;
    if self.minval is not None:
      minval = self.minval;
    else:
      minval = 0;
    if self.maxval is not None:
      maxval = self.maxval;
    else:
      maxval = 0;
    subset = [];
    names = [];
    for spec in re.split("[\s,]+",value):
      if spec:
        # single number
        match = re.match("^\d+$",spec);
        if match:
          index,minval,maxval = self._validate_index(int(spec),minval,maxval);
          subset.append(index);
          continue;
        # [number]:[number]
        match = re.match("^(\d+)?:(\d+)?$",spec);
        if not match:
          if self.allow_names:
            names.append(spec);
            continue;
          else:
            raise ValueError,"illegal %s '%s'"%(self.name,spec);
        if match.group(1):
          index1 = int(match.group(1));
        elif self.minval is not None:
          index1 = self.minval;
        else:
          raise ValueError,"illegal %s '%s'"%(self.name,spec);
        if match.group(2):
          index2 = int(match.group(2));
        elif self.maxval is not None:
          index2 = self.maxval;
        else:
          index2 = index1;
        index1,minval,maxval = self._validate_index(index1,minval,maxval);
        index2,minval,maxval = self._validate_index(index2,minval,maxval);
        # add to subset
        subset += range(index1,index2+1);
    if self.allow_names:
      return subset,names;
    else:
      return subset;

  def validator (self,value):
    self.parse_list(value);
    return True;
