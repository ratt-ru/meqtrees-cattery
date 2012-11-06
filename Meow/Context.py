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

array = None;
observation = None;
vdm = None;
mssel = None;
correlations = ["XX","XY","YX","YY" ];
active_correlations = ["XX","XY","YX","YY" ];

max_abs_w = 0;
# if this is set to true (i.e. by an imported module), then MSUtils will attempt to
# discover the maximum W coordinate and fill it
discover_max_abs_w = False;

unit_coherency = 1;
# 0.5 corresponds to the ME papers 1-4, and the "old" (up to release 1.1.1) 
# MeqTrees convention
# 1 corresponds to AIPS/TMS/etc, and MeqTrees from 1.2 on

# This functions returns x*unit_coherency. Useful for specifying nodes,
# as it causes useless multiplications by 1 to be omitted,
def unitCoherency (x):
  return x*unit_coherency if unit_coherency != 1 else x;

def set (array=None,observation=None):
  """Sets the global Meow context."""
  for arg in 'array','observation':
    val = locals().get(arg,None);
    if val is not None:
      globals()[arg] = val;
  
def get_array (array1):
  """Resolves 'array1' argument: returns array1 if it is non-false, else
  returns the global array, or throws an error if none is set."""
  arr = array1 or array;
  if not arr:
    raise ValueError,"'array' must be set in global Meow.Context, or supplied explicitly";
  return arr;
  
def get_observation (obs1):
  """Resolves 'obs1' argument: returns obs1 if it is non-false, else
  returns the global observation, or throws an error if none is set."""
  obs = obs1 or observation;
  if not obs:
    raise ValueError,"'observation' must be set in global Meow.Context, or supplied explicitly";
  return obs;
  
def get_dir0 (dir0):
  """Resolves 'dir0' argument: returns dir0 if it is non-false, else
  returns the global observation's phase center, or throws an error 
  if none is set."""
  if dir0:
    return dir0;
  if not observation:
    raise ValueError,"'observation' must be set in global Meow.Context, or a 'dir0' supplied explicitly";
  return observation.phase_centre;



  
