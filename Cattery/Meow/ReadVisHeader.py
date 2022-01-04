#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Updates various nodes in the tree based on the visibility
# header. Uses the naming convetion from Contrib.OMS

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
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from Timba.meqkernel import set_state

def process_vis_header (hdr):
  """handler for the visheader""";
  try:
    # phase center
    (ra0,dec0) = hdr.phase_ref;
    # print '[ReadVisHeader] phase centre: ',ra0,dec0;
    try:
      set_state('ra',value=ra0);
      set_state('dec',value=dec0);
    except: pass;
    # antenna positions
    pos = hdr.antenna_pos;
    if pos.ndim != 2 or pos.shape[0] != 3:
	    raise ValueError('incorrectly shaped antenna_pos');
    nant = pos.shape[1];
    coords = ('x','y','z');
    for iant in range(nant):
      sn = ':num'+str(iant);
      # since some antennas may be missing from the tree, ignore errors
      try:
          for (j,label) in enumerate(coords):
              # print '[ReadVisHeader] ',label+sn, 'value = ',pos[j,iant]
              set_state(label+sn,value=pos[j,iant]);
      except: pass;
      # also set nodes for the old naming convention
      sn = ":" + str(iant+1);
      try:
          for (j,label) in enumerate(coords):
              # print '[ReadVisHeader] ',label+sn, 'value = ',pos[j,iant]
              set_state(label+sn,value=pos[j,iant]);
      except: pass;
    # array reference position
    try:
      for (j,label) in enumerate(coords):
	      set_state(label+'0',value=pos[j,0]);
    except: pass;
    # time extent
    (t0,t1) = hdr.time_extent;
    # print '[ReadVisHeader] time extent: ',t0,t1;
    try:
      set_state('time0',value=t0);
      set_state('time1',value=t1);
    except: pass;
    # freq range
    if isinstance(hdr.channel_freq,(int,float)):
      f0=f1 = hdr.channel_freq;
    else:
      f0,f1 = hdr.channel_freq[0],hdr.channel_freq[-1];
    # print '[ReadVisHeader] freq range: ',f0,f1;
    try:
      set_state('freq0',value=f0);
      set_state('freq1',value=f1);
    except: pass;
  except:
    print(traceback.print_exc());
    raise;
    
