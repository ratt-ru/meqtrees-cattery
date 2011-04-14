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

from SkyComponent import *
import Context
import Parallelization

class Patch (SkyComponent):
  def __init__(self,ns,name,direction,components=[]):
    SkyComponent.__init__(self,ns,name,direction);
    self._components = list(components);

  def add (self,*comps):
    """adds components to patch""";
    self._components += comps;

  def coherency (self,array=None,observation=None,nodes=None,**kw):
    nodes = nodes or self.ns.vis;
    array = Context.get_array(array);
    ifrs = array.ifrs();
    # no components -- use 0
    if not self._components:
      [ nodes(*ifr) << 0.0 for ifr in ifrs ];
    else:
      # each component is either a SkyComponent, else a visibility basenode
      # for up list of (component,visibility_node) pairs, where component=None if component is directly a node
      cv_list = [ (comp,comp.visibilities(array,observation)) if isinstance(comp,SkyComponent) else (None,comp) for comp in self._components ];
      # determine which components are now solvable
      solvables    = [ vis for comp,vis in cv_list if comp and comp.get_solvables() ];
      nonsolvables = [ vis for comp,vis in cv_list if not (comp and comp.get_solvables()) ];
      # if both types are present, add separately for optimum cache reuse
      if solvables and nonsolvables:
        solv = nodes('solv');
        nonsolv = nodes('nonsolv');
        # use the intelligence in Parallelization to add in a smart way, depending on out
        # parallelization settings
        Parallelization.add_visibilities(solv,solvables,ifrs);
        Parallelization.add_visibilities(nonsolv,nonsolvables,ifrs);
        for ifr in ifrs:
          nodes(*ifr) << Meq.Add(solv(*ifr),nonsolv(*ifr));
      else:
        # use the intelligence in Parallelization to add in a smart way, depending on out
        # parallelization settings
        Parallelization.add_visibilities(nodes,solvables+nonsolvables,ifrs);
    return nodes;
