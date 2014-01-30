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
from SkyComponent import *

class CorruptComponent(SkyComponent):
  """A CorruptComponent represents an SkyComponent, plus a set of
  associated Jones matrices that corrupt the SkyComponent's visibilities.
  """;
  def __init__(self,ns,skycomp,label=None,station_jones=None,jones=None):
    """Initializes a corrupt component. skycomp is a SkyComponent
    object.
    'label' is used to qualify visibilities, this needs to be unique
    (for the given corruption at least). If not supplied, the jones or
    station_jones name is used.
    """;
    if label is None:
      if jones:
        label = getattr(jones,'name','corrupt');
      elif station_jones:
        # bit of a dirty trick here -- qualify station_jones with 0,
        # which doesn't even have to be a real station, to get at the
        # actual node. This allows us to support generic callables
        # for station_jones
        label = getattr(station_jones(0),'basename','corrupt');
      else:
        label = "corrupt";
    SkyComponent.__init__(self,ns,skycomp.name+':'+label,skycomp.direction);
    self.label    = label;
    self.skycomp  = skycomp;
    self._jones = [];
    if station_jones is not None:
      self.add_station_jones(station_jones);
    if jones is not None:
      self.add_jones(jones);

  def is_smeared (self):
    return self.skycomp.is_smeared();

  def add_station_jones (self,jones,prepend=False):
    """adds a per-station image-plane effect represented by a Jones
    matrix. 'jones' should be a callable object (e.g. an unqualified node)
    such that jones(x) returns the Jones matrix node for station x, or None
    if no matrix is to be applied""";
    # prepend to list
    if prepend:
      self._jones.insert(0,jones);
    else:
      self._jones.append(jones);

  def add_jones (self,jones,prepend=False):
    """adds an station-independent image-plane effect represented by the
    given Jones matrix. Argument should be an valid node.""";
    # emulate a per-station jones so that we may jones matrices uniformly
    # in visibility()
    self.add_station_jones(lambda sta : jones,prepend);

  def jones_list (self):
    return self._jones;

  def apply_jones (self,vis,vis0,ifr_list):
    """Creates nodes to apply the Jones chain associated with this
    component.
    'ifr_list' should be a list of (sta1,sta2) pairs
    'vis' is the output node which will be qualified with (sta1,sta2)
    'vis0' is an input visibility node which will be qualified with (sta1,sta2)
    """;
    # multiply input visibilities by our jones list
    for (sta1,sta2) in ifr_list:
      # collect list of per-source station-qualified Jones terms
      terms = [ jones(sta1) for jones in self.jones_list() if jones(sta1) is not None ];
      # reverse list since they are applied in reverse order
      # first (J2*J1*C*...)
      terms.reverse();
      terms.append(vis0(sta1,sta2));
      # collect list of conjugate terms. The '**' operator
      # is for init-if-not-initialized
      terms += [ jones(sta2)('conj') ** Meq.ConjTranspose(jones(sta2))
                 for jones in self.jones_list() if jones(sta2) is not None ];
      # create multiplication node
      vis(sta1,sta2) << Meq.MatrixMultiply(*terms);
    return vis;

  def coherency (self,array=None,observation=None,nodes=None,**kw):
    coh0 = self.skycomp.coherency(array,observation,**kw);
    coh = nodes or (coh0(self.label) if is_node(coh0) else self.ns.coh);
    ifrs = Context.get_array(array).ifrs();
    # do we have extra jones terms?
    if self.jones_list():
      if not coh(*ifrs[0]).initialized():
        self.apply_jones(coh,coh0,ifrs);
    # no jones terms, use nominal visibilities directly
    else:
      for ifr in ifrs:
        coh(*ifr) << Meq.Identity(coh0(*ifr));
    return coh;
