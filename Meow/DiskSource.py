#
#% $Id: DiskSource.py 5418 2007-07-19 16:49:13Z sarod $
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

import math
from Timba.TDL import *
from Timba.Meq import meq
from PointSource import *
import Context

STOKES = ("I","Q","U","V");

class DiskSource(PointSource):
  def __init__(self,ns,name,direction,
               I=0.0,Q=None,U=None,V=None,
               spi=None,freq0=None,
               RM=None,
               size=None,order=0,symmetric=False):
    PointSource.__init__(self,ns,name,direction,I,Q,U,V,
                        spi,freq0,RM);
    self._symmetric = symmetric;
    #size=radius in radians
    # order=0, ring, order=1, disk
    self._order=order
    self._add_parm('sigma',size,tags='shape');

  def is_symmetric (self):
    return self._symmetric;

  def sigma (self):
    """Returns the size for this source,
    """
    return self._parm('sigma');

  def phi (self):
    """Returns the orientation node for this source""";
    return self._parm("phi");

  def transformation_matrix (self):
    pass

  def coherency (self,array=None,observation=None,nodes=None,**kw):
    coherency = nodes or self.ns.coh;
    array = Context.get_array(array);
    if not coherency(*array.ifrs()[0]).initialized():
      observation = Context.get_observation(observation);
      dir0 = observation.phase_centre;
      radec0 = dir0.radec();
      # 1/wl = freq/c
      iwl = self.ns.inv_wavelength << ((self.ns.freq<<Meq.Freq) / 2.99792458e+8);
      # 1/(wl^2): scaling factor applied to exp() argument below
      m_iwlsq = self.ns.m_inv_wavelength_sq << Meq.Sqr(iwl);
      # scaling factor of gaussian for unit flux
      gscale = self.ns.diskcomponent_scale << Meq.Constant(1.0);
      # baseline UVs
      uv_ifr = array.uv_ifr(dir0);
      # flux scale -- coherency multiplied by scale constant above
      fluxscale = self.ns.fluxscale << self.brightness(observation) * gscale;
      # transformed uv's (rotated and scaled)
      uv1sq = self.ns.uv1sq;
      u1sq = self.ns.u1sq;
      v1sq = self.ns.v1sq;
      sumuvsq = self.ns.sumuvsq;
      sigma = self._parm("sigma");
      # now generate nodes
      for ifr in array.ifrs():
        # rotate uvs and scale to wavelength
        uv1s = uv1sq(*ifr) << Meq.Sqr(uv_ifr(*ifr));
        u1s = u1sq(*ifr) << Meq.Selector(uv1s,index=0);
        v1s = v1sq(*ifr) << Meq.Selector(uv1s,index=1);
        rr  = sumuvsq(*ifr) <<Meq.Sqrt((u1s+v1s)*m_iwlsq);
        coherency(*ifr) << fluxscale * Meq.Bessel(rr*2*math.pi*sigma,order=self._order)/rr;
    return coherency;
