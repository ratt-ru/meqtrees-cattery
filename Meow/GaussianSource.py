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

import math
from Timba.TDL import *
from Timba.Meq import meq
from PointSource import *
import Context
  
STOKES = ("I","Q","U","V");

class GaussianSource(PointSource):
  def __init__(self,ns,name,direction,
               I=0.0,Q=None,U=None,V=None,
               spi=None,freq0=None,
               RM=None,
               size=None,phi=0):
    PointSource.__init__(self,ns,name,direction,I,Q,U,V,
                        spi,freq0,RM);
    # create polc(s) for size
    if isinstance(size,(list,tuple)) and len(size) == 2:
      self._symmetric = False;
      # setup orientation
      # note: the orientation angle, phi, of the major axis
      # is defined in the direction East through South; i.e.
      # an angle of zero defines a Gaussian oriented east-west
      self._add_parm('phi',phi,tags='shape');
      self._add_parm('sigma1',size[0],tags="shape");
      self._add_parm('sigma2',size[1],tags="shape");
    else:
      self._symmetric = True;
      self._add_parm('sigma',size,tags='shape');
    
  def is_symmetric (self):
    return self._symmetric;
    
  def sigma (self):
    """Returns the size for this source (single node for symmetric,
    two-pack for elliptic).""";
    if self._symmetric:
      return self._parm('sigma');
    else:
      return self.ns.sigma ** Meq.Composer(self._parm("sigma1"),self._parm("sigma2"));
    return size;
      
  def phi (self):
    """Returns the orientation node for this source""";
    return self._parm("phi");
    
  def transformation_matrix (self):
    # for a symmetric case, the transformation matrix is just multiplication 
    # by sigma
    if self.is_symmetric():
      return self._parm("sigma");
    # else build up full rotation-scaling matrix
    xfm = self.ns.xfmatrix;
    if not xfm.initialized():
      phi,a,b = [ self._get_constant(name) for name in "phi","sigma1","sigma2" ]; 
      # if phi,a and b are all constants, make constant matrix, else make matrix nodes
      if phi is not None and a is not None and b is not None:
        xfm << Meq.Matrix22(a*math.cos(phi),-a*math.sin(phi),b*math.sin(phi),b*math.cos(phi));
      else:
        phi,a,b = [ self._parm(name) for name in "phi","sigma1","sigma2" ]; 
        # get direction sin/cos
        cos_phi = self.ns.cos_phi << Meq.Cos(phi);
        sin_phi = self.ns.sin_phi << Meq.Sin(phi);
        xfm << Meq.Matrix22(
            a*cos_phi,Meq.Negate(a*sin_phi),
            b*sin_phi,b*cos_phi);
    return xfm;

  def make_visibilities (self,nodes,array,observation,**kw):
    array = Context.get_array(array);
    observation = Context.get_observation(observation);
    dir0 = observation.phase_centre;
    radec0 = dir0.radec();
    # 1/wl = freq/c
    iwl = self.ns0.inv_wavelength << ((self.ns0.freq<<Meq.Freq) / 2.99792458e+8);
    # -1/2*(wl^2): scaling factor applied to exp() argument below, since we want coordinates in wavelengths
    m_iwlsq = self.ns0.m_inv_wavelength_sq << -2*Meq.Sqr(iwl);
    # scaling factor of gaussian for unit flux
    # gscale = self.ns0.gaussiancomponent_scale << Meq.Constant(0.5*math.pi);
    gscale = self.ns0.gaussiancomponent_scale << Meq.Constant(1);
    # baseline UVs
    uv_ifr = array.uv_ifr(dir0);
    # rotation matrix
    xfm = self.transformation_matrix();
    # flux scale -- coherency multiplied by scale constant above
    fluxscale = self.ns.fluxscale.qadd(radec0()) \
          << self.coherency(observation) * gscale;
    # transformed uv's (rotated and scaled)
    uv1sq = self.ns.uv1sq.qadd(radec0);
    u1sq = self.ns.u1sq.qadd(radec0);
    v1sq = self.ns.v1sq.qadd(radec0);
    # gaussian coherencies go here
    gcoh = self.ns.gauss_coh.qadd(radec0);
    # now generate nodes
    for ifr in array.ifrs():
      # rotate uvs and scale to wavelength
      uv1s = uv1sq(*ifr) << Meq.Sqr(Meq.MatrixMultiply(xfm,uv_ifr(*ifr)));
      u1s = u1sq(*ifr) << Meq.Selector(uv1s,index=0); 
      v1s = v1sq(*ifr) << Meq.Selector(uv1s,index=1); 
      gcoh(*ifr) << fluxscale * Meq.Exp((u1s+v1s)*m_iwlsq);
    # phase shift to source position
    self.direction.make_phase_shift(nodes,gcoh,array,dir0,smearing=self.smearing);
