#
#% $Id: GaussianSource.py 5418 2007-07-19 16:49:13Z oms $ 
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

class Shapelet(PointSource):
  def __init__(self,ns,name,direction,
               I=0.0,Q=None,U=None,V=None,
               spi=None,freq0=None,
               RM=None,
               size=None,phi=0,symmetric=False,
               modefile=None):
    PointSource.__init__(self,ns,name,direction,I,Q,U,V,
                        spi,freq0,RM);
    # read mode file
    infile=open(modefile,'r');
    all=infile.readlines()
    infile.close()
    thisline=all[0].split()
    n0=int(thisline[0]);
    beta=float(thisline[1]);
    self._children=[]
    #self._children+=[Meq.ToComplex(self._parm('beta_'+str(n0),beta,tags='beta'),0)];
    self._children+=[self._parm('beta_'+str(n0),beta,tags='beta')];
    for ci in range(1,n0*n0+1):
      thisline=all[ci].split()
      modeval=float(thisline[1])
      #self._children+=[Meq.ToComplex(self._parm('mode_'+str(ci),modeval,tags='modes'),0)];
      self._children+=[self._parm('mode_'+str(ci),modeval,tags='modes')];

    clist=self.ns.clist;
    if not clist.initialized():
     self.ns.clist<<Meq.Composer(children=self._children)
    self._children=clist
    # create private function
    #self.shp=self.ns.shp
    #if not self.shp.initialized():
    #  self.shp<<Meq.PrivateFunction(children=self._children, lib_name="/home/sarod/shapelet/src/lib/shapeletpriv.so",function_name="test", dep_mask=0xff);


    # create polc(s) for size
    self._symmetric = symmetric;
    if symmetric:
      self._add_parm('sigma',size,tags='shape');
    else:
      # setup orientation
      # note: the orientation angle, phi, of the major axis
      # is defined in the direction East through South; i.e.
      # an angle of zero defines a Gaussian oriented east-west
      self._add_parm('phi',phi,tags='shape');
      if isinstance(size,(int,float)):
        s1 = s2 = size;
      elif isinstance(size,(tuple,list)):
        if len(size) != 2:
          raise TypeError,"size: one or two numeric values expected";
        s1,s2 = size;
      else:
        raise TypeError,"size: one or two numeric values expected";
      self._add_parm('sigma1',s1,tags="shape");
      self._add_parm('sigma2',s2,tags="shape");
    
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
      phi = self.phi();
      # get direction sin/cos
      cos_phi = self.ns.cos_phi << Meq.Cos(phi);
      sin_phi = self.ns.sin_phi << Meq.Sin(phi);
      # get sigma parameters
      (a,b) = (self._parm("sigma1"),self._parm("sigma2"));
      xfm << Meq.Matrix22(
          a*cos_phi,Meq.Negate(a*sin_phi),
          b*sin_phi,b*cos_phi);
    return xfm;

  def make_visibilities (self,nodes,array,observation):
    array = Context.get_array(array);
    observation = Context.get_observation(observation);
    dir0 = observation.phase_centre;
    radec0 = dir0.radec();
    # baseline UVs
    uv1 = self.ns.uv1.qadd(radec0);
    u1=self.ns.u1.qadd(radec0);
    v1=self.ns.v1.qadd(radec0);

    shp=self.ns.shp.qadd(radec0)
    #if not self.shp.initialized():


    uv_ifr = array.uv_ifr(dir0);
    # inverse wavelength 1/wl = freq/c
    iwl = self.ns0.inv_wavelength << ((self.ns0.freq<<Meq.Freq) / 2.99792458e+8);

    coherency = self.coherency(observation);
    # gaussian coherencies go here
    # scale 
    fscale=832
    gcoh = self.ns.shapelet_coh.qadd(radec0);
    for ifr in array.ifrs():
      uv=uv_ifr(*ifr)*iwl 
      #u1s = u1(*ifr) << Meq.Selector(uv_ifr(*ifr),index=0);
      #v1s = v1(*ifr) << Meq.Selector(uv_ifr(*ifr),index=1);
      #uvs1 = uv1(*ifr) << Meq.Composer(u1s,v1s)*iwl;

      #shptf= shp(*ifr)<<Meq.PrivateFunction(children=self._children, lib_name="/home/sarod/shapelet/src/lib/shapeletpriv.so",function_name="test", dep_mask=0xff);
      shptf= shp(*ifr)<<Meq.ShapeletVisTf(modes=self._children,dep_mask=0xff);
      gcoh(*ifr) << coherency * fscale * Meq.Compounder(children=[uv,shptf],common_axes=[hiid('L'),hiid('M')]);

    # phase shift to source position
    self.direction.make_phase_shift(nodes,gcoh,array,dir0);
