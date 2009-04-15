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
import Context,Parm
  
STOKES = ("I","Q","U","V");

class Shapelet(PointSource):
  def __init__(self,ns,name,direction,
               I=0.0,Q=None,U=None,V=None,
               spi=None,freq0=None,
               RM=None,
               scale=1,phi=0,solvable=False,
               modefile=None):
    PointSource.__init__(self,ns,name,direction,I,Q,U,V,
                        spi,freq0,RM);
    # read mode file
    infile=open(modefile,'r');
    all=infile.readlines()
    infile.close()
    ### remove first line because it only has RA,Dec
    all.pop(0)
    thisline=all[0].split()
    n0=int(thisline[0]);
    beta=float(thisline[1]);
    self._children=[]
    self._phi=phi
    self._scale=scale

    if not solvable:
      self._children+=[self._parm('beta_'+str(n0),beta,tags='beta')];
    else:
      self._children+=[self._parm('beta_'+str(n0),Meow.Parm(beta),tags='beta')];
    for ci in range(1,n0*n0+1):
      thisline=all[ci].split()
      modeval=float(thisline[1])
      if not solvable:
       self._children+=[self._parm('mode_'+str(ci),modeval,tags='modes')];
      else:
       self._children+=[self._parm('mode_'+str(ci),Meow.Parm(modeval),tags='modes')];

    clist=self.ns.clist;
    if not clist.initialized():
     self.ns.clist<<Meq.Composer(children=self._children)
    self._children=clist

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
    # flux scale 
    fscale=832
    fscale=1.0
    gcoh = self.ns.shapelet_coh.qadd(radec0);
    for ifr in array.ifrs():
      uv=uv_ifr(*ifr)*iwl 
      shptf= shp(*ifr)<<Meq.ShapeletVisTf(modes=self._children,phi=self._phi,scale=self._scale,dep_mask=0xff);
      gcoh(*ifr) << coherency * fscale * Meq.Compounder(children=[uv,shptf],common_axes=[hiid('L'),hiid('M')]);

    # phase shift to source position
    self.direction.make_phase_shift(nodes,gcoh,array,dir0);
