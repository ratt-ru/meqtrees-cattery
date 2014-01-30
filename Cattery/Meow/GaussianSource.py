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

import math
from Timba.TDL import *
from Timba.Meq import meq
from PointSource import *
import Context

STOKES = ("I","Q","U","V");

class GaussianSource(PointSource):
  
  
  def __init__(self,ns,name,direction,
               # standard flux parameters
               I=0.0,Q=None,U=None,V=None,
               spi=None,freq0=None,
               RM=None,
               # new-style extent paramaters: projections of the major axis FWHM
               # onto the l and m axes, plus an minor/major ratio (None for symmetric)
               lproj=None,mproj=None,ratio=None,
               # old-style (deprecated) extent parameters
               size=None,phi=0):
    PointSource.__init__(self,ns,name,direction,I,Q,U,V,
                        spi,freq0,RM);
    # check for old interface, convert to new-style arguments
    if size is not None:
      print "WARNING: using deprecated interface to Meow.GaussianSource.";
      if isinstance(size,(list,tuple)) and len(size) == 2:
        try:
          smaj,smin = map(float,size);
          pa = -float(phi);
        except:
          raise TypeError,"when using deprecated interface to Meow.Gaussian, extents must be scalar"
        ratio = smin/smaj;
        fwhm_lproj = math.sin(pa)*smaj;
        fwhm_mproj = math.cos(pa)*smaj;
      else:
        try:
          smaj = float(size);
        except:
          raise TypeError,"when using deprecated interface to Meow.Gaussian, extents must be scalar"
        fwhm_lproj,fwhm_mproj,ratio = 0,smaj,1;
      # convert to normal parameters
    # create polc(s) for size
    self._add_parm('lproj',lproj,tags="shape");
    self._add_parm('mproj',mproj,tags="shape");
    self._add_parm('ratio',ratio,tags="shape");
    
  def shape_static (self):
    """If shape is defined statically, returns tuple of lproj,mproj,ratio.
    Else returns None.""";
    parms = [ self._get_constant(x) for x in 'lproj','mproj','ratio' ];
    if all([ p is not None for p in parms ]):
      return parms;
    return None;

  def shape (self):
    """Returns tuple of lproj,mproj,ratio.""";
    return [ self._parm(x) for x in 'lproj','mproj','ratio' ];

  def coherency (self,array=None,observation=None,nodes=None,**kw):
    """Returns the Gaussian coherency (at l=m=0). We'll use PSVTensor to generate this here.
    """;
    coherency = nodes or self.ns.coh;
    array = Context.get_array(array);
    uvw = Meow.Context.array.uvw_ifr();
    # build up tensors for source
    if not self.ns.BT.initialized():
      B_static = self.brightness_static();
      if B_static:
        self.ns.BT << Meq.Constant(B_static,dims=[1,2,2]);
      else:
        self.ns.BT << Meq.Composer(self.brightness(),dims=[1,2,2]);
    # lmn tensor is 0,0,1, since this tree is supposed to predict coherency at phase center
    self.ns.lmnT ** Meq.Constant([0,0,1],dims=[1,3]);
    # shape tensor
    if not self.ns.shapeT.initialized():
      shape_static = self.shape_static();
      if shape_static:
        self.ns.shapeT << Meq.Constant(shape_static,dims=[1,3]);
      else:
        self.ns.shapeT << Meq.Composer(dims=[1,3],*self.shape());
    # coherency nodes
    if not coherency(*array.ifrs()[0]).initialized():
      for p,q in array.ifrs():
        coherency(p,q) << Meq.PSVTensor(self.ns.lmnT,self.ns.BT,uvw(p,q),self.ns.shapeT);
    return coherency;
