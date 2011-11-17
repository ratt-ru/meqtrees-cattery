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

from Timba.TDL import *
from Timba.Meq import meq
from Parameterization import *
import Jones
import Context
from math import cos,sin,acos,asin,atan2,sqrt,pi

def radec_to_lmn (ra,dec,ra0,dec0):
  """Returns l,m,n corresponding to direction ra,dec w.r.t. direction ra0,dec0""";
## our old formula, perhaps unjustly suspected by me
## See purrlog for 3C147_spw0, entries of Nov 21.
## Doesn't this break down at the pole (l always 0)?
  l = cos(dec) * sin(ra-ra0);
  m = sin(dec) * cos(dec0) - cos(dec) * sin(dec0) * cos(ra-ra0);
## Sarod's formula from LSM.common_utils. Doesn't seem to work right!
## (that's because it's for NCP lm coordinates used in NEWSTAR sky models)
#  l = sin(ra-ra0)*math.cos(dec);
#   sind0 = sin(dec0);
#   if sind0 != 0:
#     m = -(cos(ra-ra0)*cos(dec)-cos(dec0))/math.sin(dec0);
#   else:
#     m = 0
  n = sqrt(1-l*l-m*m);
  return l,m,n;

def lm_to_radec (l,m,ra0,dec0):
  """Returns ra,dec corresponding to l,m w.r.t. direction ra0,dec0""";
  # see formula at http://en.wikipedia.org/wiki/Orthographic_projection_(cartography)
  rho = sqrt(l**2+m**2);
  if rho == 0.0:
    ra = ra0
    dec = dec0
  else:
    cc = asin(rho);
    ra = ra0 + atan2( l*sin(cc),rho*cos(dec0)*cos(cc)-m*sin(dec0)*sin(cc) );
    dec = asin( cos(cc)*sin(dec0) + m*sin(cc)*cos(dec0)/rho );

  return ra,dec;

class Direction (Parameterization):
  """A Direction represents an absolute direction on the sky, in ra,dec (radians).
  'name' may be None, this usually identifies the phase centre.
  A direction may be static, in which case it is known and fixed at compile-time.
  """;
  def __init__(self,ns,name,ra,dec,static=False,
               quals=[],kwquals={}):
    Parameterization.__init__(self,ns,name,
                              quals=quals,kwquals=kwquals);
    self._add_parm('ra',ra,tags="direction solvable");
    self._add_parm('dec',dec,tags="direction solvable");
    if static and self._is_constant('ra') and self._is_constant('dec'):
      self.static = ra,dec;
      self.static_lmn = {};
    else:
      self.static = None;

  def radec (self):
    """Returns ra-dec 2-vector node for this direction.
    """;
    radec = self.ns.radec;
    if not radec.initialized():
      ra = self._parm('ra');
      dec = self._parm('dec');
      radec << Meq.Composer(ra,dec);
    return radec;

  def ra (self):
    """Returns R.A. node for this direction.""";
    return self.ns.ra ** Meq.Selector(self.radec(),index=0);

  def dec (self):
    """Returns dec node for this direction.""";
    return self.ns.dec ** Meq.Selector(self.radec(),index=1);

  def radec_static (self):
    """Returns ra-dec tuple for this direction if static, None if not.
    """;
    return self.static;

  def pa (self,xyz=None):
    """Returns parallactic angle relative to given xyz position
    """;
    xyz = xyz or Context.array.xyz0();
    pa = self.ns.pa(*xyz.quals,**xyz.kwquals);
    if not pa.initialized():
      pa << Meq.ParAngle(radec=self.radec(),xyz=xyz);
    return pa;

  def pa_rot (self,xyz=None):
    """Returns rotation matrix for PA relative to given xyz position
    """;
    pa = self.pa(xyz);
    pa_rot = pa('rot');
    if not pa_rot.initialized():
      rsin = pa('sin') ** Meq.Sin(pa);
      rcos = pa('cos') ** Meq.Cos(pa);
      pa_rot << Meq.Matrix22(rcos,-rsin,rsin,rcos);
    return pa_rot;

  def pa_invrot (self,xyz=None):
    """Returns inverse rotation matrix for PA relative to given xyz position
    """;
    pa_rot = self.pa_rot(xyz);
    pa_invrot = pa_rot('T') ** Meq.Transpose(pa_rot);
    return pa_invrot;

  def azel (self,xyz=None):
    """Returns az-el 2-vector node for this direction.
    """;
    xyz = xyz or Context.array.xyz0();
    azel = self.ns.azel(*xyz.quals,**xyz.kwquals);
    if not azel.initialized():
      azel << Meq.AzEl(self.radec(),xyz);
    return azel;

  def az (self,xyz=None):
    """Returns az node for this direction.""";
    azel = self.azel(xyz);
    return self.ns.az(*azel.quals,**azel.kwquals) ** Meq.Selector(azel,index=0);

  def el (self,xyz=None):
    """Returns el node for this direction.""";
    azel = self.azel(xyz);
    return self.ns.el(*azel.quals,**azel.kwquals) ** Meq.Selector(azel,index=1);

  def lmn (self,dir0=None):
    """Returns LMN 3-vector node for this direction, given a reference
    direction dir0, or using the global phase center if not supplied.
    Qualifiers from dir0 are added in.
    All other lmn-related methods below call on this one to get
    the basic lmn 3-vector.
    """;
    dir0 = Context.get_dir0(dir0);
    lmn = self.ns.lmn(*dir0._quals,**dir0._kwquals);
    if not lmn.initialized():
      if self is dir0:
        lmn << Meq.Constant(value=Timba.array.array((0.,0.,1.)));
      else:
        lmnst = self.lmn_static(dir0);
        if lmnst:
          lmn << Meq.Constant(value=Timba.array.array(lmnst));
        else:
          lmn << Meq.LMN(radec_0=dir0.radec(),radec=self.radec());
    return lmn;

  def lmn_static (self,dir0=None):
    """Returns static LMN tuple, given a reference direction dir0, or using the global phase
    center if not supplied.
    Both this direction and the reference direction must be static, otherwise None is returned.
    """;
    dir0 = Context.get_dir0(dir0);
    if not self.static or not dir0.static:
      return None;
    ra0,dec0 = dir0.radec_static();
    # see if we have already computed this lmn
    lmn = self.static_lmn.get((ra0,dec0),None);
    if lmn is None:
      ra,dec = self.radec_static();
      lmn = self.static_lmn[(ra0,dec0)] = radec_to_lmn(ra,dec,ra0,dec0);
    return lmn;
    
  def is_phase_centre (self):
    dir0 = Context.observation.phase_centre;
    if self is dir0:
      return True;
    if not self.static or not dir0.static:
      return False;
    return dir0.radec_static() == self.radec_static();
    
  def _lmn_component (self,name,dir0,index):
    """Helper method for below, returns part of the LMN vector.""";
    lmn = self.lmn(dir0);
    comp = self.ns[name].qadd(lmn);  # use ns0: all qualifiers are in lmn already
    # if we used self.ns, we'd have duplicate qualifiers
    if not comp.initialized():
      comp << Meq.Selector(lmn,index=index,multi=True);
    return comp;

  def lm (self,dir0=None):
    """Returns an LM 2-vector node for this direction. All args are as
    per lmn().""";
    dir0 = Context.get_dir0(dir0);
    lmnst = self.lmn_static(dir0);
    if lmnst:
      lm = self.ns.lm(*dir0._quals,**dir0._kwquals);
      return lm ** Meq.Constant(value=Timba.array.array(lmnst[0:2]));
    else:
      return self._lmn_component('lm',dir0,[0,1]);
  def l (self,dir0=None):
    """Returns an L node for this direction. All args are as per lmn().""";
    return self._lmn_component('l',dir0,0);
  def m (self,dir0=None):
    """Returns an M node for this direction. All args are as per lmn().""";
    return self._lmn_component('m',dir0,1);
  def n (self,dir0=None):
    """Returns an N node for this direction. All args are as per lmn().""";
    return self._lmn_component('n',dir0,2);

  def lmn_1 (self,dir0=None):
    """Returns L,M,N-1 3-vector node for this direction.
     All args are as per lmn().""";
    dir0 = Context.get_dir0(dir0);
    lmn_1 = self.ns.lmn_minus1(*dir0._quals,**dir0._kwquals);
    if not lmn_1.initialized():
      lmnst = self.lmn_static(dir0);
      if lmnst:
        lmn_1 << Meq.Constant(value=Timba.array.array((lmnst[0],lmnst[1],lmnst[2]-1)));
      else:
        lmn = self.lmn(dir0);
        lmn_1 << Meq.Paster(lmn,self.n(dir0)-1,index=2);
    return lmn_1;

  def KJones (self,array=None,dir0=None):
    """makes and returns a series of Kjones (phase shift) nodes
    for this direction, given a reference direction dir0, or using
    the global phase center if not supplied.
    Return value is an under-qualified node K, which should be
    qualified with a station index.
    If 'smearing' is True, then uses an alternative implementation of K --
    makes a separate node K('arg'), which holds the complex argument of K. This argument node
    can be used to compute smearing factors.
    """;
    # if direction is same, K is identity for all stations
    if self is dir0:
      Kj = self.ns.K << 1;
      return lambda p: Kj;
    else:
      Kj = self.ns.K;
      if dir0:
        Kj = Kj.qadd(dir0.radec());
      array = Context.get_array(array);
      stations = array.stations();
      Kjarg = Kj('arg');
      if not Kj(stations[0]).initialized():
        if Kjarg(stations[0]).initialized():
          for p in stations:
            Kj(p) << Meq.VisPhaseShift(Kjarg(p));
        else:
          lmn_1 = self.lmn_1(dir0);
          uvw = array.uvw(dir0);
          for p in stations:
            Kj(p) << Meq.VisPhaseShift(lmn=lmn_1,uvw=uvw(p));
      return Kj;

  def _KJonesArg (self,Kj,array,dir0):
    # if direction is same, K is identity for all stations
    Kjarg = Kj('arg');
    if self is dir0:
      Kjarg << 0;
      return lambda p: Kjarg;
    else:
      stations = array.stations();
      if not Kjarg(stations[0]).initialized():
        lmn_1 = self.lmn_1(dir0);
        uvw = array.uvw(dir0);
        for p in stations:
          Kjarg(p) << Meq.MatrixMultiply(uvw(p),lmn_1);
      return Kjarg;

  def smear_factor (self,array=None,dir0=None):
    """Returns smearing factor associated with this direction.
    Returns something that must be qualified with p,q, or can be None if there's no smearing.
    By default, uses the Direction implementation."""
    if self is dir0:
      return None;
    else:
      Kj = self.ns.K;
      if dir0:
        Kj = Kj.qadd(dir0.radec());
      Kjsm = Kj('smear');
      array = Context.get_array(array);
      ifrs = array.ifrs();
      if not Kjsm(*ifrs[0]).initialized():
        Kjarg = self._KJonesArg(Kj,array,dir0);
        for p,q in ifrs:
          Kjsm(p,q) << Meq.TFSmearFactor(Kjarg(p),Kjarg(q));
      return Kjsm;

  def make_phase_shift (self,vis,vis0,array=None,dir0=None):
    """phase-shifts visibilities given by vis0(p,q) from dir0
    (the global phase center by default) to the given direction.
    Shifted visibilities are created as vis(p,q).
    """;
    dir0 = Context.get_dir0(dir0);
    array = Context.get_array(array);
    # if direction is the same, use an Identity transform
    if self is dir0:
      for p,q in array.ifrs():
        vis(p,q) << Meq.Identity(vis0(p,q));
    # else apply KJones and smear, if need be
    else:
      Kj = self.KJones(array=array,dir0=dir0);
      Jones.apply_corruption(vis,vis0,Kj,array.ifrs());
