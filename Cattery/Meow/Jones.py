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

from Timba.TDL import *
from Timba.Meq import meq
from Meow import Parameterization
from . import Context


def gain_ap_matrix (jones,ampl=1.,phase=0.,tags=[],series=None,solvable=True):
  """Creates an amplitude-phase gain matrix of the form:
        ( ax*e^{i*px}       0      )
        (     0       ay*e^{i*p_y} )
  'jones' should be a node stub.
  'series' can be a list of qualifiers to make a series of matrices.
  The x/y amplitudes and phases are created as Meq.Parms (with initial
  values 1 and 0, respectively) by qualifying 'jones' with
  'xa','xp','ya','yp'. These Parm nodes will be tagged with the given set
  of 'tags', plus 'ampl' and 'phase'.
  """
  if Context.observation.circular():
    X,Y,XA,XP,YA,YP = 'r','l','ra','rp','la','lp';
  else:
    X,Y,XA,XP,YA,YP = 'x','y','xa','xp','ya','yp';
  # create matrix per-station, or just a single matrix
  if series:
    for p in series:
      xa = Parameterization.resolve_parameter("ampl",jones(p,XA),ampl,tags=tags,solvable=solvable);
      #print "tags here",tags;
      xp = Parameterization.resolve_parameter("phase",jones(p,XP),phase,tags=tags,solvable=solvable);
      ya = Parameterization.resolve_parameter("ampl",jones(p,YA),ampl,tags=tags,solvable=solvable);
      yp = Parameterization.resolve_parameter("phase",jones(p,YP),phase,tags=tags,solvable=solvable);
      jones(p) << Meq.Matrix22(
        jones(p,X) << Meq.Polar(xa,xp),
        0,0,
        jones(p,Y) << Meq.Polar(ya,yp)
      );
  else:
    xa = Parameterization.resolve_parameter("ampl",jones(XA),ampl,tags=tags,solvable=solvable);
    xp = Parameterization.resolve_parameter("phase",jones(XP),phase,tags=tags,solvable=solvable);
    ya = Parameterization.resolve_parameter("ampl",jones(YA),ampl,tags=tags,solvable=solvable);
    yp = Parameterization.resolve_parameter("phase",jones(YP),phase,tags=tags,solvable=solvable);
    jones << Meq.Matrix22(
      jones(X) << Meq.Polar(xa,xp),
      0,0,
      jones(Y) << Meq.Polar(ya,yp)
    );
  return jones;


def rotation_matrix (jones,rot=0.,tags=[],series=None):
  """Creates a rotation matrix.
  'jones' should be a node stub.
  'series' can be a list of qualifiers to make a series of matrices.
  The rotation angle is created as a Meq.Parm by qualifying 'jones' with
  'angle'. These Parm nodes will be tagged with the given set
  of 'tags', plus 'rotation'.
  """
  # create matrix per-station, or just a single matrix
  if series:
    for p in series:
      angle = Parameterization.resolve_parameter("rotation",jones('angle',p),rot,tags=tags);
      cosangle = jones("cosa",p) << Meq.Cos(angle);
      sinangle = jones("sina",p) << Meq.Sin(angle);
      jones(p) << Meq.Matrix22(cosangle,-sinangle,sinangle,cosangle);
  else:
    angle = Parameterization.resolve_parameter("rotation",jones('angle'),rot,tags=tags);
    cosangle = jones("cosa") << Meq.Cos(angle);
    sinangle = jones("sina") << Meq.Sin(angle);
    jones << Meq.Matrix22(cosangle,-sinangle,sinangle,cosangle);
  return jones;

def decoupled_rotation_matrix (jones,rot=0.,tags=[],series=None):
  """Creates a decoupled rotation matrix (independent angles for x and y)
  'jones' should be a node stub.
  'series' can be a list of qualifiers to make a series of matrices.
  The rotation angle is created as a Meq.Parm by qualifying 'jones' with
  'angle'. These Parm nodes will be tagged with the given set
  of 'tags', plus 'rotation'.
  """
  # create matrix per-station, or just a single matrix
  if series:
    for p in series:
      ax = Parameterization.resolve_parameter("rotation",jones('angle','x',p),rot,tags=tags);
      ay = Parameterization.resolve_parameter("rotation",jones('angle','y',p),rot,tags=tags);
      jones(p) << Meq.Matrix22(Meq.Cos(ax),-Meq.Sin(ax),Meq.Sin(ay),Meq.Cos(ay));
  else:
    ax = Parameterization.resolve_parameter("rotation",jones('angle','x'),rot,tags=tags);
    ay = Parameterization.resolve_parameter("rotation",jones('angle','y'),rot,tags=tags);
    jones << Meq.Matrix22(Meq.Cos(ax),-Meq.Sin(ax),Meq.Sin(ay),Meq.Cos(ay));
  return jones;


def ellipticity_matrix (jones,ell=0.,tags=[],series=None):
  """Creates a dipole ellipticity matrix.
  'jones' should be a node stub.
  'series' can be a list of qualifiers to make a series of matrices.
  The ellipticity angle is created as a Meq.Parm by qualifying 'jones' with
  'ell'. These Parm nodes will be tagged with the given set
  of 'tags', plus 'ellipticity'.
  """
  # create matrix per-station, or just a single matrix
  if series:
    for p in series:
      angle = Parameterization.resolve_parameter("ellipticity",jones('angle',p),ell,tags=tags);
      cosangle = jones("cosa",p) << Meq.Cos(angle);
      isinangle = jones("sina",p) << Meq.ToComplex(0,Meq.Sin(angle));
      jones(p) << Meq.Matrix22(cosangle,isinangle,isinangle,cosangle);
  else:
    angle = Parameterization.resolve_parameter("ellipticity",jones('angle'),ell,tags=tags);
    cosangle = jones("cosa") << Meq.Cos(angle);
    isinangle = jones("sina",p) << Meq.ToComplex(0,Meq.Sin(angle));
    jones << Meq.Matrix22(cosangle,isinangle,isinangle,cosangle);
  return jones;

def decoupled_ellipticity_matrix (jones,ell=0.,tags=[],series=None):
  """Creates a decoupled ellipticity matrix (independent angles for x and y).
  'jones' should be a node stub.
  'series' can be a list of qualifiers to make a series of matrices.
  The ellipticity angle is created as a Meq.Parm by qualifying 'jones' with
  'ell'. These Parm nodes will be tagged with the given set
  of 'tags', plus 'ellipticity'.
  """
  # create matrix per-station, or just a single matrix
  if series:
    for p in series:
      ax = Parameterization.resolve_parameter("ellipticity",jones('angle','x',p),ell,tags=tags);
      ay = Parameterization.resolve_parameter("ellipticity",jones('angle','y',p),ell,tags=tags);
      jones(p) << Meq.Matrix22(Meq.Cos(ax),Meq.ToComplex(0,Meq.Sin(ax)),
                               Meq.ToComplex(0,Meq.Sin(ay)),Meq.Cos(ay));
  else:
    ax = Parameterization.resolve_parameter("ellipticity",jones('angle','x'),ell,tags=tags);
    ay = Parameterization.resolve_parameter("ellipticity",jones('angle','y'),ell,tags=tags);
    jones << Meq.Matrix22(Meq.Cos(ax),Meq.ToComplex(0,Meq.Sin(ax)),
                          Meq.ToComplex(0,Meq.Sin(ay)),Meq.Cos(ay));
  return jones;


def define_rotation_matrix (angle):
  """Returns node definition for a rotation matrix, given an 'angle' node
  or value. Since this is a node _definition_, it must be bound to a node
  name with '<<', e.g.:
    ns.pa << Meq.ParAngle(...);
    ns.P << Jones.define_rotation_matrix(ns.pa);
  """
  cos = angle("cos") << Meq.Cos(angle);
  sin = angle("sin") << Meq.Sin(angle);
  return Meq.Matrix22(cos,-sin,sin,cos);


def apply_corruption (vis,vis0,jones,extra_term=None,ifrs=None):
  """Creates nodes to corrupt with a set of Jones matrices.
  'vis' is the output node which will be qualified with (p,q)
  'vis0' is an input visibility node which will be qualified with (p,q2)
  'jones' is either one unqualified Jones matrix, or else a list/tuple of
    Jones matrices. In either case each Jones term will be qualified with
    the station index (p).
  'ifrs' should be a list of p,q pairs; by default Meow.Context is used.
  """;
  if not isinstance(jones,(list,tuple)):
    jones = (jones,);
  # multiply input visibilities by our jones list
  for p,q in (ifrs or Context.array.ifrs()):
    terms = [vis0(p,q)];
    # collect list of per-source station-qualified Jones terms
    for J in jones:
      J2c = J(q)('conj') ** Meq.ConjTranspose(J(q));
      terms = [J(p)] + terms + [J2c];
    # create multiplication node
    vis(p,q) << Meq.MatrixMultiply(*terms);
  return vis;


def apply_correction (vis,vis0,jones,ifrs=None):
  """Creates nodes to apply the inverse of a set of Jones matrices.
  'vis' is the output node which will be qualified with (sta1,sta2)
  'vis0' is an input visibility node which will be qualified with (sta1,sta2)
  'jones' is either one unqualified Jones matrix, or else a list/tuple of
    Jones matrices. In either case each Jones term will be qualified with
    the station index (p).
  'ifrs' should be a list of p,q pairs; by default Meow.Context is used.
  """;
  if not isinstance(jones,(list,tuple)):
    jones = (jones,);
  # multiply input visibilities by our jones list
  for p,q in (ifrs or Context.array.ifrs()):
    terms = [vis0(p,q)];
    # collect list of per-source station-qualified Jones terms
    for J in jones:
      J1i = J(p)('inv') ** Meq.MatrixInvert22(J(p));
      J2i = J(q)('inv') ** Meq.MatrixInvert22(J(q));
      J2ci = J2i('conj') ** Meq.ConjTranspose(J2i);
      terms = [J1i] + terms + [J2ci];
    # create multiplication node
    vis(p,q) << Meq.MatrixMultiply(*terms);
  return vis;
