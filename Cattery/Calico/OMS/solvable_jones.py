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

from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

"""This implements generic solvable Jones modules.
  DiagAmplPhase is a diagonal matrix, with separate solvable groups for
      amplitudes and phases
  FullRealImag is a full 2x2 matric, with separate solvable groups for
      diagonal and off-diagonal terms.
DiagAmplPhase is currently untested. In any case the solver shows poor convergence
for this type of matrix decomposition, so FullRealImag is preferred.
""";

from Timba.TDL import *
import Meow
from Meow import Context
from Meow import Jones,ParmGroup,Bookmarks
from Meow.Parameterization import resolve_parameter

class DiagAmplPhase (object):
  def __init__ (self,label):
    self.options = [];

  def runtime_options (self):
    return self.options;

  def compute_jones (self,nodes,stations=None,tags=None,label='',**kw):
    stations = stations or Context.array.stations();
    g_ampl_def = Meow.Parm(1);
    g_phase_def = Meow.Parm(0);
    nodes = Jones.gain_ap_matrix(nodes,g_ampl_def,g_phase_def,tags=tags,series=stations);

    # make parmgroups for phases and gains
    self.pg_phase = ParmGroup.ParmGroup(label+"_phase",
                    nodes.search(tags="solvable phase"),
                    table_name="%s_phase.fmep"%label,bookmark=4);
    self.pg_ampl  = ParmGroup.ParmGroup(label+"_ampl",
                    nodes.search(tags="solvable ampl"),
                    table_name="%s_ampl.fmep"%label,bookmark=4);

    # make solvejobs
    ParmGroup.SolveJob("cal_"+label+"_phase","Calibrate %s phases"%label,self.pg_phase);
    ParmGroup.SolveJob("cal_"+label+"_ampl","Calibrate %s amplitudes"%label,self.pg_ampl);

    return nodes;

class FullRealImag (object):
  """This implements a full 2x2 Jones matrix with solvable real and imaginary parts""";
  def __init__ (self,label):
    self.tdloption_namespace = label+".fullrealimag";
    self.options = [
      TDLOption("init_diag","Initial value, diagonal",[0,1],default=1,more=complex,namespace=self),
      TDLOption("init_offdiag","Initial value, off-diagonal",[0,1],default=0,more=complex,namespace=self),
    ];
    self._offdiag = True;

  def compile_options (self):
    return self.options;

  def compute_jones (self,jones,stations=None,tags=None,label='',**kw):
    stations = stations or Context.array.stations();
    # set up qualifier labels depending on polarization definition
    if Context.observation.circular():
      x,y,X,Y = 'r','l','R','L';
    else:
      x,y,X,Y = 'x','y','X','Y';
    xx,xy,yx,yy = x+x,x+y,y+x,y+y;
    rxx,rxy,ryx,ryy = [ q+":r" for q in (xx,xy,yx,yy) ];
    ixx,ixy,iyx,iyy = [ q+":i" for q in (xx,xy,yx,yy) ];
    # create parm definitions for each jones element
    tags = NodeTags(tags) + "solvable";
    diag_real = Meq.Parm(complex(self.init_diag).real,tags=tags+"diag real");
    diag_imag = Meq.Parm(complex(self.init_diag).imag,tags=tags+"diag imag");
    offdiag_real = Meq.Parm(complex(self.init_offdiag).real,tags=tags+"offdiag real");
    offdiag_imag = Meq.Parm(complex(self.init_offdiag).imag,tags=tags+"offdiag imag");
    # now loop to create nodes
    for p in stations:
      jones(p) << Meq.Matrix22(
        jones(p,xx) << Meq.ToComplex(
            jones(p,rxx) << diag_real,
            jones(p,ixx) << diag_imag
        ),
        (jones(p,xy) << Meq.ToComplex(
            jones(p,rxy) << offdiag_real,
            jones(p,ixy) << offdiag_imag
        )) if self._offdiag else 0,
        (jones(p,yx) << Meq.ToComplex(
            jones(p,ryx) << offdiag_real,
            jones(p,iyx) << offdiag_imag
        )) if self._offdiag else 0,
        jones(p,yy) << Meq.ToComplex(
            jones(p,ryy) << diag_real,
            jones(p,iyy) << diag_imag
        )
      );
    # make parmgroups for diagonal and off-diagonal terms
    subgroups = [
      ParmGroup.Subgroup(X+X,[jones(p,zz) for zz in (rxx,ixx) for p in stations]),
      ParmGroup.Subgroup(Y+Y,[jones(p,zz) for zz in (ryy,iyy) for p in stations]),
      ParmGroup.Subgroup("real part",[jones(p,zz) for zz in (rxx,ryy) for p in stations]),
      ParmGroup.Subgroup("imaginary part",[jones(p,zz) for zz in (ixx,iyy) for p in stations])
    ];
    subgroups += [
      ParmGroup.Subgroup("station %s"%p,[jones(p,zz) for zz in (rxx,ixx,ryy,iyy) ])
      for p in stations
    ];
    self.pg_diag  = ParmGroup.ParmGroup(label+"_diag",
            [ jones(p,zz) for zz in (rxx,ixx,ryy,iyy) for p in stations ],
            subgroups = subgroups,
            table_name="%s_diag.fmep"%label,bookmark=False);
    # make bookmarks
    Bookmarks.make_node_folder("%s diagonal terms"%label,
      [ jones(p,zz) for p in stations for zz in [xx,yy] ],sorted=True);
    # make solvejobs
    ParmGroup.SolveJob("cal_"+label+"_diag","Calibrate %s diagonal terms"%label,self.pg_diag);
    
    if self._offdiag:
      subgroups = [
        ParmGroup.Subgroup(X+Y,[jones(p,zz) for zz in (rxy,ixy) for p in stations]),
        ParmGroup.Subgroup(Y+X,[jones(p,zz) for zz in (ryx,iyx) for p in stations]),
        ParmGroup.Subgroup("real part",[jones(p,zz) for zz in (rxy,ryx) for p in stations]),
        ParmGroup.Subgroup("imaginary part",[jones(p,zz) for zz in (ixy,iyx) for p in stations])
      ];
      subgroups += [
        ParmGroup.Subgroup("station %s"%p,[jones(p,zz) for zz in (rxy,ixy,ryx,iyx) ])
        for p in stations
      ];
      self.pg_offdiag  = ParmGroup.ParmGroup(label+"_offdiag",
              [ jones(p,zz) for zz in (rxy,ixy,ryx,iyx) for p in stations ],
              subgroups = subgroups,
              table_name="%s_offdiag.fmep"%label,bookmark=False);
      # make bookmarks
      Bookmarks.make_node_folder("%s off-diagonal terms"%label,
        [ jones(p,zz) for p in stations for zz in [xy,yx] ],sorted=True);
      # make solvejobs
      ParmGroup.SolveJob("cal_"+label+"_offdiag","Calibrate %s off-diagonal terms"%label,self.pg_offdiag);

    return jones;

class DiagRealImag (FullRealImag):
  """This implements a diagonal 2x2 Jones matrix with solvable real and imaginary parts""";
  # note that this is not really proper OO design, but it's convenient
  def __init__ (self,label):
    self.tdloption_namespace = label+".diagrealimag";
    self.options = [
      TDLOption("init_diag","Initial value",[0,1],default=1,more=complex,namespace=self),
    ];
    self._offdiag = False;
    self.init_offdiag = 0;

