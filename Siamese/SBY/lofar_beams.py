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
import math
import os.path
import Meow
from Meow import Context

dirname = os.path.dirname(__file__);

TDLCompileOption('beam_library_path',
      "Beam library path",TDLDirSelect(default=dirname));
TDLCompileOption('station_config_path',
      "Station configuration file",TDLFileSelect(default=os.path.join(dirname,'AntennaCoords')));
TDLCompileOption('beam_stations',"Use station model for",[None,"all",1],more=str,
    doc="""<P>If 'None' is selected, a dipole beam model is used for all stations.</P>
<P>If 'all' is selected, a station beam model is used for all stations.</P>
<P>Otherwise enter a station ID to use a station beam for that station, and a dipole beam for the rest.</P>
""");

TDLCompileMenu("Dipole parameters",
  TDLOption('beam_L',"L",[0.9758],more=float),
  TDLOption('beam_phi0',"phi0",[0],more=float),
  TDLOption('beam_h',"h",[1.706],more=float),
  TDLOption('beam_alpha',"alpha",[math.pi/4.001],more=float),
  TDLOption('beam_az0x',"X dipole orientation",[math.pi/4],more=float),
  TDLOption('beam_az0y',"Y dipole orientation",[3*math.pi/4],more=float),
  TDLOption('beam_scale',"Rescaling factor",[88],more=float)
);
TDLCompileOption('station_phi0',"Station orientation",[math.pi/4],more=float);

def make_dipole_parameters (ns):
  ns.L << beam_L + 0j;
  ns.phi('x') << beam_phi0 + 0j;
  ns.phi('y') << beam_phi0 - math.pi/2 + 0j
  ns.h << beam_h + 0j;
  ns.alpha << beam_alpha + 0j;
  ns.az0('x') << beam_az0x;
  ns.az0('y') << beam_az0y - math.pi/2;
  ns.beam_scale << beam_scale;

def make_dipole_beam (ns,E,azX,azY,el):
  E('unscaled') << Meq.Matrix22(
    E('xx') <<  Meq.PrivateFunction(children=(ns.h,ns.L,ns.alpha,ns.phi('x'),azX,el),
                lib_name=beam_library_path+"/beam_dr_theta.so",function_name="test"),
    E('xy') <<  Meq.PrivateFunction(children=(ns.h,ns.L,ns.alpha,ns.phi('x'),azX,el),
                lib_name=beam_library_path+"/beam_dr_phi.so",function_name="test"),
    E('yx') <<  Meq.PrivateFunction(children=(ns.h,ns.L,ns.alpha,ns.phi('y'),azY,el),
                lib_name=beam_library_path+"/beam_dr_theta.so",function_name="test"),
    E('yy') <<  Meq.PrivateFunction(children=(ns.h,ns.L,ns.alpha,ns.phi('y'),azY,el),
                lib_name=beam_library_path+"/beam_dr_phi.so",function_name="test")
  );
  E << E('unscaled')/ns.beam_scale;

def make_station_beam (ns,S,azel0,azel,ref_freq):
  S << Meq.StationBeam(filename=station_config_path,
              azel_0=azel0,azel=azel,
              phi_0=station_phi0,
              ref_freq=ref_freq);


def compute_jones (Jones,sources,stations=None,pointing_offsets=None,**kw):
  ns = Jones.Subscope();
  stations = stations or Context.array.stations();
  obs = Context.get_observation(None);
  freq0 = obs.freq0();
  freq1 = obs.freq1();

  make_dipole_parameters(ns);
  # get reference frequency
  ns.ref_freq << (freq1 + freq0)/2;

  # create per-direction, per-station E Jones matrices
  for src in sources:
    for station in stations:
      Ej = Jones(src,station);
      # az,el of source relative to this station
      az = src.direction.az(Context.array.xyz(station));
      el = src.direction.el(Context.array.xyz(station));
      # azimuth relative to X/Y dipoles
      azX = Ej('azX') << az - ns.az0('x');
      azY = Ej('azY') << az - ns.az0('y');
      # make either a station beam, or a dipole beam
      if beam_stations == "all" or str(station) == str(beam_stations):
        B = Ej('B');
        S = Ej('S');
        # dipole beam
        make_dipole_beam(ns,B,azX,azY,el);
        # station gains: scalar matrix
        # if independent X/Y beams are formed, then this will become a disgonal
        # matrix (and B*S will need to be a Meq.MatrixMultiply)
        S << Meq.StationBeam(filename=station_config_path,
                              azel_0=obs.phase_center.azel(Context.array.xyz(station)),
                              azel=src.direction.azel(Context.array.xyz(station)),
                              phi_0=station_phi0,
                              ref_freq=ns.ref_freq);
        Ej << B*S;
      else:
        # dipole beam
        make_dipole_beam(ns,Ej,azX,azY,el);

  return Jones;