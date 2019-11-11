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

from Timba.TDL import *
import math
import os.path
import Meow
from Meow import Context

dirname = os.path.dirname(__file__);
libname = os.path.join(dirname,"lofar_beams_lib.so");
DEG = math.pi/180;

DIPOLES = "dipoles";
STATIONS = "stations";
MIX = "station vs. dipoles";

array_opt = TDLCompileOption('array_composition',"Array composition",[DIPOLES,STATIONS,MIX],
    doc="""<P>Three options are available: an array of dipoles, an array of proper stations,
or one proper station with the rest being dipoles.</P>""");
station_id_opt = TDLCompileOption('station_id',"Antenna ID of station",[1],more=str,
    doc="""<P>ID of antenna corresponding to the one proper station.</P>""");

TDLCompileMenu("Dipole configuration",
  TDLMenu("LOFAR LBA (droopy dipoles)",
      TDLOption('lba_L',"Dipole length (L, meters)",[0.9758],more=float),
      TDLOption('lba_h',"Dipole height (h, meters)",[1.706],more=float),
      TDLOption('lba_alpha',"Dipole slant (alpha, deg)",[45.00001],more=float),
      TDLOption('lba_phi0_x',"X dipole orientation, deg",[45],more=float),
      TDLOption('lba_phi0_y',"Y dipole orientation, deg",[135],more=float),
      TDLOption('lba_beam_scale',"Beam rescaling factor",[88],more=float,
        doc="""This is just a scale factor applied to the complex beam gain, to keep the
        numbers reasonable"""),
    toggle="lba_model"),
  TDLMenu("LOFAR HBA (bowtie)",
      TDLOption('hba_phi0_x',"X dipole orientation, deg",[45],more=float),
      TDLOption('hba_phi0_y',"Y dipole orientation, deg",[135],more=float),
      TDLOption('hba_beam_scale',"Beam rescaling factor",[600],more=float,
        doc="""This is just a scale factor applied to the complex beam gain, to keep the
        numbers reasonable"""),
    toggle="hba_model"),
  exclusive='dipole_model',
);

station_menu = TDLCompileMenu("Station configuration",
  TDLCompileOption('station_config_path',
    "Station configuration files",['%s.coords','AntennaCoords.lba','AntennaCoords.hba'],more=str,
    doc="""This file specifies the layout of dipoles within a station."""
  ),
  TDLCompileOption('station_phi0',"Station orientation, deg",[45],more=float,
    doc="""This parameter can be used to rotate the station layout."""
  )
);

def _set_array_composition (value):
  station_menu.show(value!=DIPOLES);
  station_id_opt.show(value==MIX);
array_opt.when_changed(_set_array_composition);

TDLCompileOption('beam_library_path',"Beam library",TDLFileSelect(default=libname),
  doc="""This is a dynamic library module (usually called lofar_beams_lib.so) conatining the beam computation routines.
  If you don't see this library, you might need to type "make" in the directory containing this module to build it.
  """
);

def make_hba_parameters (ns):
  ns.beam_scale << hba_beam_scale;

def lba_dipole_beam (ns,E,az,el):
  # init beam constants
  if not ns.lba_beam_scale.initialized():
    ns.L << lba_L + 0j;
    ns.phi0_x << -lba_phi0_x*DEG + 0j;
    ns.phi0_y << -lba_phi0_y*DEG + 0j
    ns.h << lba_h + 0j;
    ns.alpha << lba_alpha*DEG + 0j;
    ns.lba_beam_scale << lba_beam_scale;
  # make dipole beam
  E('unscaled') << Meq.Matrix22(
    E('xx') <<  Meq.PrivateFunction(az,el,ns.phi0_x,ns.h,ns.L,ns.alpha,
                lib_name=beam_library_path,function_name="lba_theta"),
    E('xy') <<  Meq.PrivateFunction(az,el,ns.phi0_x,ns.h,ns.L,ns.alpha,
                lib_name=beam_library_path,function_name="lba_phi"),
    E('yx') <<  Meq.PrivateFunction(az,el,ns.phi0_y,ns.h,ns.L,ns.alpha,
                lib_name=beam_library_path,function_name="lba_theta"),
    E('yy') <<  Meq.PrivateFunction(az,el,ns.phi0_y,ns.h,ns.L,ns.alpha,
                lib_name=beam_library_path,function_name="lba_phi")
  );
  E << E('unscaled')/ns.lba_beam_scale;

def hba_dipole_beam (ns,E,az,el):
  ns.hba_beam_scale ** hba_beam_scale;
  ns.phi0_x ** (-lba_phi0_x*DEG + 0j);
  ns.phi0_y ** (-lba_phi0_y*DEG + 0j);
  E('unscaled') << Meq.Matrix22(
    E('xx') <<  Meq.PrivateFunction(az,el,ns.phi0_x,
                lib_name=beam_library_path,function_name="hba_theta"),
    E('xy') <<  Meq.PrivateFunction(az,el,ns.phi0_x,
                lib_name=beam_library_path,function_name="hba_phi"),
    E('yx') <<  Meq.PrivateFunction(az,el,ns.phi0_y,
                lib_name=beam_library_path,function_name="hba_theta"),
    E('yy') <<  Meq.PrivateFunction(az,el,ns.phi0_y,
                lib_name=beam_library_path,function_name="hba_phi")
  );
  E << E('unscaled')/ns.hba_beam_scale;

def compute_jones (Jones,sources,stations=None,pointing_offsets=None,**kw):
  ns = Jones.Subscope();
  stations = stations or Context.array.stations();
  obs = Context.get_observation(None);
  freq0 = obs.freq0();
  freq1 = obs.freq1();

  if hba_model:
    dipole_model = hba_dipole_beam;
  else:
    dipole_model = lba_dipole_beam;
  # get array composition
  global station_id;
  if array_composition != MIX:
    station_id = None;
  # get reference frequency
  ns.ref_freq << (freq1 + freq0)/2;

  if "%" in station_config_path:
    make_filename = lambda p:station_config_path%p;
  else:
    make_filename = lambda p:station_config_path;

  # create per-direction, per-station E Jones matrices
  for src in sources:
    for station in stations:
      Ej = Jones(src,station);
      # az,el of source relative to this station
      az = src.direction.az(Context.array.xyz(station));
      el = src.direction.el(Context.array.xyz(station));
      # make either a station beam, or a dipole beam
      if array_composition == STATIONS or str(station) == str(station_id):
        B = Ej('B');
        S = Ej('S');
        # dipole beam
        dipole_model(ns,B,az,el);
        # station gains: scalar matrix
        # if independent X/Y beams are formed, then this will become a disgonal
        # matrix (and B*S will need to be a Meq.MatrixMultiply)
        filename = make_filename(station);
        if not os.path.exists(filename):
          raise RuntimeError("Cannot find station layout file %s"%filename);
        S << Meq.StationBeam(filename=filename,
                              azel_0=obs.phase_center.azel(Context.array.xyz(station)),
                              azel=src.direction.azel(Context.array.xyz(station)),
                              phi_0=station_phi0,
                              ref_freq=ns.ref_freq);
        Ej << B*S;
      else:
        # dipole beam
        dipole_model(ns,Ej,az,el);

  return Jones;