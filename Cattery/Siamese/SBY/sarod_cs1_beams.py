# NB: this file has been superceded by lofar_beams.py
#

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
from Meow.Utils import *
from Meow.Direction import *
from Meow.PointSource import *
from Meow.GaussianSource import *
from Meow import Context
import Meow



#################### beams beams ...
def makebeam_droopy_phi(ns,station=0,pol='X',L=0.9758, phi0=0, alpha=math.pi/4.001, h=1.706, solvable=False,meptable=None):
    p_L=ns.p_L(pol,station)
    if not p_L.initialized():
       p_L<<Meq.ToComplex(Meq.Parm(L,node_groups='Parm', table_name=meptable,auto_save=True,tags="Ephi"),0);
    p_phi=ns.p_phi(pol,station)
    if pol=='Y':
      # add extra pi/2 
      if not p_phi.initialized():
        p_phi<<Meq.ToComplex(Meq.Parm(phi0-math.pi/2,node_groups='Parm', table_name=meptable,auto_save=True,tags="Ephi"),0);
    else:
      if not p_phi.initialized():
        p_phi<<Meq.ToComplex(Meq.Parm(phi0,node_groups='Parm', table_name=meptable,auto_save=True,tags="Ephi"),0);
    p_h=ns.p_h(pol,station)
    if not p_h.initialized():
      p_h<<Meq.ToComplex(Meq.Parm(h,node_groups='Parm', table_name=meptable,auto_save=True,tags="Ephi"),0);
    p_alpha=ns.p_alpha(pol,station)
    if not p_alpha.initialized():
      p_alpha<<Meq.ToComplex(Meq.Parm(alpha,node_groups='Parm', table_name=meptable,auto_save=True,tags="Ephi"),0);

    beam = ns.Ephi(pol,station) << Meq.PrivateFunction(children =(p_h,p_L,p_alpha,p_phi),
        lib_name=beam_library_path+"/beam_dr_phi.so",function_name="test");
    return beam;

def makebeam_droopy_theta(ns,station=0,pol='X',L=0.9758, phi0=0, alpha=math.pi/4.001, h=1.706, solvable=False,meptable=None):
    p_L=ns.p_L(pol,station)
    if not p_L.initialized():
       p_L<<Meq.ToComplex(Meq.Parm(L,node_groups='Parm', table_name=meptable,auto_save=True,tags="Etheta"),0);
    p_phi=ns.p_phi(pol,station)
    if pol=='Y':
      # add extra pi/2 
      if not p_phi.initialized():
        p_phi<<Meq.ToComplex(Meq.Parm(phi0-math.pi/2,node_groups='Parm', table_name=meptable,auto_save=True,tags="Etheta"),0);
    else:
      if not p_phi.initialized():
        p_phi<<Meq.ToComplex(Meq.Parm(phi0,node_groups='Parm', table_name=meptable,auto_save=True,tags="Etheta"),0);
    p_h=ns.p_h(pol,station)
    if not p_h.initialized():
      p_h<<Meq.ToComplex(Meq.Parm(h,node_groups='Parm', table_name=meptable,auto_save=True,tags="Etheta"),0);
    p_alpha=ns.p_alpha(pol,station)
    if not p_alpha.initialized():
      p_alpha<<Meq.ToComplex(Meq.Parm(alpha,node_groups='Parm', table_name=meptable,auto_save=True,tags="Etheta"),0);

    beam = ns.Etheta(pol,station) << Meq.PrivateFunction(children =(p_h,p_L,p_alpha,p_phi),
        lib_name=beam_library_path+"/beam_dr_theta.so",function_name="test");
    return beam;




#### even more complex droopy dipole - with station beam
def CS1_LBA_beam(Jones,sources,stations=None,pointing_offsets=None,**kw):
  Bx_phi={}
  Bx_theta={}
  By_phi={}
  By_theta={}
  
  ns=Jones.Subscope()
  stations = stations or Context.array.stations();
  radec0= Context.get_observation(None).radec0();
  freq0=Context.get_observation(None).freq0();
  freq1=Context.get_observation(None).freq1();
  array=Context.array;
  meptable=None
  solvable=None

  for station in stations:
   Bx_phi[station] = makebeam_droopy_phi(ns,station=station,meptable=meptable,solvable=solvable);
   Bx_theta[station] = makebeam_droopy_theta(ns,station=station,meptable=meptable,solvable=solvable);
   By_phi[station] = makebeam_droopy_phi(ns,pol='Y',station=station,meptable=meptable,solvable=solvable);
   By_theta[station] = makebeam_droopy_theta(ns,pol='Y',station=station,meptable=meptable,solvable=solvable);

  Ej0 = Jones;

  # get array xyz
  xyz=array.xyz();

  Xstatbeam=ns.Xstatbeam<<Meq.StationBeam(filename=beam_library_path+"/AntennaCoords",radec=radec0,xyz=array.xyz0(),phi0=Meq.Constant(-math.pi/4),ref_freq=Meq.Add(freq0,freq1)/2)
  Ystatbeam=ns.Ystatbeam<<Meq.StationBeam(filename=beam_library_path+"/AntennaCoords",radec=radec0,xyz=array.xyz0(),phi0=Meq.Constant(-math.pi/4),ref_freq=Meq.Add(freq0,freq1)/2)

  # create per-direction, per-station E Jones matrices
  for src in sources:
    dirname = src.direction.name;
    radec=src.direction.radec()
    Ej = Ej0(dirname);

    # create Az,El per source, using station 1
    azelnode=ns.azel(dirname)<<Meq.AzEl(radec=src.direction.radec(),xyz=xyz(1))
    # make shifts
    az=ns.az(dirname)<<Meq.Selector(azelnode,multi=True,index=[0])
    azX=ns.azX(dirname)<<az-math.pi/4
    azY=ns.azY(dirname)<<az-math.pi/4
    el=ns.el(dirname)<<Meq.Selector(azelnode,multi=True,index=[1])
   

    for station in stations:
        azelX =ns.azelX(dirname,station)<<Meq.Composer(azX,el)
        azelY =ns.azelY(dirname,station)<<Meq.Composer(azY,el)
        Xediag_phi = ns.Xediag_phi(dirname,station) << Meq.Compounder(children=[azelX,Bx_phi[station]],common_axes=[hiid('l'),hiid('m')])
        Xediag_theta= ns.Xediag_theta(dirname,station) << Meq.Compounder(children=[azelX,Bx_theta[station]],common_axes=[hiid('l'),hiid('m')])
        Yediag_phi = ns.Yediag_phi(dirname,station) << Meq.Compounder(children=[azelY,By_phi[station]],common_axes=[hiid('l'),hiid('m')])
        Yediag_theta = ns.Yediag_theta(dirname,station) << Meq.Compounder(children=[azelY,By_theta[station]],common_axes=[hiid('l'),hiid('m')])
        # create E matrix, normalize for zenith at 60MHz
        if beam_stations == "all" or str(station) == str(beam_stations):
          Xstatgain=ns.Xstatgain(dirname,station)<<Meq.Compounder(children=[azelX,Xstatbeam],common_axes=[hiid('l'),hiid('m')])
          Ystatgain=ns.Ystatgain(dirname,station)<<Meq.Compounder(children=[azelY,Ystatbeam],common_axes=[hiid('l'),hiid('m')])
          Ej(station) <<Meq.Matrix22(Xstatgain*Xediag_theta,Xstatgain*Xediag_phi,Ystatgain*Yediag_theta,Ystatgain*Yediag_phi)/88.00
        else:
          Ej(station) <<Meq.Matrix22(Xediag_theta,Xediag_phi,Yediag_theta,Yediag_phi)/88.00

  return Ej0;


def compute_jones (Jones,sources,stations=None,pointing_offsets=None,**kw):
  beam_model(Jones,sources,pointing_offsets=pointing_offsets,stations=stations)

  return Jones;


TDLCompileOption('beam_model',"Beam model",[CS1_LBA_beam]);
TDLCompileOption('beam_library_path',
      "Beam library path",TDLDirSelect());
TDLCompileOption('beam_stations',"Use station model for",[None,"all",1],more=str,
    doc="""<P>If 'None' is selected, a dipole beam model is used for all stations.</P>
<P>If 'all' is selected, a station beam model is used for all stations.</P>
<P>Otherwise enter a station ID to use a station beam for that station, and a dipole beam for the rest.</P>
""");
