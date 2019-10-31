from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from Timba.TDL import *
from Timba.Meq import meq
from Lions.PiercePoints import PiercePoints
import Meow
import math
import numpy

# Include 2 wave TID, giving for each wave the wavelength, direction, speed
# and relative amplitude. The phase offset for both waves is set to 0.


def compile_options():
    return [TDLCompileOption("TEC0", "Zero offset TEC value", [1., 5., 10.], more=float,),
            TDLCompileOption("height", "Altitude", [200., 300., 400.], more=float),
            TDLCompileOption("Wavelength_1", "Wavelength first wave (km)", [250., 500., 1000.], more=float),
            TDLCompileOption("Wavelength_2", "Wavelength second wave (km)", [250., 500., 1000.], more=float),
            TDLCompileOption("Speed_1", "Speed first wave (km/h)", [300., 600., 1200.], more=float),
            TDLCompileOption("Speed_2", "Speed second wave (km/h)", [300., 600., 1200.], more=float),
            TDLCompileOption("Theta_1", "Direction first wave (degrees)", [15., 30., 45., 60., 75., 90.], more=float, doc="""Angle of
 propagation, counter clock-wise from East"""),
            TDLCompileOption("Theta_2", "Direction second wave (degrees)", [15., 30., 45., 60., 75., 90.], more=float, doc="""Angle o
f propagation, counter clock-wise from East"""),
            TDLCompileOption("Amp_1", "Relative amplitude first wave", [0.1, 0.01, 0.02, 0.05, 0.1], more=float),
            TDLCompileOption("Amp_2", "Relative amplitude second wave", [0.0, 0.01, 0.02, 0.05, 0.1], more=float),
            TDLCompileOption("use_lonlat", "Use Longitude/Lattitude of PP instead of projected (x,y)", False)]
# def compile_options():
# return [
# TDLCompileOption("TEC0","Zero offset TEC value",[1.,5.,10.],more=float,
# doc="""TEC0"""),
# TDLCompileOption("Wavelength_x","Wavelength X or Longitude (km)",[250.,500.,1000.],more=float),
# TDLCompileOption("Wavelength_y","Wavelength Y or Lattitude (km)",[250.,500.,1000.],more=float),
# TDLCompileOption("Speed_x","Phase speed X or Longitude (km/h)",[300.,600.,1200.],more=float),
# TDLCompileOption("Speed_y","Phase speed Y or Lattitude (km/h)",[300.,600.,1200.],more=float),
# TDLCompileOption("Amp_x","Relative amplitude X or Longitude",[0.0,0.01,0.02,0.05,0.1],more=float),
# TDLCompileOption("Amp_y","Relative amplitude Y or Lattitude",[0.0,0.01,0.02,0.05,0.1],more=float),
# TDLCompileOption("use_lonlat","Use Longitude/Lattitude of PP instead of projected (x,y)",False)];


class MIM(PiercePoints.PiercePoints):

    """Create MIM_model with travelling waves as function of the pierc points"""

# def __init__(self,ns,name,sources,stations=None,height=300,ref_station=1,tags="iono",make_log=False):
# PiercePoints.PiercePoints.__init__(self,ns,name,sources,stations,height,make_log);
# self.ref_station=ref_station;
# self._add_parm(name="TEC0",value=Meow.Parm(TEC0),tags=tags)
# self._add_parm(name="Amp_x",value=Meow.Parm(Amp_x),tags=tags)
# self._add_parm(name="Amp_y",value=Meow.Parm(Amp_y),tags=tags)
# self._add_parm(name="Wavelength_x",value=Meow.Parm(Wavelength_x),tags=tags)
# self._add_parm(name="Wavelength_y",value=Meow.Parm(Wavelength_y),tags=tags)
# self._add_parm(name="Speed_x",value=Meow.Parm(Speed_x),tags=tags)
# self._add_parm(name="Speed_y",value=Meow.Parm(Speed_y),tags=tags)
    def __init__(self, ns, name, sources, stations=None, ref_station=None, tags="iono", make_log=False):
        PiercePoints.PiercePoints.__init__(self, ns, name, sources, stations, height, make_log)

        if use_lonlat:
            earth_radius = 6365
            # convert km to rad and km/h to rad/s
            W1 = Wavelength_1 / (earth_radius + height)
            W2 = Wavelength_2 / (earth_radius + height)
            Sp1 = (Speed_1 / (earth_radius + height)) / 3600.
            Sp2 = (Speed_2 / (earth_radius + height)) / 3600.
        else:
            # convert km to m and km/h to m/s
            W1 = Wavelength_1 * 1000
            W2 = Wavelength_2 * 1000
            Sp1 = Speed_1 / 3.6
            Sp2 = Speed_2 / 3.6
        self.ref_station = ref_station
        self._add_parm(name="TEC0", value=Meow.Parm(TEC0), tags=tags)
        self._add_parm(name="height", value=Meow.Parm(height), tags=tags)
        self._add_parm(name="Amp_1", value=Meow.Parm(Amp_1), tags=tags)
        self._add_parm(name="Amp_2", value=Meow.Parm(Amp_2), tags=tags)
        self._add_parm(name="Wavelength_1", value=Meow.Parm(W1), tags=tags)
        self._add_parm(name="Wavelength_2", value=Meow.Parm(W2), tags=tags)
        self._add_parm(name="Speed_1", value=Meow.Parm(Sp1), tags=tags)
        self._add_parm(name="Speed_2", value=Meow.Parm(Sp2), tags=tags)
        self._add_parm(name="Theta_1", value=Meow.Parm(Theta_1), tags=tags)
        self._add_parm(name="Theta_2", value=Meow.Parm(Theta_2), tags=tags)
        self.make_display_grid(W1, W2, Theta_1, Theta_2, Amp_1, Amp_2, TEC0)

    def make_display_grid(self, W1, W2, Theta_1, Theta_2, Amp_1, Amp_2, TEC0):

        T0 = TEC0
        A1 = Amp_1
        A2 = Amp_2
        TH1 = Theta_1
        TH2 = Theta_2
        PP = numpy.meshgrid(numpy.arange(0, 10000e3, 100e3), numpy.arange(0, 10000e3, 100e3))
        data = A1 * T0 * numpy.cos((2 * math.pi / (W1)) * (numpy.cos(TH1) * PP[0])) +\
            A1 * T0 * numpy.cos((2 * math.pi / (W1)) * (numpy.sin(TH1) * PP[1])) +\
            A2 * T0 * numpy.cos((2 * math.pi / (W2)) * (numpy.cos(TH2) * PP[0])) +\
            A2 * T0 * numpy.cos((2 * math.pi / (W2)) * (numpy.sin(TH2) * PP[1]))
        numpy.save("phase_screen_display_TID", data)

    def make_time(self):
        if not self.ns['time'].initialized():
            self.ns['time'] << Meq.Time()
        return self.ns['time']

    def make_tec(self, tags=[]):
        time = self.make_time()
        ns = self.ns
        # Get the pp positions
        if use_lonlat:
            self.make_longlat_pp(ref_station=self.ref_station)
            PP_x = ns['pp']('lon')
            PP_y = ns['pp']('lat')
        else:
            self.make_xy_pp(ref_station=self.ref_station)
            PP_x = ns['pp']('x')
            PP_y = ns['pp']('y')

# T0= self._parm("TEC0");
# Ax= self._parm("Amp_x");
# Ay= self._parm("Amp_y");
# Wx= self._parm("Wavelength_x");
# Wy= self._parm("Wavelength_y");
# Vx= self._parm("Speed_x");
# Vy= self._parm("Speed_y");
        T0 = self._parm("TEC0")
        A1 = self._parm("Amp_1")
        A2 = self._parm("Amp_2")
        W1 = self._parm("Wavelength_1")
        W2 = self._parm("Wavelength_2")
        V1 = self._parm("Speed_1")
        V2 = self._parm("Speed_2")
        # convert direction angle to radians
        TH1 = self._parm("Theta_1") * math.pi / 180
        TH2 = self._parm("Theta_2") * math.pi / 180

        for station in self.stations:
            for src in self.src:
                tec = ns['tec'](src, station)
                sec = ns['sec'](src, station)
                if not tec.initialized():
# Work-around to get x,y positions for pierce-point
# tec << T0 + sec*\
# (Ax*T0*Meq.Sin((2*math.pi/(1000*Wx))*(PP_x(src,station)-Vx*time/3.6)) + \
# Ay*T0*Meq.Sin((2*math.pi/(1000*Wy))*(PP_y(src,station)-Vy*time/3.6)))
                    tec << T0 + sec * (
                        A1 * T0 * Meq.Cos((2 * math.pi / (W1)) * (Meq.Cos(TH1) * PP_x(src, station) - V1 * time)) +
                        A1 * T0 * Meq.Cos((2 * math.pi / (W1)) * (Meq.Sin(TH1) * PP_y(src, station) - V1 * time)) +
                        A2 * T0 * Meq.Cos((2 * math.pi / (W2)) * (Meq.Cos(TH2) * PP_x(src, station) - V2 * time)) +
                        A2 * T0 * Meq.Cos((2 * math.pi / (W2)) * (Meq.Sin(TH2) * PP_y(src, station) - V2 * time)))

        return ns['tec']
