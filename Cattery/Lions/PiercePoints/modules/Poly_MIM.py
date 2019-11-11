from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from Timba.TDL import *
from Timba.Meq import meq
from Lions.PiercePoints import PiercePoints
import Meow


def compile_options():
    """TODO: include options to set the default values of the parms (See TIDMIM how). The number of Parms depend on the rank, so in total N_long*N_lat parameters. Typical values for the 01 and 10 parameter depend on the choice of use_lonlat. If False, typically use 1*e-5, otherwise 1*e-3."""
    return [
        TDLCompileOption("N_long", "Rank of Polynomial in X/Long direction", [1, 2, 3], more=int,
                         doc="""TIC"""),
        TDLCompileOption("N_lat", "Rank of Polynomial in Y/Lattitude direction", [1, 2, 3], more=int,
                         doc="""TEC"""),
        TDLCompileOption("height", "Altitude", [300.], more=float),
        TDLCompileOption("use_lonlat", "Use Longitude/Lattitude of PP instead of projected (x,y)", False)]


class MIM(PiercePoints.PiercePoints):

    """A Poly_MIM is a PiercePoints that creates a polynomial as functrion of the pierc points"""

    def __init__(self, ns, name, sources, stations=None, ref_station=None, tags="iono", make_log=False):
        PiercePoints.PiercePoints.__init__(self, ns, name, sources, stations, height, make_log)
        self.ref_station = ref_station

        for nlo in range(N_long):
            for nla in range(N_lat):
                name = "N:" + str(nlo) + ":" + str(nla)
                self._add_parm(name=name, value=Meow.Parm(0.), tags=tags)

    def make_tec(self):
        self.make_longlat_pp(ref_station=self.ref_station)
        ns = self.ns
        if use_lonlat:
            lon = ns['pp']('lon')
            lat = ns['pp']('lat')
        else:
            lon = ns['pp']('x')
            lat = ns['pp']('y')

        N = []
        for nlo in range(N_long):
            N.append([])
            for nla in range(N_lat):
                name = "N:" + str(nlo) + ":" + str(nla)

                N[nlo].append(self._parm(name))

        for station in self.stations:
            for src in self.src:
                tec = ns['tec'](src, station)
                sec = ns['sec'](src, station)
                if not tec.initialized():
                    Mim = None
                    for nlo in range(N_long - 1, -1, -1):
                        Miml = N[nlo][N_lat - 1]
                        for nla in range(N_lat - 2, -1, -1):
                            Miml = Miml * lat(src, station)
                            Miml = Miml + N[nlo][nla]
                        if Mim:
                            Mim = Mim + Miml
                        else:
                            Mim = Miml
                        if nlo > 0:
                            Mim = Mim * lon(src, station)
                    tec << Mim * sec

        return ns['tec']
