from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from Timba.TDL import *
from Timba.Meq import meq
from Meow.Parameterization import *
import Meow
from Meow import SkyComponent
from Meow import ParmGroup
from Lions.xyzComponent import xyzComponent


class MIM_model(Parameterization):

    """A MIM_model creates TEC values/phase shifts given a source and a station.
    The source can either be a Meow.Skycomponent or a xyzComponent.
    """

    def __init__(self, ns, name, sources, array=None):
        Parameterization.__init__(self, ns, name)

        if isinstance(sources, Meow.SkyComponent) or isinstance(sources, xyzComponent):
            self.src = [sources]

        else:
            if isinstance(sources, list):
                self.src = sources

            else:
                raise TypeError("sources must be a (list of) Sky/xyz-Component")
        self.array = array
        if not self.array:
            # get from Context:
            self.array = Meow.Context.get_array(False)

        self.stations = self.array.stations()

    def make_tec(self):
        """abstract method, creates nodes to compute the TEC-values"""
        raise TypeError(type(self).__name__ + ".make_tec() not defined")

    def make_phase_error(self):
        """create phase errors, default = exp(i lambda * 25 * tec)"""
        tec = self.make_tec()
        freq = self.make_freq()
        for src in self.src:
            for station in self.stations:
                if not self.ns['mim_phase'](src, station).initialized():
                    self.ns['mim_phase'](src, station) << (75e8 / freq) * tec(src, station)

        return self.ns['mim_phase']

    def make_freq(self):
        if not self.ns['freq'].initialized():
            self.ns['freq'] << Meq.Freq()
        return self.ns['freq']

    def get_az(self):
        self.make_azel()
        return self.ns['az']

    def get_el(self):
        self.make_azel()
        return self.ns['el']

    def get_azel(self):
        return self.make_azel()

    def make_azel(self):
        # create azel nodes for all stat src combinations
        xyz = self.array.xyz()
        for station in self.array.stations():
            for src in self.src:
                if not self.ns['azel'](src, station).initialized():
                    # if meow.skycomponent get from there
                    if isinstance(src, Meow.SkyComponent):
                        azel = self.ns['azel'](src, station) << Meq.AzEl(src.radec(), xyz(station))

                    self.ns['el'](src, station) << Meq.Selector(azel, index=1)
                    self.ns['az'](src, station) << Meq.Selector(azel, index=0)

        return self.ns['azel']

    def compute_jones(self, Jones, **kw):
        phase = self.make_phase_error()
        for src in self.src:
            for station in self.stations:
                if not Jones(src, station).initialized():
                    Jones(src, station) << Meq.Matrix22(Meq.Polar(1, -phase(src, station)), 0, 0, Meq.Polar(1, -phase(src, station)))

        return Jones

# def make_solve_jobs(self,nodes,tags =["iono"],label=''):
# make parmgroups for solvables

# tags = tags + " solvable"
# print "Search",tags,nodes.search(tags=tags)
# self.iono_pg = ParmGroup.ParmGroup(label+"_iono",
# nodes.search(tags=tags),
# table_name="%s_iono"%label,bookmark=4);
# make solvejobs
# print "SolveJobn???"
# ParmGroup.SolveJob("cal_"+label+"_iono","Calibrate %s ionosphere terms"%label,self.iono_pg);
