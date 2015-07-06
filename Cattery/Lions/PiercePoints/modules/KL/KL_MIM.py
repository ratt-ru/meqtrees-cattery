from Timba.TDL import *
from Timba.Meq import meq
from Lions.PiercePoints.PiercePoints import *
from Lions.PiercePoints.modules.KL import KLNode
import Meow


def compile_options():
    """KL Model is only used for solving!!! Use Kolmogorov_MIM for equivalent Simulation"""
    return [TDLCompileOption("rank", "Rank", [1, 5, 10], more=int, doc="""Cut off"""),
            TDLCompileOption("height", "Altitude", [200., 300., 400.], more=float, doc="""Altitude of the ionospheric layer""")]


class MIM(PiercePoints):

    """Create MIM_model with KL transform"""

    def __init__(self, ns, name, sources, stations=None, ref_station=None, tags="iono", make_log=False):
        PiercePoints.__init__(self, ns, name, sources, stations, height, make_log)
        self.ref_station = ref_station
        for i in range(rank):
            name = "KLParm:" + str(i)
            self._add_parm(name=name, value=Meow.Parm(0, tiling=record(time=1)), tags=tags)

# def make_phase_error(self):
# """create phase errors"""
# pp=self.make_pp(ref_station=self.ref_station);
# KL_node=self.create_KL_node();
# n=0;
# for src in self.src:
# for station in self.stations:
# if not self.ns['mim_phase'](src,station).initialized():
# sel=self.ns['select'](n)<<Meq.Selector(KL_node,index=n);
# self.ns['mim_phase'](src,station)<<Meq.Polar(1,-sel);
# n=n+1;
# return self.ns['mim_phase'];
    def make_tec(self):
        # fit a virtual TEC value, this makes life easier (eg. include freq. and sec dependence)
        pp = self.make_pp(ref_station=self.ref_station)
        KL_node = self.create_KL_node()
        n = 0
        for src in self.src:
            for station in self.stations:
                sec = self.ns['sec'](src, station)
                if not self.ns['tec'](src, station).initialized():
                    sel = self.ns['tec'](n) << Meq.Selector(KL_node, index=n)
                    self.ns['tec'](src, station) << sel * sec
                n = n + 1

        return self.ns['tec']

    def combine_parms(self):
        parms = self.ns['parms']()
        if parms.initialized():
            return parms
        parmlist = ()
        for i in range(rank):
            name = "KLParm:" + str(i)
            parmlist += (self._parm(name),)
        parms << Meq.Composer(children=parmlist)
        return parms

    def combine_pps(self):
        pps = self.ns['pps']()
        if pps.initialized():
            return pps
        pplist = ()
        pp = self.make_pp(ref_station=self.ref_station)
        dims = 0
        for src in self.src:
            for station in self.stations:
                pplist += (pp(src, station),)
                dims = dims + 1

        pps << Meq.Composer(children=pplist, dims=(dims, 3))
        return pps

    def create_KL_node(self):
        KL_node = self.ns['KL_node']()
        if KL_node.initialized():
            return KL_node
        parms = self.combine_parms()
        pps = self.combine_pps()
        kl = self.ns['KL_node']() << Meq.PyNode(children=(parms, pps), class_name="KLNode", module_name="Lions.PiercePoints.modules.KL.KLNode",
                                                rank=rank)
        return kl
