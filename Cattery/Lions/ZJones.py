from Timba.TDL import *
import Meow
from Meow import Context
from Meow import Jones, ParmGroup, Bookmarks
from Lions.PiercePoints.modules.KL import KL_MIM
from Lions.PiercePoints.modules import TID_MIM
from Lions.PiercePoints.modules import Poly_MIM
from Lions.PiercePoints.modules import Kolmogorov_MIM
from Lions.PiercePoints.modules import VLSS_MIM


def _modname(obj):
    if hasattr(obj, 'name'):
        name = obj.name
    elif hasattr(obj, '__name__'):
        name = obj.__name__
    else:
        name = obj.__class__.__name__
    return name


def _modopts(mod, opttype='compile'):
    """for the given module, returns a list of compile/runtime options suitable for passing
    to TDLMenu. If the module implements a compile/runtime_options() method, uses that,
    else simply uses the module itself."""
    modopts = getattr(mod, opttype + '_options', None)
    # if module defines an xx_options() callable, use that
    if callable(modopts):
        return list(modopts())
    # else if this is a true module, it may have options to be stolen, so insert as is
    elif inspect.ismodule(mod):
        return [mod]
    # else item is an object emulating a module, so insert nothing
    else:
        return []

modules = [TID_MIM, Kolmogorov_MIM, VLSS_MIM, Poly_MIM, KL_MIM]
submenus = [TDLMenu("Use '%s' module" % _modname(mod), name=_modname(mod),
                    toggle=_modname(mod).replace('.', '_'), namespace={},
                    *_modopts(mod, 'compile'))
            for mod in modules]
mainmenu = TDLMenu("MIM model", exclusive="selname",
                   *(submenus))


class ZJones(object):

    def __init__(self):

        self.options = []

    def runtime_options(self):
        return self.options

    def compute_jones(self, jones, sources, stations=None, tags=None, label='', inspectors=[], **kw):
        stations = stations or Context.array.stations()
        print "selected", selname
        for mod in modules:
            if _modname(mod).replace('.', '_') == selname:
                self.mim_model = mod
                break
        mim = self.mim_model.MIM(jones.Subscope(), None, sources, Context.array, tags=tags, ref_station=ref_station, make_log=make_log)
        mim.compute_jones(jones, tags=tags)
        inspectors += mim.inspectors()
        qual_list = [stations]
        if sources:
            qual_list.insert(0, [src for src in sources if jones(src, stations[0]).initialized()])
        inspectors.append(
            jones.Subscope()['inspector'] << Meow.StdTrees.define_inspector(jones, *qual_list))
        return jones

    def compile_options(self):
        return [TDLOption("make_log", "Create Log Nodes", False), TDLOption("ref_station", "Rotate frame to reference station", [None, 1], more=int), mainmenu, ]
