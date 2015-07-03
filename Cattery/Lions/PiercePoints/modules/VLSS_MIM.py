from Timba.TDL import *
from Timba.Meq import meq
from Lions.PiercePoints.PiercePoints import *
import Meow


def compile_options():
    return [TDLCompileOption("night", "VLSS nighttime Kolmogorov model", True),
            TDLCompileOption("day", "VLSS daytime Kolmogorov model", False)]


class MIM(PiercePoints):

    """Create MIM_model with Kolmogorov phase screen, default settings from VLSS observations"""

    def __init__(self, ns, name, sources, stations=None, ref_station=None, tags="iono", make_log=False):
        global beta
        global N
        global speedx
        global speedy
        global scale
        global amp_scale
        global TEC0
        global height
        global seed_nr
        global use_lonlat
        if night:
            beta = 5. / 3.
            N = 1024
            speedx = 400.
            speedy = 50.
            scale = 0.5
            amp_scale = 0.5
            TEC0 = 5.0
            height = 400.
            seed_nr = 310573
            use_lonlat = True
        if day:
            beta = 5. / 3.
            N = 1024
            speedx = 400.
            speedy = 50.
            scale = 0.5
            amp_scale = 1.0
            TEC0 = 5.0
            height = 300.
            seed_nr = 3418
            use_lonlat = True

        PiercePoints.__init__(self, ns, name, sources, stations, height, make_log)

        self.ref_station = ref_station
        global starttime
        starttime = 4453401230.76  # hardcode for now to make time run from ~0, get from MS later
        #        init_phasescreen(N,beta);

        # convert velocity and size to sensible values if PP in lonlat
        global vx
        global vy
        global pixscale
        if use_lonlat:
            print 'Using longlat'
            R = height + 6365.0  # hardcoded Earth radius to same value as other MIM
            # convert velocity in [km/h] to angular velocity
            vx = ((speedx / 3600) / R) + 7.272205E-05  # correct for PP velocity in longitude
            vy = (speedy / 3600) / R
            # convert size in [km] to angular size in radians
            pixscale = scale / R
        else:
            print 'Using ECEF'
            # for ECEF coordinates the velocity is not yet properly implemented
            vx = speedx / 3.6  # [m/s]
            vy = speedy / 3.6
            pixscale = scale * 1000.0  # [m]

        print 'V_x =', vx
        print 'V_y =', vy
        print 'Scale = ', pixscale

    def make_tec(self):
        # fit a virtual TEC value, this makes life easier (eg. include freq. and sec dependence)
        if use_lonlat:
            pp = self.make_longlat_vector_pp(ref_station=self.ref_station)
            pp = pp('longlat')
        else:
            pp = self.make_pp(ref_station=self.ref_station)
        n = 0
        for src in self.src:
            for station in self.stations:
                Kol_node = self.create_Kol_node(pp, src, station)
                sec = self.ns['sec'](src, station)
                if not self.ns['tec'](src, station).initialized():
                    self.ns['tec'](src, station) << Kol_node * sec
                n = n + 1

        return self.ns['tec']

    def create_Kol_node(self, pp, src, station):
        Kol_node = self.ns['Kol_node'](src, station)
        if Kol_node.initialized():
            return Kol_node
        kl = self.ns['Kol_node'](src, station) << Meq.PyNode(children=(pp(src, station),), class_name="KolmogorovNode", module_name="Lions.PiercePoints.modules.KolmogorovNode",
                                                             grid_size=N, beta=beta, scale=pixscale, speedx=vx, speedy=vy, amp_scale=amp_scale, seed_nr=seed_nr, tec0=TEC0, starttime=starttime)
        return kl
