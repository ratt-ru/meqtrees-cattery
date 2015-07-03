from Timba.TDL import *
from Timba.Meq import meq
from Lions.MIM_model import *

import Meow


class PiercePoints(MIM_model):

    """A PiercePoints is a MIM_model that uses pierce points at a given height (default=300km) which can be a (solvable) parameter. get_tec is a function of the calculated piercepoints"""

    def __init__(self, ns, name, sources, stations=None, height=None, make_log=False):
        MIM_model.__init__(self, ns, name, sources, stations)
        if isinstance(height, Meow.Parm):
            self._height = height
        else:
            self._height = Meow.Parm(height)

        self.ns.h << self._height .make()
        self._make_log = make_log

    def make_rot_matrix(self, ref_station=1):
        ns = self.ns
        xyz = self.array.xyz()
        # create rotation matrix"
        longlat = ns['ll'](ref_station)
        rot_matrix = ns['rot'](ref_station)
        if not rot_matrix.initialized():
            if not longlat.initialized():
                longlat << Meq.LongLat(xyz(ref_station), use_w=True)
                lon = ns['lon'](ref_station) << Meq.Selector(longlat, index=0)
                lat = ns['lat'](ref_station) << Meq.Selector(longlat, index=1)
                h_stat = ns['h_stat'](ref_station) << Meq.Selector(longlat, index=2)  # hm..not sure of we can use this if below earth surface!!Maybe find another way to find the earth_radius
                cosl = ns['coslon'](ref_station) << Meq.Cos(lon)
                sinl = ns['sinlon'](ref_station) << Meq.Sin(lon)
                cosphi = ns['coslat'](ref_station) << Meq.Cos(lat)
                sinphi = ns['sinlat'](ref_station) << Meq.Sin(lat)
                # create rotation matrix to convert to/from ENU coordinates
                rot_matrix << Meq.Composer(-1 * sinl, cosl, 0,
                                           -1 * sinphi * cosl, -1 * sinphi * sinl, cosphi,
                                           cosphi * cosl, cosphi * sinl, sinphi,
                                           dims=[3, 3])
        return rot_matrix

    def make_pp(self, ref_station=None):
        # returns xyz position of pierce point, assumes spherical earth, if elliptical some other method should be chosen. azel seems to be related to spherical earth plus alpha_prime and the last formula calculating the scale.
        ns = self.ns
        if not ns.pi.initialized():
            ns.pi << Meq.Constant(math.pi)

        xyz = self.array.xyz()
        for station in self.stations:
            rot_matrix = ns['rot'](station)
            if not rot_matrix.initialized():
                self.make_rot_matrix(station)
        if ref_station:
            ref_rot_matrix = ns['rot'](ref_station)

        for station in self.stations:
            norm_xyz = create_inproduct(ns, xyz(station), xyz(station))
            earth_radius = ns.earth_radius(station) << Meq.Subtract(norm_xyz, ns['h_stat'](station))
            rot_matrix = ns['rot'](station)
            for src in self.src:
                az = self.get_az()(src, station)
                el = self.get_el()(src, station)
                cos_az = Meq.Cos(az)
                sin_az = Meq.Sin(az)
                cos_el = Meq.Cos(el)
                sin_el = Meq.Sin(el)
                diff = ns['diff'](src, station)
                diff_vector = diff
                if not diff.initialized():
                    diff_vector = diff << Meq.Composer(Meq.MatrixMultiply(Meq.Composer(cos_el * sin_az,
                                                                          cos_el * cos_az,
                                                                          sin_el, dims=[1, 3]),
                                                                          rot_matrix))

                alpha_prime = Meq.Asin(cos_el * norm_xyz / (earth_radius + ns.h * 1000.))
                # correct. draw triangle from piercepoint to center earth, and along the line connecting pp and antenna with 90 degree angle between center earth and that line.
                sec = ns.sec(src, station) << 1. / Meq.Cos(alpha_prime)
                sin_beta = ns.sin_beta(src, station) << Meq.Sin((0.5 * ns.pi - el) - alpha_prime)
                # angle at center earth
                scale = ns.scale(src, station) << (earth_radius + ns.h * 1000.) * sin_beta / cos_el
                # correct. draw triangle from center earth to pierce point and along vector a (=xyz_station), with a 90 degree angle between that line and piercepoint

                if ref_station:
                    ns['pp'](src, station) << Meq.MatrixMultiply(ref_rot_matrix, xyz(station) + diff_vector * scale)
                else:
                    ns['pp'](src, station) << xyz(station) + diff_vector * scale

        return ns['pp']

    def make_longlat_pp(self, ref_station=None):
        '''make longitude and lattitude of piercepoints'''
        pp = self.make_xyz_pp(ref_station=ref_station)
        for station in self.stations:
            for src in self.src:
                x = pp('x', src, station)
                y = pp('y', src, station)
                z = pp('z', src, station)
                pp('lon', src, station) << Meq.Atan(x / y)
                pp('lat', src, station) << Meq.Asin(z / Meq.Sqrt(x * x + y * y + z * z))
        return pp

    def make_longlat_vector_pp(self, ref_station=None):
        '''combine the long lat of the piercepoint in a vector, to use them in the same way as the xy(z) vector that is pp(src,station)'''
        pp = self.make_longlat_pp(ref_station=ref_station)
        for station in self.stations:
            for src in self.src:
                longlatvector = pp('longlat', src, station)
                if not longlatvector.initialized():
                    longlatvector << Meq.Composer(pp('lon', src, station), pp('lat', src, station))
        return pp

    def make_xy_pp(self, ref_station=None):
        '''make xy of piercepoints'''
        pp = self.make_pp(ref_station=ref_station)
        for station in self.stations:
            for src in self.src:
                x = pp('x', src, station)
                y = pp('y', src, station)
                if not x.initialized():
                    x << Meq.Selector(pp(src, station), index=0)
                if not y.initialized():
                    y << Meq.Selector(pp(src, station), index=1)
        return pp

    def make_xyz_pp(self, ref_station=None):
        '''make xyz of piercepoints'''
        pp = self.make_xy_pp(ref_station=ref_station)
        for station in self.stations:
            for src in self.src:
                z = pp('z', src, station)
                if not z.initialized():
                    z << Meq.Selector(pp(src, station), index=2)

        return pp

    def create_log_nodes(self, xy=False):  # create nodes to log phase per x,y (or long/lat if xy=False)
        if xy:
            pp = self.make_xy_pp()
        else:
            pp = self.make_longlat_pp()

        phases = self.make_phase_error()
        for station in self.stations:
            for src in self.src:
                filename = "phases_" + str(station) + "_" + src.name + ".dat"
                log = pp('log', src, station)
                if not log.initialized():
                    if xy:
                        log << Meq.PyNode(children=(pp('x', src, station), pp('y', src, station), phases(src, station)),
                                          class_name="PrintPyNode", module_name="Lions.PrintPyNode", filename=filename)
                    else:
                        log << Meq.PyNode(children=(pp('lon', src, station), pp('lat', src, station)),
                                          class_name="PrintPyNode", module_name="Lions.PrintPyNode", filename=filename)
        return pp('log')

    def create_station_log_nodes(self):  # create nodes to log phase of all stations per source
        phases = self.make_phase_error()
        ns = self.ns
        for src in self.src:
            direction = src.direction
            pd = direction._parmdefs
            if isinstance(direction, Meow.LMDirection):
                nm = 'l.' + str(pd['l'][0]) + '.m.' + str(pd['m'][0])
            else:
                nm = 'ra.' + str(pd['ra'][0]) + '.dec.' + str(pd['dec'][0])
            print "creating file for", nm
            for station in self.stations:
                filename = "phases_" + nm + "_" + str(station) + ".dat"
                log = ns['station_log'](src, station)
                if not log.initialized():
                    log << Meq.PyNode(children=[phases(src, station), ],
                                      class_name="PrintPyNode", module_name="Lions.PrintPyNode", filename=filename)
        return ns['station_log']

    def inspectors(self):
        inspectors = []
        if self._make_log:
            # ip=self.ns.inspector<<Meow.StdTrees.define_inspector(self.create_log_nodes(),self.src,self.stations);
            ip = self.ns['inspector']('lognodes') << Meow.StdTrees.define_inspector(self.create_station_log_nodes(), self.src, self.stations)
            inspectors.append(ip)
        print "INSPECTORS:::::", inspectors
        return inspectors


def create_inproduct(ns, a, b, length=0):
    """Computes the dot product of two vectors of arbitrary length"""
    # Definition of dot product: multiply vector elements separately?
    # I.e. A=(a1,a2,a3), B=(b1,b2,b3) ==> A dot B = (a1*b1 + a2*b2 + a3*b3)
    # For now assume that A and B have equal number of elements
    # The Meq.NElements returns a node_stub which is not accepted by the range()
    # so for now only use 3-element vectors.
    # n_elements = ns.Aelements << Meq.NElements(a)
    # we can try to get length from the number of children but this only works if you dotproduct composers...
    if ns.dot(a.name, b.name).initialized():
        return ns.dot(a.name, b.name)

    if length == 0:
        length = a.num_children()
        length_b = b.num_children()
        if length != length_b:
            print "Vectors must be of same length!!!, using smallest"
            length = min(length, length_b)
    sumdot = 0
#    print "creating dot product of length",length;
    for i in range(length):  # n_elements is a node_stub, cannot be used for loop counts
        if not a('select', i).initialized():
            a('select', i) << Meq.Selector(a, index=i)
        first = a('select', i)
        if not b('select', i).initialized():
            b('select', i) << Meq.Selector(b, index=i)
        second = b('select', i)
        sumdot = sumdot + (first * second)
    dot = ns.dot(a.name, b.name) << Meq.Sqrt(sumdot)
    return dot
