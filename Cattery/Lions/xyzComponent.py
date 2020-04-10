from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from Timba.TDL import *
from Timba.Meq import meq
from Meow.Parameterization import *
import Meow


class xyzComponent (Parameterization):

    """A xyzComponent represents an abstract ITRF element.
    xyzComponents have an name and an associated xyz vector in ITRF coordinates. These can be eg. satellite or station positions.
    """

    def __init__(self, ns, name, position):
        Parameterization.__init__(self, ns, name)
        self.xyz = position

    def make_rot_matrix(self):
        """returns the rotation matrix to ENU with this xyz position as referencepoint. NOTE for the moment this can only be done with constant xyz positions (restriction of longlat)"""
        if not self.ns.rot_matrix(self.name).isdefined():
            lon, lat = self.make_longlat()
            cosl = self.ns.cosl(self.name) << Meq.Cos(lon)
            sinl = self.ns.sinl(self.name) << Meq.Sin(lon)
            cosphi = self.ns.cosphi(self.name) << Meq.Cos(lat)
            sinphi = self.ns.sinphi(self.name) << Meq.Sin(lat)

            self.ns.rot_matrix(self.name) << MeqComposer(-1 * sinl, cosl, 0,
                                                         -1 * sinphi * cosl, -1 * sinphi * sinl, cosphi,
                                                         cosphi * cosl, cosphi * sinl, sinphi,
                                                         dims=[3, 3])

        rot_matrix = self.ns.rot_matrix(self.name)
        return rot_matrix

    def make_longlat(self, use_w=1):
        if not self.ns.longlat(self.name).isdefined():

            ll = self.ns.longlat(self.name) << Meq.LongLat(self.xyz, use_w=use_w)
            # use_w=1: use elliptical model for earth

            lon = self.ns.lon(self.name) << Meq.Selector(ll, index=0)
            lat = self.ns.lat(self.name) << Meq.Selector(ll, index=1)

        lon = self.ns.lon(self.name)
        lat = self.ns.lat(self.name)
        return lon, lat
