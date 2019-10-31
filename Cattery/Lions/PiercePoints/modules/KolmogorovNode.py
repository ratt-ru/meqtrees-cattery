# -*- coding: utf-8 -*-
# standard preamble
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from Timba.TDL import *
from Timba import pynode
from Timba import dmi
from Timba import utils
from Timba.Meq import meq
from Lions.PiercePoints.modules import PhaseScreen

Settings.forest_state.cache_policy = 100

_dbg = utils.verbosity(0, name='test_pynode')
_dprint = _dbg.dprint
_dprintf = _dbg.dprintf

initialized = False

# This class is meant to illustrate the pynode interface. All pynodes
# need to be derived from the Timba.pynode.PyNode class.


class KolmogorovNode (pynode.PyNode):
    # An __init__ method is not necessary at all. If you do define your
    # own for some reason, make sure you call pynode.PyNode.__init__()

    def __init__(self, *args):
        pynode.PyNode.__init__(self, *args)
        self.set_symdeps('domain')

    def update_state(self, mystate):
        global initialized
        # mystate is a magic "state helper" object is used both to set up
        # initial/default state, and to update state on the fly later.

        # This does the following:
        #  - checks the incoming state record for field 'a'
        #  - if present, sets self.a = staterec.a
        #  - if not present but we're initializing the node, sets self.a = 3,
        #    and also sets staterec.a = 3 on the C++ side
        #  - if not present and not initializing, does nothing
        mystate('scale', 200)
        # scale in km.
        mystate('speedx', 0.)
        # speed in m/s.
        mystate('speedy', 0.)
        # speed in m/s.
        mystate('grid_size', 10)
        # scale in km.
        mystate('tec0', 5.0)
        # beta.
        mystate('beta', 5. / 3.)
        # beta.
        mystate('amp_scale', 1e-5)
        # scale of amplitude.
        mystate('seed_nr', None)
        # scale of amplitude.
        mystate('starttime', 0.0)
        # scale of amplitude.

        if not initialized:
            # divide gridsize by 2 here, phasescreen produces twice the number of pixels
            PhaseScreen.init_phasescreen(self.grid_size, self.beta, self.seed_nr)
            initialized = True

    def get_result(self, request, *children):
        if len(children) < 1:
            raise TypeError("this is NOT a leaf node, At least 1  child with piercepoints expected!")
        res1 = children[0]
        vs1 = res1.vellsets
        # pierce_points, vector of length 2 or 3 (x,y(,z))
        vector_size = len(vs1)
        # for now use fist two:
        if vector_size < 2:
            raise TypeError("vector size of child 1 too small, at leat x/y expected")
        xv = vs1[0].value[0]
        yv = vs1[1].value[0]

        if vs1[0].has_field('shape'):
            shapex = vs1[0].shape
        else:
            shapex = (1,)
        if vs1[1].has_field('shape'):
            shapey = vs1[1].shape
        else:
            shapey = (1,)

        cells = request.cells
        seg = cells.segments.time
#        print '************************************************************'
#        print cells
#        print '------------------------------------------------------------'
#        print seg
#        print '************************************************************'
        # the startt and endt are the timeslots when tiling is set > 1
        if type(seg.start_index) == type(1):
            startt = seg.start_index
            endt = seg.end_index
        else:
            startt = seg.start_index[0]
            endt = seg.end_index[-1]
#        print '************************************************************'
#        print startt, endt
#        print '************************************************************'

        # make time a lot smaller to prevent precision errors for int
        # the actual value of the constant
        time = cells.grid.time - self.starttime

        if startt >= endt:
            time = [time, ]
        val = []
        for it in range(startt, endt + 1):
            if shapex[0] > 1:
                xv = vs1[0].value[it]
            if shapey[0] > 1:
                yv = vs1[1].value[it]

            xshift = (time[it]) * self.speedx
            yshift = (time[it]) * self.speedy
            xn = (xv + xshift) / self.scale
            yn = (yv + yshift) / self.scale
            xn = int(xn) % (self.grid_size)  # zorg dat xn een integer tussen 0 en grid_size is
            yn = int(yn) % (self.grid_size)  # zorg dat xn een integer tussen 0 en grid_size is

            val.append(PhaseScreen.phasescreen[xn][yn] * self.amp_scale + self.tec0)
        # fill result
        res = meq.result(None, cells)
        # print startt,endt,seg,val;
        val2 = meq.vells(shape=meq.shape(endt + 1,))
        val2[:] = val
        vs = meq.vellset(val2)
        res.vellsets = [vs, ]
        return res
