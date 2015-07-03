# standard preamble
from Timba.TDL import *
from Timba import pynode
from Timba import dmi
from Timba import utils
from Timba.Meq import meq
from numpy import *
Settings.forest_state.cache_policy = 100

_dbg = utils.verbosity(0, name='test_pynode')
_dprint = _dbg.dprint
_dprintf = _dbg.dprintf

# This class is meant to illustrate the pynode interface. All pynodes
# need to be derived from the Timba.pynode.PyNode class.


class PrintPyNode (pynode.PyNode):

    # An __init__ method is not necessary at all. If you do define your
    # own for some reason, make sure you call pynode.PyNode.__init__()
    # as done below
    def __init__(self, *args):
        pynode.PyNode.__init__(self, *args)
        self.set_symdeps('domain')

    def update_state(self, mystate):
        # mystate is a magic "state helper" object is used both to set up
        # initial/default state, and to update state on the fly later.

        # This does the following:
        #  - checks the incoming state record for field 'a'
        #  - if present, sets self.a = staterec.a
        #  - if not present but we're initializing the node, sets self.a = 3,
        #    and also sets staterec.a = 3 on the C++ side
        #  - if not present and not initializing, does nothing
        mystate('filename', 'phases.dat')

    # leaf node), and a request object
    def get_result(self, request, *children):
        filen = open(self.filename, "wb")
        # add data to existing file
        for ch in children:
            for result in ch.vellsets:
                f = array(result.value)
                f.tofile(file=filen, format='%.6f ', sep=' : ')
                 # f.write(result.value);

        # f.close();
        return children[0]
