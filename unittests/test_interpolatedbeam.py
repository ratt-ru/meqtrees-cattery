from nose.tools import *
import numpy as np
import time
from subprocess import Popen
import os

import Cattery.Siamese.OMS.Utils as Utils
import Cattery.Siamese as Siamese
import Cattery.Siamese.OMS.InterpolatedBeams as InterpolatedBeams

PREFIX = "/test/test_beam"
PATTERN = PREFIX + "_$(corr)_$(reim).fits"

def __run(cmdargs, timeout=600):
    cmd = " ".join(cmdargs)
    print(f"Running {cmd} with timeout {timeout} seconds")
    p = Popen(cmdargs, 
                env=os.environ.copy())

    x = timeout
    delay = 1.0
    timeout = int(x / delay)
    while p.poll() is None and timeout > 0:
        time.sleep(delay)
        timeout -= delay
    #timeout reached, kill process if it is still rolling
    ret = p.poll()
    if ret is None:
        p.kill()
        ret = 99

    if ret == 99:
        raise RuntimeError("Test timeout reached. Killed process.")
    elif ret != 0:
        raise RuntimeError("{} exited with non-zero return code {}".format(cmdargs[0], ret))

REALIMAG = dict(re="real",im="imag")

def __make_beam_filename (filename_pattern,corr,reim,stationtype):
    """Makes beam filename for the given correlation and real/imaginary component (one of "re" or "im")"""
    return Utils.substitute_pattern(filename_pattern,
                corr=corr.lower(),xy=corr.lower(),CORR=corr.upper(),XY=corr.upper(),
                reim=reim.lower(),REIM=reim.upper(),ReIm=reim.title(),
                realimag=REALIMAG[reim].lower(),REALIMAG=REALIMAG[reim].upper(),
                RealImag=REALIMAG[reim].title(), 
                stype=stationtype.lower(), STYPE=stationtype.upper())

def testinterpolatedbeam():
    print("creating beams...")
    eidosargs = (f"eidos -d 0.5 -r 0.015625 -f 950 1050 20 -P {PREFIX} -o8").split(" ")
    __run(eidosargs, timeout=900)
    filenames = __make_beam_filename(PATTERN, "xy", "re", "default"), \
                __make_beam_filename(PATTERN, "xy", "im", "default")
    print("loading beam patterns %s %s" % filenames)
    vb = InterpolatedBeams.LMVoltageBeam(
        verbose=10000,
        l_axis="px", 
        m_axis="py"
    )
    vb.read(*filenames)
    print("successfully loaded beam patterns %s %s" % filenames)
    # approx small angle approximation
    vel = vb.interpolate(np.deg2rad(0.05),np.deg2rad(0.05),freq=1000e6,freqaxis=1)
    assert np.isclose(vel, np.array([-0.00091718+0.00020009j])), \
        "Last known good value from Eidos beams interpolated value does not " \
        "match current interpolated value"