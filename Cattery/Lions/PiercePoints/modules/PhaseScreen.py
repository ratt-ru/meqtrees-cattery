# just a trick to import the phasescreen
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

import math
import numpy
from numpy.fft import *
from numpy.random import *

phasescreen = []


def init_phasescreen(N=10, beta=5., seed_nr=None):
    global phasescreen
    # X,Y are matrices containing the x and y coordinates
    X = numpy.matrix(numpy.ones((2 * N, 1))) * numpy.matrix(list(range(-N, N))) * 1.0 / N
    Y = numpy.matrix(list(range(-N, N))).T * numpy.matrix(numpy.ones((1, 2 * N))) * 1.0 / N
    # Q is the distance from the origin
    Q = numpy.sqrt(numpy.power(X, 2) + numpy.power(Y, 2))

    # start from the same random seed to get identical screen
    if seed_nr is not None:
        # print "seeding with",seed_nr;
        seed(seed_nr)

    # W is complex white noise
    W = (standard_normal((2 * N, 2 * N))) * numpy.exp(1j * uniform(-math.pi, math.pi, (2 * N, 2 * N)))

    # Shape white noise by multiplying it by Q^(-2-beta)
    Q[N, N] = 1
    S = numpy.multiply(W, numpy.sqrt(numpy.power(Q, -2 - beta)))
    S[N, N] = 0

    # Compute the inverse real ffts in both directions
    screen = numpy.real(ifft2(fftshift(S)))

    # Normalize phasescreen
    # print "in PhaseScreen",len(screen);
    minscreen = numpy.min(screen)
    maxscreen = numpy.max(screen)

    # get values between 0 and 1:
    # phasescreen = (screen / (2*maxscreen)) + 0.5 ;
    # alternative get values between 0 and 1
    # screen = screen - minscreen;
    # phasescreen = screen / maxscreen;
    # get values between -1 and 1
    phasescreen = (screen / maxscreen)
    numpy.save('phase_screen_display', phasescreen)
