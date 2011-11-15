# -*- coding: utf-8 -*-

import math
import numpy
from scipy.ndimage import interpolation
import Kittens.utils

_verbosity = Kittens.utils.verbosity(name="vb");
dprint = _verbosity.dprint;
dprintf = _verbosity.dprintf;

DEG = math.pi/180;

def expand_axis (x,axis,n):
  """Expands an array to N elements along the given axis. Array must have
  1 element along the given axis (or have fewer axes)"""
  # if array has fewer axes, pad with size-1 axes
  if x.ndim <= axis:
    x = x.reshape(list(x.shape)+[1]*(axis-x.ndim+1));
  if x.shape[axis] == n:
    return x;
  elif x.shape[axis] != 1:
    raise TypeError,"array must have length 1 along axis %d, it has %d"%(axis,x.shape[axis]);
  return numpy.concatenate([x]*n,axis);

def unite_shapes (a,b):
  """Makes two arrays have the same shape, as follows:
  - each axis must be the same shape,
  - or if one array has shape-1, and the other shape-N, then the shape-1 is expanded.
  """
  if a.shape == b.shape:
    return a,b;
  # promote shapes to same number of dinensions, by padding missing dimensions with 1's
  sa = [1]*max(a.ndim,b.ndim);
  sb = list(sa);
  sa[:a.ndim] = a.shape;
  sb[:b.ndim] = b.shape;
  # reshape arrays
  a = a.reshape(sa);
  b = b.reshape(sb);
  # now loop, and expand missing axes
  for axis,(na,nb) in enumerate(zip(sa,sb)):
    if na != nb:
      if na == 1:
        a = expand_axis(a,axis,nb);
      elif nb == 1:
        b = expand_axis(b,axis,na);
      else:
        raise TypeError,"error: trying to unite incompatible shapes %s and %s"%(sa,sb);
  return a,b;



class InterpolatedVoltageBeam (object):
  """This class implements a complex (interpolated) voltage beam as a function of LM."""
  def __init__ (self,spline_order=2,ampl_interpolation=True,l0=0,m0=0):
    self._spline_order = spline_order;
    self._ampl_interpolation = ampl_interpolation
    self.l0,self.m0 = l0,m0;

  def hasFrequencyAxis (self):
    return bool(self._freqToBeam);
    
  def freqGrid (self):
    return self._freq_grid;
    
  def setFreqGrid (self,freqs):
    self._freq_grid = freqs;
    
  def freqToBeam (self,freq):
    """Maps frequencies (could be a list or an array) to "natural" coordinate system of the beam""";
    return self._freqToBeam(freq);
    
  def lmToBeam (self,l,m):
    """Maps l/m coordinates onto the "natural" coordinate system of the beam""";
    return self._lmToBeam(l,m);
    
  def setBeamCube (self,beam,lmap=None,mmap=None,lmmap=None,freqmap=None):
    """Initializes beam cube. Beam is a 2D (or 3D) cube. lmap, mmap, freqmap
    are callables used to convert l/m/freq coordinates into coordinates within the cube""";
    if not lmmap:
      if lmap and mmap:
        lmmap = lambda l,m,lm=lmap,mm=mmap:lm(m),mm(m);
      else:
        raise ValueError,"Either lmmap, or both lmap and mmap need to be supplied";
    self._lmToBeam = lmmap;
    self._freqToBeam = freqmap;
    # init cube 
    self._beam = beam;
    if self._ampl_interpolation:
      self._beam_ampl = abs(beam);
    else:
      self._beam_ampl = None;
    # replace with pre-filtrered arrays if needed 
    if self._spline_order > 1:
      self._beam_real = interpolation.spline_filter(beam.real,order=self._spline_order);
      self._beam_imag = interpolation.spline_filter(beam.imag,order=self._spline_order);
      if self._beam_ampl is not None:
        self._beam_ampl = interpolation.spline_filter(self._beam_ampl,order=self._spline_order);
    else:
      self._beam_real = beam.real;
      self._beam_imag = beam.imag;
    
  def transformCoordinates (self,l,m,freq=None,time=None,freqaxis=None,timeaxis=None):
    """Transforms sets of l/m (and freq/time, etc.) coordinates into an array of 
    local in-beam coordinates. 
    l,m may be arrays (both must be the same shape, or will be promoted to the same shape)
    If beam has a freq dependence, then an array of frequency coordinates (freq) must be given,
    and freqaxis must be set to the number of the frequency axis. Then the following
    possibilities apply:

    (A) len(freq)==1  (freqaxis need not be set)
        l/m is interpolated at the same frequency point.
        Output array is same shape as l/m.

    (B) len(freq)>1 and l.shape[freqaxis] == 1:
        The same l/m is interpolated at every point in freq.
        Output array is same shape as l/m, plus an extra frequency axis (number freqaxis)

    (C) len(freq)>1 and l.shape[freqaxis] == len(freq)
        l/m has its own freq dependence, so a different l/m/freq is interpolated at every point.
        Output array is same shape as l/m.

    And finally (D):

    (D) No dependence on frequency in the beam.
        We simply interpolate every l/m value as is. Output array is same shape as l/m.
        
    Return value is a tuple of coords,shape, where 'coords' is an (N,M) array suitable for 
    feeding to e.g. map_coordinates() (N=2 for an lm-only beam, 3 for an l-m-freq beam, etc.),
    and 'shape' is the shape of the output array (where the product of all elements in 'shape'
    is M).

    'time' is currently ignored -- provided for later compatibility (i.e. beams with time planes)
    """
    # make sure inputs are arrays
    l = numpy.array(l) + self.l0;
    m = numpy.array(m) + self.m0;
    freq = numpy.array(freq);
    # promote l,m to the same shape
    l,m = unite_shapes(l,m);
    dprint(3,"input l/m shape is",l.shape,m.shape);
    dprint(3,"input l/m[0] is",l.ravel()[0],m.ravel()[0]);
#    dprint(4,"input l/m is",l,m);
    l,m = self.lmToBeam(l,m);
    dprint(3,"in beam coordinates this is",l.ravel()[0],m.ravel()[0]);
#    dprint(4,"in beam coordinates this is",l,m);
    # now we make a 3xN coordinate array for map_coordinates
    # lm[0,:] will be flattened L array, lm[1,:] will be flattened M array
    # lm[2,:] will be flattened freq array (if we have a freq dependence)
    # Do we have a frequency axis in the beam? (case A,B,C):
    if self.hasFrequencyAxis():
      if freq is None:
        raise ValueError,"frequencies not specified, but beam has a frequency dependence";
      freq = numpy.array(freq);
      if not freq.ndim:
        freq = freq.reshape(1);
      dprint(4,"input freq is",freq);
      freq = self._freqToBeam(freq);
      dprint(4,"in beam coordinates this is",freq);
      # case (A): reuse same frequency for every l/m point
      if len(freq) == 1:
        lm = numpy.vstack((l.ravel(),m.ravel(),[freq[0]]*l.size));
      # case B/C:
      else:
        # first turn freq vector into an array of the proper shape
        if freqaxis is None:
          raise ValueError,"frequency axis not specified, but beam has a frequency dependence";
        freqshape = [1]*(freqaxis+1);
        freqshape[freqaxis] = len(freq);
        freq = freq.reshape(freqshape);
        # now promote freq to same shape as l,m. This takes care of cases B and C.
        # (the freq axis of l/m will be expanded, and other axes of freq will be expanded)
        l,freq = unite_shapes(l,freq);
        m,freq = unite_shapes(m,freq);
        lm = numpy.vstack((l.ravel(),m.ravel(),freq.ravel()));
    # case (D): no frequency dependence in the beam
    else:
      lm = numpy.vstack((l.ravel(),m.ravel(),[0]*l.size));
    return lm,l.shape;

  def interpolate (self,l,m,time=None,freq=None,freqaxis=None,output=None):
    """Interpolates beam into the given l/m (and freq/time, etc.) coordinates.
    See transformCoordinates() above for a description of the parameters. An output array
    may be specified (it will be resized to the correct shape just in case).
    Returns array of interpolated coordinates.
    """;
    # transform coordinates
    dprint(3,"transforming coordinates");
    coords,output_shape = self.transformCoordinates(l,m,time=time,freq=freq,freqaxis=freqaxis);
    # prepare output array
    if output is None:
      output = numpy.zeros(output_shape,complex);
    elif output.shape != output_shape:
      output.resize(output_shape);
    dprint(3,"interpolating %s coordinate points to output shape %s"%(coords.shape,output_shape));
    # interpolate real and imag parts separately
    output.real = interpolation.map_coordinates(self._beam_real,coords,order=self._spline_order,
                  prefilter=(self._spline_order==1)).reshape(output_shape);
    output.imag = interpolation.map_coordinates(self._beam_imag,coords,order=self._spline_order,
                  prefilter=(self._spline_order==1)).reshape(output_shape);
    if not self._beam_ampl is None:
      output_ampl = interpolation.map_coordinates(self._beam_ampl,coords,order=self._spline_order,
                  prefilter=(self._spline_order==1)).reshape(output_shape);
      phase_array = numpy.arctan2(output.imag,output.real)
      output.real = output_ampl * numpy.cos(phase_array)
      output.imag = output_ampl * numpy.sin(phase_array)
    dprint(3,"interpolated value [0] is",output.ravel()[0]);
    # dprint(4,"interpolated value is",output);
    return output;
