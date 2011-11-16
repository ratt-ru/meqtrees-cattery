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
  def __init__ (self,hier_interpol=True,spline_order=2,l0=0,m0=0):
    self._spline_order = spline_order;
    self.l0,self.m0 = l0,m0;
    self._hier_interpol = hier_interpol;
    self.interpolate = self.interpolate_linfreq if hier_interpol else self.interpolate_3d;

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
    self._beam = beam.copy();
    self._beam_ampl = abs(beam);
    self._beam_real = beam.real;
    self._beam_imag = beam.imag;
    # replace with pre-filtered arrays if needed 
    if beam.ndim > 2 and self._spline_order > 1:
      if not self._hier_interpol:
        self._beam_real = interpolation.spline_filter(beam.real,order=self._spline_order);
        self._beam_imag = interpolation.spline_filter(beam.imag,order=self._spline_order);
        self._beam_ampl = interpolation.spline_filter(self._beam_ampl,order=self._spline_order);
      else:
        # prefilter per frequency plane in per-frequency mode
        print "pre-filtering per plane",beam.shape;
        for i in range(beam.shape[2]):
          self._beam_real[...,i] = interpolation.spline_filter(beam.real[...,i],order=self._spline_order);
          self._beam_imag[...,i] = interpolation.spline_filter(beam.imag[...,i],order=self._spline_order);
          self._beam_ampl[...,i] = interpolation.spline_filter(self._beam_ampl[...,i],order=self._spline_order);
    
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
        freq = [0];
      else:
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

  def interpolate_3d (self,l,m,time=None,freq=None,freqaxis=None,output=None):
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
    global _first;
    if not _first:
      l0,m0,freq0 = coords[:,0];
      print "First interpolation point is at",l0,m0,freq0;
      print "Frequencies are",self._freq_grid;
      for l in int(l0),int(l0)+1:
        for m in int(m0),int(m0)+1:
          print l,m,"amplitudes",abs(self._beam[l,m,:]);
          print l,m,"phases",numpy.angle(self._beam[l,m,:])/DEG;
    
    dprint(3,"interpolating %s coordinate points to output shape %s"%(coords.shape,output_shape));
    # interpolate real and imag parts separately
    output.real = interpolation.map_coordinates(self._beam_real,coords,order=self._spline_order,
                    prefilter=(self._spline_order==1)).reshape(output_shape);
    output.imag = interpolation.map_coordinates(self._beam_imag,coords,order=self._spline_order,
                    prefilter=(self._spline_order==1)).reshape(output_shape);
    output_ampl = interpolation.map_coordinates(self._beam_ampl,coords,order=self._spline_order,
                    prefilter=(self._spline_order==1)).reshape(output_shape);
    phase_array = numpy.angle(output);
    output.real = output_ampl * numpy.cos(phase_array);
    output.imag = output_ampl * numpy.sin(phase_array);
    
    dprint(3,"interpolated value [0] is",output.ravel()[0]);
    # dprint(4,"interpolated value is",output);
    if not _first:
      _first = l0,m0;
      for i in range(coords.shape[1]):
        if coords[0,i] == l0 and coords[1,i] == m0:
          print "Interpolated amplitude at frequency %f is %f"%(coords[2,i],abs(output.ravel()[i]));
    
    return output;
    
  def _interpolate_freqplane (self,ifreq,coords,verbose=False):
    val = self._freqplanes.get(ifreq);
    if val is None:
      output = numpy.zeros(self._freqplane_shape,complex);
      # interpolate real and imag parts separately
      output.real = interpolation.map_coordinates(self._beam_real[...,ifreq],coords[:2,...],
                      order=self._spline_order,prefilter=(self._spline_order==1)).reshape(self._freqplane_shape);
      output.imag = interpolation.map_coordinates(self._beam_imag[...,ifreq],coords[:2,...],
                      order=self._spline_order,prefilter=(self._spline_order==1)).reshape(self._freqplane_shape);
      output_ampl = interpolation.map_coordinates(self._beam_ampl[...,ifreq],coords[:2,...],
                      order=self._spline_order,prefilter=(self._spline_order==1)).reshape(self._freqplane_shape);
      if verbose:
        l0,m0,freq = coords[:,0];
        for l in int(l0),int(l0)+1:
          for m in int(m0),int(m0)+1:
            print l,m,freq,"beam plane amplitude",abs(self._beam[l,m,ifreq]),self._beam_ampl[l,m,ifreq];
        print "interpolated amplitude",output_ampl[0];
      phase_array = numpy.angle(output);
      output.real = output_ampl * numpy.cos(phase_array);
      output.imag = output_ampl * numpy.sin(phase_array);
      val = self._freqplanes[ifreq] = output,output_ampl;
    return val;
    

  def interpolate_linfreq (self,l,m,freq,time=None,freqaxis=None,output=None):
    """Interpolates beam into the given l/m (and time, etc.) coordinates, with a fixed frequency grid.
    Interpolates each frequency plane separately, then does a linear interpolation between the 
    frequency planes.
    See transformCoordinates() above for a description of the parameters. An output array
    may be specified (it will be resized to the correct shape just in case).
    Returns array of interpolated coordinates.
    """;
    # transform coordinates
    dprint(3,"transforming coordinates");
    coords,output_shape = self.transformCoordinates(l,m,time=time,freq=None,freqaxis=freqaxis);
    freqcoord = self._freqToBeam(freq); 
    self._freqplanes = {};
    self._freqplane_shape = output_shape;
    # this will hold the interpolated per-frequency result
    output_shape = list(output_shape) + [len(freqcoord)];
    # prepare output array
    if output is None:
      output = numpy.zeros(output_shape,complex);
    elif output.shape != output_shape:
      output.resize(output_shape);
    global _print_first;
    verbose = _print_first and _verbosity.get_verbose()>0;
    if verbose:
      l0,m0,freq0 = coords[:,0];
      verbose = l0,m0;
      print "First interpolation point is at",l0,m0;
      print "Frequencies are",self._freq_grid;
      for l in int(l0),int(l0)+1:
        for m in int(m0),int(m0)+1:
          print l,m,"amplitudes",abs(self._beam[l,m,:]);
          print l,m,"phases",numpy.angle(self._beam[l,m,:])/DEG;
    
    dprint(3,"interpolating %s coordinate points to output shape %s"%(coords.shape,output_shape));
    for ifreq,freq in enumerate(freqcoord):
      if freq <= 0:
        output[...,ifreq] = self._interpolate_freqplane(0,coords)[0];
      elif freq >= len(self._freq_grid)-1:
        output[...,ifreq] = self._interpolate_freqplane(len(self._freq_grid)-1,coords,)[0];
      else:
        f0,f1 = int(freq),int(freq)+1;
        (c0,abs0),(c1,abs1) = self._interpolate_freqplane(f0,coords,
                  verbose and not ifreq),self._interpolate_freqplane(f1,coords,verbose and not ifreq);
        if verbose and not ifreq:
          print "per-plane weights",f0,f1;
          print "per-plane amplitudes",abs0[0],abs1[0];
          print "per-plane phases",numpy.angle(c0)/DEG,numpy.angle(c1)/DEG;
        ax = abs0*(f1-freq) + abs1*(freq-f0);
        cx = c0*(f1-freq) + c1*(freq-f0);
        px = numpy.angle(cx);
        output[...,ifreq].real = ax*numpy.cos(px);
        output[...,ifreq].imag = ax*numpy.sin(px);
    
    dprint(3,"interpolated value [0] is",output.ravel()[0]);
    # dprint(4,"interpolated value is",output);
    if verbose:
      _print_first = None;
      for ifreq,freq in enumerate(freqcoord):
        print "Frequencies",freqcoord;
        print "Interpolated amplitudes",abs(output[0,...]);
        print "Interpolated phases",numpy.angle(output[0,...])/DEG;
    
    return output;

# set to True to print the first interpolation
_print_first = True;