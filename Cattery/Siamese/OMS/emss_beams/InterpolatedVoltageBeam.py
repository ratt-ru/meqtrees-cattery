# -*- coding: utf-8 -*-

import math
import numpy
from scipy.ndimage import interpolation
import Kittens.utils

_verbosity = Kittens.utils.verbosity(name="vb");
#_verbosity.set_verbose(3)
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
  # promote shapes to same number of dimensions, by padding missing dimensions with 1's
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

def unite_multiple_shapes (*arr0):
  """Makes two or more arrays have the same shape, as follows:
  - each axis must be the same shape,
  - or if one array has shape-1, and another shape-N, then the shape-1 is expanded.
  Argument list may contain Nones, which will be returned as-is.
  """
  # less than two shapes, or all shapes are equal: return as is
  if len(arr0) < 2 or all([x.shape == arr0[0].shape for x in arr0[1:] if x is not None ]):
    return arr0;
  # promote shapes to same number of dimensions, by padding missing dimensions with 1's
  ndim = max([x.ndim for x in arr0 if x is not None]);
  sh0 = [1]*ndim;
  arr = [];      # output arrays
  for x in arr0:
    if x is None:
      arr.append(None);
    else:
      sh = list(sh0);
      sh[:x.ndim] = x.shape;
      arr.append(x.reshape(sh));
  # now loop over axes, and expand missing ones
  for axis in range(ndim):
    nx = [ (x.shape[axis] if x is not None else 0) for x in arr ];
    nx0 = max(nx);
    if nx0 != 1:
      for i,x in enumerate(arr):
        if x is not None:
          if nx[i] != nx0:
            if nx[i] == 1:
              arr[i] = expand_axis(arr[i],axis,nx0);
            else:
              raise TypeError,"error: trying to unite incompatible shapes: %s"%" ".join(
                [ "[%s]"%(",".join(map(str,x.shape))) if x is not None else "[]" for x in arr0 ]);
  return arr;

class InterpolatedVoltageBeam (object):
  """This class implements a complex (interpolated) voltage beam as a function of LM."""
  def __init__ (self,hier_interpol=True,spline_order=2,l0=0,m0=0):
    self._spline_order = spline_order;
    self.l0,self.m0 = l0,m0;
    self._hier_interpol = hier_interpol;
    self.interpolate = self.interpolate_linfreq if hier_interpol else self.interpolate_3d;

  def hasFrequencyAxis (self):
    return False;
    
  def freqGrid (self):
    return self._freq_grid;
    
  def setFreqGrid (self,freqs):
    self._freq_grid = freqs;
    
  def freqToBeam (self,freq):
    """Maps frequencies (could be a list or an array) to "natural" coordinate system of the beam.
    Returns array of coordinates. Negative values such as indicate extrapolation. 
    For a value of -X, X<1 indicates extrapolation to lower freqs, X>1 to higher freqs
    """;
    raise TypeError,"frequency interpolation not available with this beam gridder";
    
  def lmToBeam (self,l,m,rotate=None):
    """Maps l/m coordinates (i.e. direction cosines) onto the "natural" coordinate system of the beam""";
    raise TypeError,"l/m coordinate interpolation not available with this beam gridder";
    
  def thetaPhiToBeam (self,theta,phi,rotate=None):
    """Maps spherical theta/phi coordinates onto the "natural" coordinate system of the beam""";
    raise TypeError,"theta/phi coordinate interpolation not available with this beam gridder";
    
  def setBeamCube (self,beam):
    """Initializes beam cube. Beam is a 2D (or 3D) cube.""";
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
        for i in range(beam.shape[2]):
          self._beam_real[...,i] = interpolation.spline_filter(beam.real[...,i],order=self._spline_order);
          self._beam_imag[...,i] = interpolation.spline_filter(beam.imag[...,i],order=self._spline_order);
          self._beam_ampl[...,i] = interpolation.spline_filter(self._beam_ampl[...,i],order=self._spline_order);
    
  def transformCoordinates (self,l,m,thetaphi=False,rotate=None,mask=None,freq=None,time=None,freqaxis=None,timeaxis=None,extra_axes=0):
    """Transforms sets of l/m (default) or theta/phi (thetaphi=True), plus freq/time (if supplied) coordinates into an array of 
    local in-beam coordinates. 
    If thetaphi=True, then l,m are actually theta,phi
    If extra_axis>0, then l/m arrays will have that many extra axes (at the front), i.e. the real freq
    axis is actually freqaxis+extra_axis
    An optional 'rotate' argument rotates the coordinates by the given angle (N through W), thus this is really the same thing
    as the local parallactic angle (PA being defined as the angle between local zenith and North, increasing N through E).
    
    mask is an optional array of values to mask
    l,m[,mask][,rotate] may be arrays (both must be the same shape, or will be promoted to the same shape);
    
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
        
    Return value is a tuple of coords,shape,mask where 
    'coords'  is an (N,M) array suitable for feeding to e.g. map_coordinates() (N=2 for an lm-only beam, 3 for an l-m-freq beam, etc.),
    'shape'   is the shape of the output array (where the product of all elements in 'shape' is M).
    'mask'    is a promoted mask array (of shape 'shape'), or None if mask was None
    Note that the 'time' argument is currently ignored -- provided for later compatibility (i.e. beams with time planes)
    """
    # make sure inputs are arrays
    l = numpy.array(l) + self.l0;
    m = numpy.array(m) + self.m0;
    # promote l,m,mask,rotate to the same shape
    l,m,mask,rotate = unite_multiple_shapes(l,m,mask,rotate);
    
    dprint(3,"input l/m shape is",l.shape,m.shape);
    dprint(3,"input l/m[0] is",l.ravel()[0],m.ravel()[0]);
#    dprint(4,"input l/m is",l,m);
    # convert coordinates using one or the other type of converter, depending on lm setting
    l,m = (self.thetaPhiToBeam if thetaphi else self.lmToBeam)(l,m,rotate=rotate);
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
        freq = self.freqToBeam(freq);
        dprint(4,"in beam coordinates this is",freq);
      # case (A): reuse same frequency for every l/m point
      if len(freq) == 1:
        freq = [freq[0]]*l.size
        lm = numpy.vstack((l.ravel(),m.ravel(),freq));
      # case B/C:
      else:
        # first turn freq vector into an array of the proper shape
        if freqaxis is None:
          raise ValueError,"frequency axis not specified, but beam has a frequency dependence";
        freqaxis += extra_axes;
        freqshape = [1]*(freqaxis+1);
        freqshape[freqaxis] = len(freq);
        freq = freq.reshape(freqshape);
        # now promote freq to same shape as l,m. This takes care of cases B and C.
        # (the freq axis of l/m will be expanded, and other axes of freq will be expanded)
        l,m,mask,rotate,freq = unite_multiple_shapes(l,m,mask,rotate,freq);
        freq = freq.ravel()
        lm = numpy.vstack((l.ravel(),m.ravel(),freq));
      # now extrapolate: freq entries<0 mean frequency is out of range, With -1<X<0 meaning
      # frequency is X*lowest, and X<-1 meaning frequency is X*highest. In either case lm needs to be multiplied
      # by X to achieve extrapolation, and the frequency coordinate needs to be set to the first or last.
      wlo = (freq<0)&(freq>=-1)
      whi = (freq<-1)
      lm[:2,wlo|whi] *= -freq[wlo|whi]
      freq[wlo] = 0
      freq[whi] = len(self._freq_grid)-1
    # case (D): no frequency dependence in the beam
    else:
      lm = numpy.vstack((l.ravel(),m.ravel(),[0]*l.size));
    return lm,l.shape,mask;

  def interpolate_3d (self,l,m,thetaphi=False,rotate=None,time=None,freq=None,freqaxis=None,output=None,mask=None,extra_axes=0):
    """Interpolates beam into the given l/m (default) or theta/phi (thetaphi=True) coordinates,
    and optionally freq/time, etc.
    If extra_axis>0, then l/m arrays will have that many extra axes (at the front), i.e. the real freq
    axis is actually freqaxis+extra_axis
    See transformCoordinates() above for a description of the parameters. An output array
    may be specified (it will be resized to the correct shape just in case).
    Returns array of interpolated coordinates.
    """;
    # transform coordinates
    dprint(3,"transforming coordinates");
    coords,output_shape,mask = self.transformCoordinates(l,m,thetaphi=thetaphi,
                               rotate=rotate,time=time,freq=freq,freqaxis=freqaxis,extra_axes=extra_axes,mask=mask);
    # prepare output array
    if output is None:
      output = numpy.zeros(output_shape,complex);
    elif output.shape != output_shape:
      output.resize(output_shape);
    global _print_first;
    verbose = _print_first and _verbosity.get_verbose()>1;
    if verbose:
      l0,m0,freq0 = coords[:,0];
      dprint(1,"First interpolation point is at",l0,m0,freq0);
      dprint(1,"Frequencies are",self._freq_grid);
      for l in int(l0),int(l0)+1:
        for m in int(m0),int(m0)+1:
          dprint(1,l,m,"amplitudes",abs(self._beam[l,m,:]));
          dprint(1,l,m,"phases",numpy.angle(self._beam[l,m,:])/DEG);
    
    dprint(3,"interpolating %s coordinate points to output shape %s"%(coords.shape,output_shape));
    # interpolate real and imag parts separately
    output.real = interpolation.map_coordinates(self._beam_real,coords,order=self._spline_order,mode='nearest',
                    prefilter=(self._spline_order==1)).reshape(output_shape);
    output.imag = interpolation.map_coordinates(self._beam_imag,coords,order=self._spline_order,mode='nearest',
                    prefilter=(self._spline_order==1)).reshape(output_shape);
    output[~(numpy.isfinite(output))] = 0;
    if mask is not None:
      output[mask] = 0;
    output_ampl = interpolation.map_coordinates(self._beam_ampl,coords,order=self._spline_order,mode='nearest',
                    prefilter=(self._spline_order==1)).reshape(output_shape);
    output_ampl[~(numpy.isfinite(output_ampl))] = 0;
    phase_array = numpy.angle(output);
    output.real = output_ampl * numpy.cos(phase_array);
    output.imag = output_ampl * numpy.sin(phase_array);
    
    dprint(3,"interpolated value [0] is",output.ravel()[0]);
    # dprint(4,"interpolated value is",output);
    if verbose:
      _print_first = None;
      for i in range(coords.shape[1]):
        if coords[0,i] == l0 and coords[1,i] == m0:
          dprint(0,"Interpolated amplitude at frequency %f is %f"%(coords[2,i],abs(output.ravel()[i])));
    
    return output;
    
  def _interpolate_freqplane (self,ifreq,coords,verbose=False):
    val = self._freqplanes.get(ifreq);
    if val is None:
      output = numpy.zeros(self._freqplane_shape,complex);
      # interpolate real and imag parts separately
      output.real = interpolation.map_coordinates(self._beam_real[...,ifreq],coords[:2,...],mode='nearest',
                      order=self._spline_order,prefilter=(self._spline_order==1)).reshape(self._freqplane_shape);
      output.imag = interpolation.map_coordinates(self._beam_imag[...,ifreq],coords[:2,...],mode='nearest',
                      order=self._spline_order,prefilter=(self._spline_order==1)).reshape(self._freqplane_shape);
      output_ampl = interpolation.map_coordinates(self._beam_ampl[...,ifreq],coords[:2,...],mode='nearest',
                      order=self._spline_order,prefilter=(self._spline_order==1)).reshape(self._freqplane_shape);
      output[~(numpy.isfinite(output))] = 0;
      output_ampl[~(numpy.isfinite(output_ampl))] = 0;
      if verbose:
        l0,m0,freq = coords[:,0];
        for l in int(l0),int(l0)+1:
          for m in int(m0),int(m0)+1:
            dprint(0,l,m,freq,"beam plane amplitude",abs(self._beam[l,m,ifreq]),self._beam_ampl[l,m,ifreq]);
        dprint(0,"interpolated amplitude",output_ampl[0]);
      phase_array = numpy.angle(output);
      output.real = output_ampl * numpy.cos(phase_array);
      output.imag = output_ampl * numpy.sin(phase_array);
      val = self._freqplanes[ifreq] = output,output_ampl;
    return val;
    

  def interpolate_linfreq (self,l,m,freq,thetaphi=False,rotate=None,time=None,freqaxis=1,output=None,mask=None,extra_axes=0):
    """Interpolates beam into the given l/m (default) or theta/phi (thetaphi=True) coordinates,
    with a fixed frequency grid.
    If extra_axis>0, then l/m arrays will have that many extra axes (at the front), i.e. the real freq
    axis is actually freqaxis+extra_axis
    Interpolates each frequency plane separately, then does a linear interpolation between the 
    frequency planes.
    See transformCoordinates() above for a description of the parameters. An output array
    may be specified (it will be resized to the correct shape just in case).
    Returns array of interpolated coordinates.
    """;
    # transform coordinates
    dprint(3,"transforming coordinates");
    # use freq=None here because we just want the lm coordinates transformed, and we loop over
    # frequencies explicitly below. The resulting output_shape has 1 for the frequency axis.
    coords,output_shape,mask = self.transformCoordinates(l,m,thetaphi=thetaphi,
                                    rotate=rotate,time=time,freq=None,freqaxis=freqaxis,mask=mask);
    freqcoord = self.freqToBeam(freq) if self.hasFrequencyAxis() else 0; 
    if numpy.isscalar(freqcoord) or freqcoord.ndim == 0:
      freqcoord = [float(freqcoord)];
    self._freqplanes = {};
    # now, ensure output shape has the right frequency axis
    output_shape = list(output_shape);
    freqaxis += extra_axes;
    if len(output_shape) <= freqaxis:
      output_shape += [1]*(freqaxis-len(output_shape)+1);
    output_shape[freqaxis] = len(freqcoord);
    self._freqplane_shape = list(output_shape);
    self._freqplane_shape[freqaxis] = 1;
    # prepare output array
    if output is None:
      output = numpy.zeros(output_shape,complex);
    elif output.shape != output_shape:
      output.resize(output_shape);
    global _print_first;
    verbose = _print_first and _verbosity.get_verbose()>1;
    if verbose:
      l0,m0,freq0 = coords[:,0];
      verbose = l0,m0;
      dprint(0,"First interpolation point is at",l0,m0);
      dprint(0,"Frequencies are",self._freq_grid);
      for l in int(l0),int(l0)+1:
        for m in int(m0),int(m0)+1:
          dprint(0,l,m,"amplitudes",abs(self._beam[l,m,:]));
          dprint(0,l,m,"phases",numpy.angle(self._beam[l,m,:])/DEG);
    
    dprint(3,"interpolating %s coordinate points to output shape %s"%(coords.shape,output_shape));
    freqslice = [slice(None)]*len(output_shape);
    reduced_slice = [slice(None)]*len(output_shape);
    reduced_slice[freqaxis] = 0; 
    for ifreq,freq in enumerate(freqcoord):
      freqslice[freqaxis] = ifreq;
      if freq < -1:
        coords *= - freq
        output[freqslice] = self._interpolate_freqplane(len(self._freq_grid)-1,coords)[0][reduced_slice]
      elif freq < 0:
        coords *= - freq
        output[freqslice] = self._interpolate_freqplane(0,coords)[0][reduced_slice]
      elif freq == 0:
        output[freqslice] = self._interpolate_freqplane(0,coords)[0][reduced_slice]
      elif freq >= len(self._freq_grid)-1:
        output[freqslice] = self._interpolate_freqplane(len(self._freq_grid)-1,coords)[0][reduced_slice]
      else:
        f0,f1 = int(freq),int(freq)+1;
        (c0,abs0),(c1,abs1) = self._interpolate_freqplane(f0,coords,
                  verbose and not ifreq),self._interpolate_freqplane(f1,coords,verbose and not ifreq)
        if verbose>1 and not ifreq:
          dprint(0,"per-plane weights",f0,f1);
          dprint(0,"per-plane amplitudes",abs0[0],abs1[0]);
          dprint(0,"per-plane phases",numpy.angle(c0)/DEG,numpy.angle(c1)/DEG);
        ax = abs0*(f1-freq) + abs1*(freq-f0);
        cx = c0*(f1-freq) + c1*(freq-f0);
        px = numpy.angle(cx);
        output[freqslice].real = (ax*numpy.cos(px))[reduced_slice];
        output[freqslice].imag = (ax*numpy.sin(px))[reduced_slice];
      # apply mask, if any
      if mask is not None:
        output[freqslice][mask] = 0;
    
    dprint(3,"interpolated value [0] is",output.ravel()[0]);
    # dprint(4,"interpolated value is",output);
    if verbose:
      _print_first = None;
      dprint(0,"Output shape is",output.shape);
      for ifreq,freq in enumerate(freqcoord):
        dprint(0,"Frequencies",freqcoord);
        dprint(0,"Interpolated amplitudes",abs(output[0,...]));
        dprint(0,"Interpolated phases",numpy.angle(output[0,...])/DEG);
    
    return output;

# set to True to print the first interpolation
_print_first = True;