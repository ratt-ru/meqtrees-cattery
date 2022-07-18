# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division
import os.path
import math
import numpy
from scipy.ndimage import interpolation
from scipy import interpolate

from Cattery import Meow
from Cattery.Meow import Context

import Kittens.utils
_verbosity = Kittens.utils.verbosity(name="vb");
#_verbosity.set_verbose(5)
dprint = _verbosity.dprint;
dprintf = _verbosity.dprintf;

## ugly hack to get around UGLY FSCKING ARROGNAT (misspelling fully intentional) pyfits-2.3 bug
pyfits = Kittens.utils.import_pyfits();


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
    raise TypeError("array must have length 1 along axis %d, it has %d"%(axis,x.shape[axis]));
  return numpy.concatenate([x]*n,axis);

def unite_shapes (a,b):
  """Makes two arrays have the same shape, as follows:
  - each axis must be the same shape,
  - or if one array has shape-1, and the other shape-N, then the shape-1 is expanded.
  """
  if numpy.isscalar(a):
    a = numpy.array([a]);
  if numpy.isscalar(b):
    b = numpy.array([b]);
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
        raise TypeError("error: trying to unite incompatible shapes %s and %s"%(sa,sb));
  return a,b;

class FITSAxes (object):
  """Helper class encapsulating a FITS header."""
  def __init__ (self,hdr):
    """Creates FITSAxes object from FITS header""";
    self._axis = {};
    naxis = hdr['NAXIS'];
    self._naxis = [0]*naxis;
    self._grid = [[]]*naxis;
    self._type = ['']*naxis;
    self._rpix = [0]*naxis;
    self._rval = [0]*naxis;
    self._rval0 = [0]*naxis;
    self._grid = [None]*naxis;
    self._w2p  = [None]*naxis;
    self._p2w  = [None]*naxis;
    self._delta = [1]*naxis;
    self._delta0 = [1]*naxis;
    self._unit = [None]*naxis;
    self._unit_scale = [1.]*naxis;
    # extract per-axis info
    for i in range(naxis):
      ax = str(i+1);
      nx = self._naxis[i] = hdr.get('NAXIS'+ax);
      # CTYPE is axis name
      self._type[i] = ctype = hdr.get('CTYPE'+ax,None);
      if self._type[i] is not None:
        self._axis[self._type[i]] = i;
      # axis gridding
      # use non-standard keywords GRtype1 .. GRtypeN to supply explicit grid values and ignore CRVAL/CRPIX/CDELT
      grid = [ hdr.get('G%s%d'%(ctype,j),None) for j in range(1,nx+1) ];
      if all([x is not None for x in grid]):
        self._grid[i] = numpy.array(grid);
        self._w2p[i] = interpolate.interp1d(grid,list(range(len(grid))),'linear');
        self._p2w[i] = interpolate.interp1d(list(range(len(grid))),grid,'linear');
      else:
        self._rval[i] = self._rval0[i] = rval = hdr.get('CRVAL'+ax,0);
        self._rpix[i] = rpix = hdr.get('CRPIX'+ax,1) - 1;
        self._delta[i] = self._delta0[i] = delta = hdr.get('CDELT'+ax,1);
        self._setup_grid(i);
      self._unit[i] = hdr.get('CUNIT'+ax,'').strip().upper();

  def _setup_grid (self,i):
    """Internal helper to set up the grid based on rval/rpix/delta"""
    nx,rpix,rval,delta = self._naxis[i],self._rpix[i],self._rval[i],self._delta[i];
    self._grid[i] = (numpy.arange(0.,float(nx))-rpix)*delta+rval;
    # warning: small angle approximation
    self._w2p[i] =  lambda world,rpix=rpix,rval=rval,delta=delta:rpix+(world-rval)/delta; # float delta
    self._p2w[i] =  lambda pix,rpix=rpix,rval=rval,delta=delta:(pix-rpix)*delta+rval;

  def ndim (self):
    return len(self._naxis);

  def naxis (self,axis):
    return self._naxis[self.iaxis(axis)];

  def iaxis (self,axisname):
    return axisname if isinstance(axisname,int) else self._axis.get(axisname,-1);

  def grid (self,axis):
    return self._grid[self.iaxis(axis)];

  def type (self,axis):
    return self._type[self.iaxis(axis)];

  def unit (self,axis):
    return self._unit[self.iaxis(axis)];

  def setUnitScale (self,axis,scale):
    iaxis = self.iaxis(axis);
    self._unit_scale[iaxis] = scale;
    self._rval[iaxis] = self._rval0[iaxis]*scale;
    self._delta[iaxis] = self._delta0[iaxis]*scale;
    self._setup_grid(iaxis);

  def toPixel (self,axis,world,sign=1):
    """Converts array of world coordinates to pixel coordinates""";
    return self._w2p[self.iaxis(axis)](sign*world);

  def toWorld (self,axis,pixel,sign=1):
    """Converts array of pixel coordinates to world coordinates""";
    return self._p2w[self.iaxis(axis)](pixel)*sign;


class LMVoltageBeam (object):
  """This class implements a complex voltage beam as a function of LM."""
  def __init__ (self,spline_order=2,l0=0,m0=0,l_axis="L",m_axis="M",
                ampl_interpolation=False,verbose=None):
    self._spline_order = spline_order;
    self.l0, self.m0 = l0, m0;
    # figure out axis names, and whether they should be swapped
    if l_axis[0] == '-':
      self._l_axis_sign,self._l_axis = -1, l_axis[1:]
    else:
      self._l_axis_sign,self._l_axis = 1, l_axis
    if m_axis[0] == '-':
      self._m_axis_sign,self._m_axis = -1, m_axis[1:]
    else:
      self._m_axis_sign,self._m_axis = 1, m_axis
    self.ampl_interpolation = ampl_interpolation
    if verbose:
      _verbosity.set_verbose(verbose);
      _verbosity.enable_timestamps(True);

  def read (self,filename_real,filename_imag=None):
    """Reads beam patterns from FITS files. If only one file is supplied, assumes a real-only beam.
    If two files are supplied, uses them for the real and imaginary parts.
    If 2N files are supplied, treats them as a frequency cube"""
    ff_re = pyfits.open(filename_real)[0];
    # form up complex beam
    beam = numpy.zeros(ff_re.data.shape,complex);
    beam.real = ff_re.data;
    beam_ampl = None
    # add imaginary part
    if filename_imag:
      im_data = pyfits.open(filename_imag)[0].data;
      if im_data.shape != ff_re.data.shape:
        raise TypeError("shape mismatch between FITS files %s and %s"%(filename_real,filename_imag));
      beam.imag = im_data;
      if self.ampl_interpolation:
        beam_ampl = numpy.abs(beam)
    # change order of axis, since FITS has first axis last
    beam = beam.transpose();
    if not beam_ampl is None:
      beam_ampl = beam_ampl.transpose()
    # figure out axes
    self._axes = axes = FITSAxes(ff_re.header);
    # find L/M axes
    laxis = axes.iaxis(self._l_axis);
    maxis = axes.iaxis(self._m_axis);
    if laxis<0 or maxis<0:
      raise TypeError("FITS file %s missing %s or %s axis"%(filename_real,self._l_axis,self._m_axis));
    # setup conversion functions
    self._lToPixel = Kittens.utils.curry(axes.toPixel,laxis,sign=self._l_axis_sign);
    self._mToPixel = Kittens.utils.curry(axes.toPixel,maxis,sign=self._m_axis_sign);
    # find frequency grid. self._freqToPixel will be None if no frequency axis
    freqaxis = axes.iaxis('FREQ');
    if freqaxis >= 0 and axes.naxis(freqaxis) > 1:
      dprint(1,"FREQ axis has %d points"%axes.naxis(freqaxis));
      self._freqgrid = axes.grid(freqaxis);
      self._freqToPixel = Kittens.utils.curry(axes.toPixel,freqaxis);
      used_axes = [laxis,maxis,freqaxis];
    else:
      self._freqToPixel = None;
      used_axes = [laxis,maxis];
    # used_axes is either (l,m,freq) or (l,m)
    # other_axes is all that remains, and they had better be all trivial
    other_axes = sorted(set(range(axes.ndim())) - set(used_axes));
    if any([axes.naxis(i)>1 for i in other_axes]):
      raise TypeError("FITS file %s has other non-trivial axes besides L/M"%filename_real);
    # setup units
    for ax in laxis,maxis:
      dprint(1,"%s axis unit is %s"%(axes.type(ax),axes.unit(ax)));
      if not axes.unit(ax) or axes.unit(ax).upper() == "DEG":
        axes.setUnitScale(ax,DEG);
    # transpose array into L,M order and reshape
    dprint(1,"beam array has shape",beam.shape);
    beam = beam.transpose(used_axes+other_axes);
    beam = beam.reshape(beam.shape[:len(used_axes)]);
    if not beam_ampl is None:
      beam_ampl = beam_ampl.transpose(used_axes+other_axes)
      beam_ampl = beam_ampl.reshape(beam_ampl.shape[:len(used_axes)]);

    dprint(1,"beam array has shape",beam.shape);
    dprint(2,"l grid is",axes.grid(laxis));
    dprint(2,"m grid is",axes.grid(maxis));
    if self._freqToPixel:
      dprint(2,"freq grid is",axes.grid(freqaxis));
    # prefilter beam for interpolator
    self._beam = beam;
    if self._spline_order > 1:
      self._beam_real = interpolation.spline_filter(beam.real,order=self._spline_order);
      self._beam_imag = interpolation.spline_filter(beam.imag,order=self._spline_order);
      if not beam_ampl is None:
        self._beam_ampl = interpolation.spline_filter(beam_ampl,order=self._spline_order);
      else:
        self._beam_ampl = beam_ampl
    else:
      self._beam_real = beam.real;
      self._beam_imag = beam.imag;
      self._beam_ampl = beam_ampl

  def hasFrequencyAxis (self):
    return bool(self._freqToPixel);

  def beam (self):
    return self._beam;

  def interpolate (self,l,m,time=None,freq=None,freqaxis=None,output=None):
    """Interpolates l/m coordinates in the beam.
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

        In all three cases, if any points in freq are outside the beam's frequency axis span,
        then the first/last channel in the beam will be scaled up/down with wavelength.

    And finally (D):

    (D) No dependence on frequency in the beam.
        We simply interpolate every l/m value as is. Output array is same shape as l/m.

    'time' is currently ignored -- provided for later compatibility (i.e. beams with time planes)
    """
    # make sure inputs are arrays
    l = numpy.array(l) + self.l0;
    m = numpy.array(m) + self.m0;
    freq = numpy.array(freq);
    # promote l,m to the same shape
    l,m = unite_shapes(l,m);
    dprint(3,"input l/m [0] is",l.ravel()[0],m.ravel()[0],"and shapes are",l.shape);
    # now we make a 2xN coordinate array for map_coordinates
    # lm[0,:] will be flattened L array, lm[1,:] will be flattened M array
    # lm[2,:] will be flattened freq array (if we have a freq dependence)
    # Do we have a frequency axis in the beam? (case A,B,C):
    if self.hasFrequencyAxis():
      if freq is None:
        raise ValueError("frequencies not specified, but beam has a frequency dependence");
      freq = numpy.array(freq);
      if not freq.ndim:
        freq = freq.reshape(1);
      dprint(3,"frequencies are",freq)
      ## catch out-of-bounds frequencies
      below = freq < self._freqgrid[0]
      above = freq > self._freqgrid[-1]
      # lm scale factor is 1 for in-bounds channels, <1 for below, >1 for above.
      scale = numpy.ones(len(freq),float)
      scale[below] = freq[below]/self._freqgrid[0]
      scale[above] = freq[above]/self._freqgrid[-1]
      freq[below] = self._freqgrid[0]
      freq[above] = self._freqgrid[-1]
      # convert frequency to fractional channel index
      chan = self._freqToPixel(freq)
      dprint(3,"in frequency plane coordinates we have",chan)
      # case (A): reuse same frequency for every l/m point
      if len(chan) == 1:
        chan = numpy.array([chan[0]]*l.size)
        scale = scale[0]  # single lm scale factor 
        lm = numpy.vstack((l.ravel(),m.ravel(),chan));
      # case B/C:
      else:
        # first turn chan vector into an array of the proper shape
        if freqaxis is None:
          raise ValueError("frequency axis not specified, but beam has a frequency dependence");
        freqshape = [1]*(freqaxis+1);
        freqshape[freqaxis] = len(chan);
        chan = chan.reshape(freqshape)
        scale = scale.reshape(freqshape)
        # now promote chan to same shape as l,m. This takes care of cases B and C.
        # (the chan axis of l/m will be expanded, and other axes of chan will be expanded)
        l1,chan = unite_shapes(l,chan)
        m1,chan = unite_shapes(m,chan)
        l,scale = unite_shapes(l,scale)
        m,scale = unite_shapes(m,scale)
        lm = numpy.vstack((l.ravel(),m.ravel(),chan.ravel()))
        scale = scale.ravel()
      if above.any() or below.any():
        lm[0,:] *= scale
        lm[1,:] *= scale
        dprint(3,"some points were extrapolated for OOB frequencies using scale factors")
    # case (D): no frequency dependence in the beam
    else:
      lm = numpy.vstack((l.ravel(),m.ravel()));
    # now convert lm to pixel coordinates
    lm[0,:] = self._lToPixel(lm[0,:])
    lm[1,:] = self._mToPixel(lm[1,:])
    dprint(3,"xy pixel coordinates are [0]",lm[0,0],lm[1,0]);
    # interpolate and reshape back to shape of L
    if output is None:
      output = numpy.zeros(l.shape,complex);
    elif output.shape != l.shape:
      output.resize(l.shape);
    dprint(3,"interpolating %d lm points"%lm.shape[1]);
    output.real = interpolation.map_coordinates(self._beam_real,lm,order=self._spline_order,
                  prefilter=(self._spline_order==1),mode='nearest').reshape(l.shape);
    output.imag = interpolation.map_coordinates(self._beam_imag,lm,order=self._spline_order,
                  prefilter=(self._spline_order==1),mode='nearest').reshape(l.shape);
    if not self._beam_ampl is None:
      output_ampl = interpolation.map_coordinates(self._beam_ampl,lm,order=self._spline_order,
                  prefilter=(self._spline_order==1),mode='nearest').reshape(l.shape);
      phase_array = numpy.arctan2(output.imag,output.real)
      output.real = output_ampl * numpy.cos(phase_array)
      output.imag = output_ampl * numpy.sin(phase_array)
    dprint(3,"interpolated value [0] is",output.ravel()[0]);
    dprint(4,"interpolated value is",output);
    return output;

class LMVoltageMultifreqBeam (LMVoltageBeam):
  """This class implements an LMVoltageBeam where the
  different frequency planes are read from different FITS files."""
  def read (self,filenames):
    freqs = [];
    for ifreq,(filename_real,filename_imag) in enumerate(filenames):
      """Reads beam patterns from FITS files. If only one file is supplied, assumes a real-only beam.
      If two files are supplied, uses them for the real and imaginary parts.
      If 2N files are supplied, treats them as a frequency cube"""
      ff_re = pyfits.open(filename_real)[0];
      # form up complex beam
      beam = numpy.zeros(ff_re.data.shape,complex);
      beam.real = ff_re.data;
      # add imaginary part
      if filename_imag:
        im_data = pyfits.open(filename_imag)[0].data;
        if im_data.shape != ff_re.data.shape:
          raise TypeError("shape mismatch between FITS files %s and %s"%(filename_real,filename_imag));
        beam.imag = im_data;
      # change order of axis, since FITS has first axis last
      beam = beam.transpose();
      # figure out axes
      axes = FITSAxes(ff_re.header);
      used_axes = [ axes.iaxis(x) for x in ("L","M","FREQ") ];
      if any([x<0 for x in used_axes]):
        raise TypeError("FITS file %s missing L, M or FREQ axis");
      laxis,maxis,freqaxis = used_axes;
      # check the other axes
      other_axes = sorted(set(range(axes.ndim())) - set(used_axes));
      if any([axes.naxis(i)>1 for i in other_axes]):
        raise TypeError("FITS file %s has other non-trivial axes besides L/M"%filename_real);
      # setup frequency grid
      freqgrid = axes.grid(freqaxis);
      if len(freqgrid) > 1:
        raise TypeError("FITS file %s has >1 frequency points");
      if freqs and freqgrid[0] < freqs[-1]:
        raise TypeError("FITS file %s has lower frequency than previous file -- monotonically increasing frequencies are expected");
      freqs.append(freqgrid[0]);
      # check if it matches previous image
      if not ifreq:
        baseshape = beam.shape;
        self._axes = axes;
        # setup 3D beam cube
        beamcube = numpy.zeros((axes.naxis(laxis),axes.naxis(maxis),len(filenames)),complex);
        # setup conversion functions
        self._lToPixel = Kittens.utils.curry(axes.toPixel,laxis);
        self._mToPixel = Kittens.utils.curry(axes.toPixel,maxis);
        for ax in laxis,maxis:
          dprint(1,"%s axis unit is %s"%(self._axes.type(ax),self._axes.unit(ax)));
          if not self._axes.unit(ax) or self._axes.unit(ax).upper() == "DEG":
            self._axes.setUnitScale(ax,DEG);
      else:
        if baseshape != beam.shape:
          raise TypeError("FITS file %s has differing dimensions"%filename_real);
      beam = beam.transpose(used_axes+other_axes);
      beam = beam.reshape(beam.shape[:len(used_axes)]);
      beamcube[:,:,ifreq] = beam[:,:,0];
    # done reading the beam cube
    dprint(1,"beam array has shape",beamcube.shape);
    dprint(2,"l grid is",self._axes.grid(laxis));
    dprint(2,"m grid is",self._axes.grid(maxis));
    dprint(2,"freq grid is",freqs);
    self._freqaxis = freqs;
    self._freq_interpolator = interpolate.interp1d(freqs,list(range(len(freqs))),'linear');
    # prefilter beam for interpolator
    self._beam = beamcube;
    self._beam_ampl = numpy.abs(beamcube) if self.ampl_interpolation else None;
    if self._spline_order > 1:
      self._beam_real = interpolation.spline_filter(beamcube.real,order=self._spline_order);
      self._beam_imag = interpolation.spline_filter(beamcube.imag,order=self._spline_order);
      if self._beam_ampl:
        self._beam_ampl = interpolation.spline_filter(self._beam_ampl,order=self._spline_order);
    else:
      self._beam_real = beamcube.real;
      self._beam_imag = beamcube.imag;

  def hasFrequencyAxis (self):
    return True;

  def _freqToPixel (self,freq):
    return self._freq_interpolator(freq);

try:
  from Timba import pynode
  from Timba.Meq import meq
  from Timba import mequtils
  standalone = False
except:
  print("importing InterpolatedBeams in standalone mode, ignoring MeqTrees-related stuff")
  standalone = True

if not standalone:  
  def _cells_grid (obj,axis):
    """helper function to get a grid out of the cells object. Returns None if none is found"""
    if hasattr(obj,'cells') and hasattr(obj.cells.grid,axis):
      return obj.cells.grid[axis];
    else:
      return None;

  class FITSBeamInterpolatorNode (pynode.PyNode):
    def __init__ (self,*args):
      pynode.PyNode.__init__(self,*args);
      # Maintain a global dict of VoltageBeam objects per each filename set, so that we reuse them
      # We also reset this dict each time a node is created (otherwise the objects end up being reused
      # even after the tree has been rebuilt.)
      global _voltage_beams;
      _voltage_beams = {};
      self._vb_key = None # leave until we select the closest beam in frequency

    def update_state (self,mystate):
        """Standard function to update our state""";
        mystate('filename_real',[]);
        mystate('filename_imag',[]);
        mystate('spline_order',3);
        mystate('normalize',False);
        mystate('ampl_interpolation',False);
        mystate('verbose',0);
        mystate('l_axis',"L")
        mystate('m_axis',"M")
        mystate('l_beam_offset',0.0);
        mystate('m_beam_offset',0.0);
        mystate('missing_is_null',False);
        # Check filename arguments, and init _vb_key for init_voltage_beams() below
        self._vb_keys = []
        def __check_make_list(inpset):
          if isinstance(inpset, str):
            return [inpset]
          elif (isinstance(inpset, tuple) or isinstance(inpset, list)):
            if all(map(lambda x: isinstance(x, str), inpset)): # singular beamset
              return [inpset]
            else:
              return inpset # treat as multiple beamsets
        # keep the interface the same as before -- may pass in a singular beamset
        # or a list of beamsets with different frequency coverages
        # Each beamset may be a single filename pair (scalar Jones term), or 4 filenames (full 2x2 matrix)
        self.filename_real = __check_make_list(self.filename_real)
        self.filename_imag = __check_make_list(self.filename_imag)
        for (ifilename_real, ifilename_imag) in zip(self.filename_real, self.filename_imag):
          if isinstance(ifilename_real,str) and isinstance(ifilename_imag,str):
            self._vb_keys.append(((ifilename_real, ifilename_imag),))
          elif len(ifilename_real) == 4 and len(ifilename_imag) == 4:          
            self._vb_keys.append(tuple(zip(ifilename_real, ifilename_imag)))
          else:
            raise ValueError("filename_real/filename_imag: either a single filename, "
                            "or lists of 4 re and 4 im filenames expected per pattern")
        # other init
        mequtils.add_axis('l')
        mequtils.add_axis('m')
        self._freqaxis = mequtils.get_axis_number("freq")
        _verbosity.set_verbose(self.verbose);

    def init_voltage_beams (self, grid):
        """initializes VoltageBeams for the given set of FITS files (per each _vb_key, that is).
        Returns list of 1 or 4 VoltageBeam objects."""
        obsfreq0 = grid['freq'].min()
        obsfreq1 = grid['freq'].max()
        # maintain a global dict of VoltageBeam objects per each filename set, so that we reuse them
        global _voltage_beams;
        if not '_voltage_beams' in globals():
          _voltage_beams = {};
        # get VoltageBeam object from global dict, or init new one if not already defined
        for vb_key in self._vb_keys:
          vbs, beam_max, freq_distance = \
            _voltage_beams.get(vb_key, (None, None, None))
          if vbs:
            for filename_real, filename_imag in vb_key:
              dprint(1, "beam files", filename_real, filename_imag, "already in memory, reusing")
          else:
            vbs = []
            # for each correlation if multi-correlation beams
            for filename_real, filename_imag in vb_key:
              # if files do not exist, replace with blanks
              dprint(1, "loading beam files", filename_real, filename_imag)
              if filename_real and not os.path.exists(filename_real) and self.missing_is_null:
                dprint(1,"beam pattern",filename_real, "not found, using null instead")
                filename_real = None
              if filename_imag and not os.path.exists(filename_imag) and self.missing_is_null:
                dprint(1,"beam pattern",filename_imag, "not found, using null instead")
                filename_real = None
              # now, create VoltageBeam if at least the real part still exists
              if filename_real:
                vb = LMVoltageBeam(
                      l0=self.l_beam_offset,m0=self.m_beam_offset,
                      l_axis=self.l_axis,m_axis=self.m_axis,
                      ampl_interpolation=self.ampl_interpolation,spline_order=self.spline_order,
                      verbose=self.verbose)
                vb.read(filename_real, filename_imag)
              else:
                vb = None
              vbs.append(vb)

            # work out norm of beam
            if not any(vbs):
              raise RuntimeError("no beam patterns have been loaded. Please check your filename pattern")
            if len(vbs) == 1:
              beam_max = abs(vbs[0].beam()).max()
            elif len(vbs) == 4:
              xx,xy,yx,yy = [ vb.beam() if vb else 0 for vb in vbs ]
              beam_max = math.sqrt((abs(xx)**2+abs(xy)**2+abs(yx)**2+abs(yy)**2).max()/2)
            else:
              raise RuntimeError("Either 1 (scalar correlation) or 4 beams (one per correlation) should be loaded")
            dprint(1,"beam max is",beam_max)
            
            # work out the distance from the Observation frequency range
            # assume these will have matching frequency ranges
            freq_distance = max(vbs[0]._freqgrid[0] - obsfreq0, 0) + \
                            max(obsfreq1 - vbs[0]._freqgrid[-1], 0)

            _voltage_beams.setdefault(vb_key, (vbs, beam_max, freq_distance))
        
        # now select the closest from the loaded vb set from self._vb_keys
        sel_key = sorted(filter(lambda x: x in self._vb_keys, _voltage_beams),
                         key=lambda vb_key: _voltage_beams[vb_key][2])[0]
        sel_vbs, sel_beam_max, sel_freq_distance = _voltage_beams[sel_key]
        if len(self._vb_keys) > 1:
          printfn = "{%s}" % ",".join(["('%s', '%s')" % corrset for corrset in sel_key])
          if sel_freq_distance == 0:
              dprint(1, "beam patterns %s overlap the frequency coverage" % printfn)
          else:
              dprint(1, "beam patterns %s are closest to the frequency coverage (%.1f MHz max separation)" % (
                              printfn, sel_freq_distance*1e-6))
          dprint(1,"MS coverage is %.1f to %.1f GHz, beams are %.1f to %.1f MHz"%(
              obsfreq0*1e-6, obsfreq1*1e-6, sel_vbs[0]._freqgrid[0]*1e-6, sel_vbs[0]._freqgrid[-1]*1e-6))
        return sel_vbs, sel_beam_max

    def get_result (self,request,*children):
      # now, figure out the lm and time/freq grid
      # lm may be a 2/3-vector or an Nx2/3 tensor
      lm = children[0];
      dims = getattr(lm,'dims',[len(lm.vellsets)]);
      if len(dims) == 2 and dims[1] in (2,3):
        nsrc,nlm = dims;
        tensor = True;
      elif len(dims) == 1 and dims[0] in (2,3):
        nsrc,nlm = 1,dims[0];
        tensor = False;
      else:
        print("child 0: %d vellsets, shape %s"%(len(lm.vellsets),getattr(lm,'dims',[])));
        raise TypeError("expecting a 2/3-vector or an Nx2/3 matrix for child 0 (lm)");
      # pointing offsets (child 1) are optional
      if len(children) > 1:
        dlm = children[1];
        if len(dlm.vellsets) != 2:
          raise TypeError("expecting a 2-vector for child 1 (dlm)");
        dl,dm = dlm.vellsets[0].value,dlm.vellsets[1].value;
      else:
        dl = dm = 0;
      # setup grid dict that will be passed to VoltageBeam.interpolate
      grid = dict();
      for axis in 'time','freq':
        values = _cells_grid(lm,axis);
        if values is None:
          values = _cells_grid(request,axis);
        if values is not None:
          grid[axis] = values;
      vellsets = [];
      # init voltagebeam(s) from beamset(s)
      vbs,beam_max = self.init_voltage_beams(grid)
      # now loop over sources and interpolte
      for isrc in range(nsrc):
        # put l,m into grid
        l,m = lm.vellsets[isrc*nlm].value,lm.vellsets[isrc*nlm+1].value;
        l,dl = unite_shapes(l,dl);
        m,dm = unite_shapes(m,dm);
        grid['l'] = l - dl;
        grid['m'] = m - dm;
        # interpolate
        for vb in vbs:
          if vb is None:
            vellsets.append(meq.vellset(meq.sca_vells(0.)));
          else:
            beam = vb.interpolate(freqaxis=self._freqaxis,**grid);
            if self.normalize and beam_max != 0:
              beam /= beam_max;
            vells = meq.complex_vells(beam.shape);
            vells[...] = beam[...];
            # make vells and return result
            vellsets.append(meq.vellset(vells));
      # create result object
      cells = request.cells if vb.hasFrequencyAxis() else getattr(lm,'cells',None);
      result = meq.result(vellsets[0],cells=cells);
      if len(vellsets) > 1:
        result.vellsets[1:] = vellsets[1:];
      # vbs is either length 4 or length 1. If length 4, then result needs to have its dimensions
      if len(vbs) > 1:
        result.dims = (nsrc,2,2) if tensor else (2,2);
      return result;



# test clause
if __name__ == "__main__":
  import sys
  _verbosity.set_verbose(5);


  vb = LMVoltageMultifreqBeam(spline_order=3);
  vb.read(
    [("mk_1455MHz_feed_5deg_xxreal.fits","mk_1455MHz_feed_5deg_xximag.fits"),
     ("mk_1460MHz_feed_5deg_xxreal.fits","mk_1460MHz_feed_5deg_xximag.fits")]);

  l0 = numpy.array([-2,0,2])*DEG;
  l = numpy.vstack([l0]*len(l0));

  a = vb.interpolate(l,l.T,freq=[1.456e+9],freqaxis=2);
  b = vb.interpolate(l,l.T,freq=[1.456e+9,1.457e+9,1.458e+9],freqaxis=2);
  c = vb.interpolate(l,l.T,freq=[1.455e+9,1.457e+9,1.458e+9,1.46e+9],freqaxis=2);

  print("C",c.shape,c);
  print("B",b.shape,b);
  print("A",a.shape,a);
  sys.exit(1);


  vb = LMVoltageBeam(spline_order=3);
  vb.read("beam_xx_re.fits","beam_xx_im.fits");

  l0 = numpy.array([-2,-1,0,1,2])*DEG;
  l = numpy.vstack([l0]*len(l0));

  print(vb.interpolate(l,l.T));

  vb = LMVoltageBeam(spline_order=3);
  vb.read("XX_0_Re.fits","XX_0_Im.fits");

  l0 = numpy.array([-4,-2,0,2,4])*DEG;
  l = numpy.vstack([l0]*len(l0));

  a = vb.interpolate(l,l.T,freq=[1e+9],freqaxis=2);
  b = vb.interpolate(l,l.T,freq=[1e+9,1.1e+9,1.2e+9],freqaxis=2);
  c = vb.interpolate(l,l.T,freq=[1e+9,1.1e+9,1.2e+9,1.3e+9,1.4e+9],freqaxis=1);

  print("A",a.shape,a);
  print("B",b.shape,b);
  print("C",c.shape,c);

