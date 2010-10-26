# -*- coding: utf-8 -*-

import os.path
import math
import numpy
import pyfits
from scipy.ndimage import interpolation

import Kittens.utils
_verbosity = Kittens.utils.verbosity(name="vb");
dprint = _verbosity.dprint;
dprintf = _verbosity.dprintf;

DEG = math.pi/180;

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
    self._delta = [1]*naxis;
    self._delta0 = [1]*naxis;
    self._unit = [None]*naxis;
    self._unit_scale = [1.]*naxis;
    # extract per-axis info
    for i in range(naxis):
      ax = str(i+1);
      nx = self._naxis[i] = hdr.get('NAXIS'+ax);
      # CTYPE is axis name
      self._type[i] = hdr.get('CTYPE'+ax,None);  
      if self._type[i] is not None:
        self._axis[self._type[i]] = i;
      # axis gridding
      self._rval[i] = rval = hdr.get('CRVAL'+ax,0);
      self._rpix[i] = rpix = hdr.get('CRPIX'+ax,1) - 1;
      self._delta[i] = self._delta0[i] = delta = hdr.get('CDELT'+ax,1);
      self._unit[i] = hdr.get('CUNIT'+ax,'').strip().upper();

  def ndim (self):
    return len(self._naxis);

  def naxis (self,axis):
    return self._naxis[self.iaxis(axis)];

  def iaxis (self,axisname):
    return axisname if isinstance(axisname,int) else self._axis.get(axisname,-1);

  def grid (self,axis):
    iaxis = self.iaxis(axis);
    return (numpy.arange(0.,float(self._naxis[iaxis])) - self._rpix[naxis])*self._delta[naxis] + self._rval[naxis];

  def type (self,axis):
    return self._type[self.iaxis(axis)];

  def unit (self,axis):
    return self._unit[self.iaxis(axis)];

  def setUnitScale (self,axis,scale):
    iaxis = self.iaxis(axis);
    self._unit_scale[iaxis] = scale;
    self._delta[iaxis] = self._delta0[iaxis]*scale;

  def toPixel (self,axis,world):
    """Converts array of world coordinates to pixel coordinates""";
    iaxis = self.iaxis(axis);
    return self._rpix[iaxis] + (world - self._rval[iaxis])/self._delta[iaxis];

  def toWorld (self,axis,pixel):
    """Converts array of pixel coordinates to world coordinates""";
    iaxis = self.iaxis(axis);
    return (pixel - self._rpix[iaxis])*self._delta[iaxis] + self._rval[iaxis];

  def grid (self,axis):
    iaxis = self.iaxis(axis);
    return self.toWorld(iaxis,numpy.arange(0.,float(self._naxis[iaxis])));

class LMVoltageBeam (object):
  """This class implements a complex voltage beam as a function of LM."""
  def __init__ (self,spline_order=2,verbose=None):
    self._spline_order = spline_order;
    if verbose:
      _verbosity.set_verbose(verbose);

  def read (self,filename_real,filename_imag=None):
    """Reads beam patterns from FITS files. If only one file is supplied, assumes a real-only beam.
    If two files are supplied, uses them for the real and imaginary parts."""
    ff_re = pyfits.open(filename_real)[0];
    # form up complex beam
    beam = numpy.zeros(ff_re.data.shape,complex);
    beam.real = ff_re.data;
    # add imaginary part
    if filename_imag:
      im_data = pyfits.open(filename_imag)[0].data;
      if im_data.shape != ff_re.data.shape:
        raise TypeError,"shape mismatch between FITS files %s and %s"%(filename_real,filename_imag);
      beam.imag = im_data;
    # change order of axis, since FITS has first axis last
    beam = beam.transpose();
    # figure out axes
    self._axes = axes = FITSAxes(ff_re.header);
    laxis = axes.iaxis('L');
    maxis = axes.iaxis('M');
    if laxis<0 or maxis<0:
      raise TypeError,"FITS file %s missing L or M axis"%filename_real;
    other_axes = sorted(set(range(axes.ndim())) - set([laxis,maxis]));
    if any([axes.naxis(i)>1 for i in other_axes]):
      raise TypeError,"FITS file %s has other non-trivial axes besides L/M"%filename_real;
    # setup units
    for ax in laxis,maxis:
      dprint(1,"%s axis unit is %s"%(axes.type(ax),axes.unit(ax)));
      if not axes.unit(ax) or axes.unit(ax).upper() == "DEG":
        axes.setUnitScale(ax,DEG);
    # transpose array into L,M order and reshape
    dprint(1,"beam array has shape",beam.shape);
    beam = beam.transpose([laxis,maxis]+other_axes);
    beam = beam.reshape(beam.shape[:2]);
    dprint(1,"beam array has shape",beam.shape);
    dprint(2,"l grid is",axes.grid(laxis));
    dprint(2,"m grid is",axes.grid(maxis));
    # setup conversion functions
    self._lToPixel = Kittens.utils.curry(axes.toPixel,laxis);
    self._mToPixel = Kittens.utils.curry(axes.toPixel,maxis);
    # prefilter beam for interpolator
    if self._spline_order > 1:
      self._beam_real = interpolation.spline_filter(beam.real,order=self._spline_order);
      self._beam_imag = interpolation.spline_filter(beam.imag,order=self._spline_order);
    else:
      self._beam_real = beam.real;
      self._beam_imag = beam.imag;
 
  def interpolate (self,l,m,time=None,freq=None,output=None):
    """Interpolates l/m coordinates in the beam.
    l,m may be arrays (both must be the same shape).
    Returns interpolated beam values. Return array will have the same shape as the l/m arrays.
    time/freq are ignored -- provided for later compatibility (i.e. beams with time/freq planes)
    """
    # check inputs
    if isinstance(l,(float,int)):
      l = numpy.array([l]);
    if isinstance(m,(float,int)):
      m = numpy.array([m]);
    if l.shape != m.shape:
      raise ValueError,"shapes of l/m do not match: %s and %s"%(l.shape,m.shape);
    dprint(3,"input l/m is",l,m);
    l = self._lToPixel(l);
    m = self._mToPixel(m);
    dprint(3,"in pixel coordinates this is",l,m);
    # make coordinate array for map_coordinates
    # lm[0,:] will be flattened L array, lm[1,:] will be flattened M array
    lm = numpy.vstack((l.ravel(),m.ravel()));
    # interpolate and reshape back to shape of L
    if output is None:
      output = numpy.zeros(l.shape,complex);
    output.real = interpolation.map_coordinates(self._beam_real,lm,order=self._spline_order,
                  prefilter=(self._spline_order==1)).reshape(l.shape);
    output.imag = interpolation.map_coordinates(self._beam_imag,lm,order=self._spline_order,
                  prefilter=(self._spline_order==1)).reshape(l.shape);
    dprint(3,"interpolated value is",output);
    return output;


from Timba import pynode
from Timba.Meq import meq
from Timba import mequtils

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

  def update_state (self,mystate):
    """Standard function to update our state""";
    mystate('filename_real',[]);
    mystate('filename_imag',[]);
    mystate('spline_order',3);
    mystate('verbose',0);
    mystate('missing_is_null',False);
    # Check filename arguments, and init _vb_key for init_voltage_beams() below
    # We may be created with a single filename pair (scalar Jones term), or 4 filenames (full 2x2 matrix)
    if isinstance(self.filename_real,str) and isinstance(self.filename_imag,str):
      self._vb_key = ((self.filename_real,self.filename_imag),);
    elif  len(self.filename_real) == 4 and len(self.filename_imag) == 4:
      self._vb_key = tuple(zip(self.filename_real,self.filename_imag));
    else:
      raise ValueError,"filename_real/filename_imag: either a single filename, or a list of 4 filenames expected";
    # other init
    mequtils.add_axis('l');
    mequtils.add_axis('m');
    _verbosity.set_verbose(self.verbose);

  def init_voltage_beams (self):
    """initializes VoltageBeams for the given set of FITS files (per each _vb_key, that is). 
    Returns list of 1 or 4 VoltageBeam objects."""
    # maintain a global dict of VoltageBeam objects per each filename set, so that we reuse them
    global _voltage_beams;
    if not '_voltage_beams' in globals():
      _voltage_beams = {};
    # get VoltageBeam object from global dict, or init new one if not already defined
    vbs = _voltage_beams.get(self._vb_key);
    if not vbs:
      vbs = [];
      for filename_real,filename_imag in self._vb_key:
        # if files do not exist, replace with blanks
        if not os.path.exists(filename_real) and self.missing_is_null:
          filename_real = None;
        if not os.path.exists(filename_imag) and self.missing_is_null:
          filename_real = None;
        # now, create VoltageBeam if at least the real part still exists
        if filename_real:
          vb = LMVoltageBeam(spline_order=self.spline_order,verbose=self.verbose);
          vb.read(filename_real,filename_imag);
        else:
          vb = None;
        vbs.append(vb);
      _voltage_beams[self._vb_key] = vbs;
    return vbs; 

  def get_result (self,request,*children):
    # get list of VoltageBeams
    vbs = self.init_voltage_beams();
    # now, figure out the lm and time/freq grid
    lm = children[0];
    l = lm.vellsets[0].value;
    m = lm.vellsets[1].value;
    # setup grid dict that will be passed to VoltageBeam.interpolate
    grid = dict(l=l,m=m);
    for axis in 'time','freq':
      values = _cells_grid(lm,axis);
      if values is None:
        values = _cells_grid(request,axis);
      if values is not None:
        grid[axis] = values;
    # interpolate
    vellsets = [];
    for vb in vbs:
      if vb is None:
        vellsets.append(meq.vellset(meq.sca_vells(0.)));
      else:
        vells = meq.complex_vells(l.shape);
        beam = vb.interpolate(output=vells,**grid);
        # make vells and return result
        vellsets.append(meq.vellset(vells));
    # create result object
    result = meq.result(vellsets[0],cells=getattr(lm,'cells',None));
    # if more than one vellset, then we have 2x2
    if len(vellsets) > 1:
      result.vellsets[1:] = vellsets[1:];
      result.dims = (2,2);
    return result;


# test clause
if __name__ == "__main__":
  _verbosity.set_verbose(5);

  vb = LMVoltageBeam(spline_order=3);
  vb.read("beam_xx_re.fits","beam_xx_im.fits");
    
  l0 = numpy.array([-2,-1,0,1,2])*DEG;
  l = numpy.vstack([l0]*len(l0));

  print vb.interpolate(l,l.T);
  

