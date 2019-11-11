# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division
import os.path
import math
import numpy
from scipy.ndimage import interpolation

import Kittens.utils
## ugly hack to get around UGLY FSCKING ARROGNAT (misspelling fully intentional) pyfits-2.3 bug
pyfits = Kittens.utils.import_pyfits();

import Siamese.OMS.InterpolatedBeams

from Siamese.OMS.InterpolatedBeams import _verbosity,dprint,dprintf,DEG,LMVoltageBeam

from Timba import pynode
from Timba.Meq import meq
from Timba import mequtils

def _cells_grid (obj,axis):
  """helper function to get a grid out of the cells object. Returns None if none is found"""
  if hasattr(obj,'cells') and hasattr(obj.cells.grid,axis):
    return obj.cells.grid[axis];
  else:
    return None;

class FITSCompoundBeamInterpolatorNode (pynode.PyNode):
  """This reads a list of 2N complex voltage beams (as real/imaginary parts, all Xs then all Ys), and interpolates the lm coordinates of the
  children through them. The result is [2,N] vellsets.
  """
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
    mystate('normalize',False);
    mystate('ampl_interpolation',False);
    mystate('l_0',0.0);
    mystate('m_0',0.0);
    mystate('verbose',0);
    mystate('missing_is_null',False);
    # Check filename arguments: we must be created with two identical-length lists
    if isinstance(self.filename_real,(list,tuple)) and isinstance(self.filename_imag,(list,tuple)) \
        and len(self.filename_real) == len(self.filename_imag) and not len(self.filename_real)%1:
      self._vb_key = tuple(zip(self.filename_real,self.filename_imag));
    else:
      raise ValueError("filename_real/filename_imag: two lists of filenames of 2N elements each expected");
    # other init
    mequtils.add_axis('l');
    mequtils.add_axis('m');
    self._freqaxis = mequtils.get_axis_number("freq");
    _verbosity.set_verbose(self.verbose);

  def init_voltage_beams (self):
    """initializes VoltageBeams for the given set of FITS files (per each _vb_key, that is).
    Returns list of 1 or 4 VoltageBeam objects."""
    # maintain a global dict of VoltageBeam objects per each filename set, so that we reuse them
    global _voltage_beams;
    if not '_voltage_beams' in globals():
      _voltage_beams = {};
    # get VoltageBeam object from global dict, or init new one if not already defined
    vbs,beam_max = _voltage_beams.get(self._vb_key,(None,None));
    if not vbs:
      vbs = [];
      for filename_real,filename_imag in self._vb_key:
        # if files do not exist, replace with blanks
        if not ( os.path.exists(filename_real) and os.path.exists(filename_imag) ) and self.missing_is_null:
          filename_real = None;
          print("No beam pattern %s or %s, assuming null beam"%(filename_real,filename_imag));
        # now, create VoltageBeam if at least the real part still exists
        if filename_real:
          vb = LMVoltageBeam(
                  l0=self.l_0,m0=self.m_0,
                  ampl_interpolation=self.ampl_interpolation,spline_order=self.spline_order,
                  verbose=self.verbose);
          vb.read(filename_real,filename_imag);
        else:
          vb = None;
        # work out norm of beam
        vbs.append(vb);
      xx = [ vb.beam() if vb else numpy.array([0]) for vb in vbs[:len(vbs)//2] ];
      yy = [ vb.beam() if vb else numpy.array([0]) for vb in vbs[len(vbs)//2:] ];
      beam_max = math.sqrt(max([ (abs(x)**2+abs(y)**2).max() for x,y in zip(xx,yy)]));
      _voltage_beams[self._vb_key] = vbs,beam_max;
    return vbs,beam_max;

  def get_result (self,request,*children):
    # get list of VoltageBeams
    vbs,beam_max = self.init_voltage_beams();
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
    hasfreq = False;
    for vb in vbs:
      if vb is None:
        vellsets.append(meq.vellset(meq.sca_vells(0.)));
      else:
        beam = vb.interpolate(freqaxis=self._freqaxis,**grid);
        hasfreq = hasfreq or vb.hasFrequencyAxis();
        if self.normalize and beam_max != 0:
          beam /= beam_max;
        vells = meq.complex_vells(beam.shape);
        vells[...] = beam[...];
        # make vells and return result
        vellsets.append(meq.vellset(vells));
    # create result object
    cells = request.cells if hasfreq else getattr(lm,'cells',None);
    result = meq.result(vellsets[0],cells=cells);
    result.vellsets[1:] = vellsets[1:];
    result.dims = (2,len(vellsets)//2);
    return result;



# test clause
if __name__ == "__main__":
  _verbosity.set_verbose(5);

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

