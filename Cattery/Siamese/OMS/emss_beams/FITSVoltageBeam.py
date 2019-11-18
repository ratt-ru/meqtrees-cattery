# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division
import os.path

## ugly hack to get around UGLY FSCKING ARROGNAT (misspelling fully intentional) pyfits-2.3 bug
import Kittens.utils
pyfits = Kittens.utils.import_pyfits

import warnings
print(pyfits.formatwarning,warnings.formatwarning)

from scipy import interpolate

from .InterpolatedVoltageBeam import *
from .InterpolatedVoltageBeam import _verbosity

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
    if iaxis < 0:
      raise TypeError("missing axis '%s'"%axis);
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
    pix = self._rpix[iaxis] + (world - self._rval[iaxis])/self._delta[iaxis]; #float px as per FITS def?
#    print "toPixel",axis,world,pix;
    return pix

  def toWorld (self,axis,pixel):
    """Converts array of pixel coordinates to world coordinates""";
    iaxis = self.iaxis(axis);
    return (pixel - self._rpix[iaxis])*self._delta[iaxis] + self._rval[iaxis];

  def grid (self,axis):
    iaxis = self.iaxis(axis);
    return self.toWorld(iaxis,numpy.arange(0.,float(self._naxis[iaxis])));

class FITSVoltageBeam (InterpolatedVoltageBeam):
  """This class implements a complex voltage beam that is read from a FITS file."""

  def hasFrequencyAxis (self):
    return self._freqToPixel is not None;

  def lmToBeam (self,l,m,rotate=None):
    if rotate is not None:
      c,s = numpy.cos(rotate),numpy.sin(rotate);
      return self._lToPixel(l*c + m*s),self._mToPixel(-l*s + m*c);
    else:
      return self._lToPixel(l),self._mToPixel(m);

  def freqToBeam (self,freq):
    if self._freqToPixel is None:
      raise RuntimeError("attempting to interpolated in frequency, but frequency map is not set. This is a bug!");
    return self._freqToPixel(freq);

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
    laxis = axes.iaxis('L');
    maxis = axes.iaxis('M');
    if laxis<0 or maxis<0:
      raise TypeError("FITS file %s missing L or M axis"%filename_real);
    # setup conversion functions
    self._lToPixel = Kittens.utils.curry(axes.toPixel,laxis);
    self._mToPixel = Kittens.utils.curry(axes.toPixel,maxis);
    # find frequency grid. self._freqToPixel will be None if no frequency axis
    freqaxis = axes.iaxis('FREQ');
    if freqaxis >= 0 and axes.naxis(freqaxis) > 1:
      dprint(1,"FREQ axis has %d points"%axes.naxis(freqaxis));
      self.setFreqGrid(axes.grid(freqaxis));
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
    dprint(1,"beam array has shape",beam.shape);
    dprint(2,"l grid is",axes.grid(laxis));
    dprint(2,"m grid is",axes.grid(maxis));
    if freqToPixel:
      dprint(2,"freq grid is",axes.grid(freqaxis));
    # set the beam cube
    self.setBeamCube(beam);
    
    
class MultifreqFITSVoltageBeam (FITSVoltageBeam):
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
        lmap = Kittens.utils.curry(axes.toPixel,laxis);
        mmap = Kittens.utils.curry(axes.toPixel,maxis);
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
    self.setFreqGrid(freqs);
    freqmap = interpolate.interp1d(freqs,list(range(len(freqs))),'linear');
    # prefilter beam for interpolator
    self.setBeamCube(beamcube,lmap=lmap,mmap=mmap,freqmap=freqmap);
    

if __name__ == "__main__":
  import sys
  _verbosity.set_verbose(5);

  vb = MultifreqFITSVoltageBeam(spline_order=3);
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

  vb = FITSVoltageBeam(spline_order=3);
  vb.read("mk_1455MHz_feed_5deg_xxreal.fits","mk_1455MHz_feed_5deg_xxreal.fits");

  l0 = numpy.array([-2,-1,0,1,2])*DEG;
  l = numpy.vstack([l0]*len(l0));

  print(vb.interpolate(l,l.T));

