# -*- coding: utf-8 -*-
import os.path
import re

from scipy import interpolate

import InterpolatedVoltageBeam
from InterpolatedVoltageBeam import *
from InterpolatedVoltageBeam import _verbosity,dprint,dprintf

def _loadPattern (filename,grid=True,rotate_xy=True,proj_theta=False):
    """Load EMSS pattern from the specified file.
    returns tuple of Ex,Ey,phi,theta,freq,gainOffset, where Ex/Ey are give the complex amplitudes in the Stokes xy frame.
    The shapes are either:
      grid=true: we have a regular grid in phi/theta space, given by the phi/theta vectors.
                 Ex/Ey have shape (nphi,ntheta).
      grid=False: irregular grid. Ex,Ey,phi,theta are all vectors of the same length
    freq is the frequency at which the beam is defined, and gainOffset is the normalized offset at center.
    """
    # Read file
    lines = file(filename).readlines()
    # Setup regexp to parse one line of pat file
    floatNum = r'([0-9.+\-E]+)'
    complexNum = r'\(\s*' + floatNum + r',\s*' + floatNum + r'\)'
    row = re.compile(r'^\s*' + r'\s+'.join([floatNum, floatNum, complexNum, complexNum]) + r'\s*$', \
                     flags=(re.IGNORECASE | re.MULTILINE))
    # Parse file data and convert to floating-point
#    for line in lines:
#      if not row.match(line):
#        print line;
    data = numpy.array(re.findall(row, ''.join(lines)), dtype='float64')
    # Find frequency
    freqExp = re.findall(r'frequency = ' + floatNum + ' MHz', ''.join(lines[:20]))
    if len(freqExp) > 0:
        freq = float(freqExp[0])*1e+6;
    else:
        freq = None;
    # Find gain expression somewhere in first few lines, if it exists
    gainExp = re.findall(r'Gain = 20 x log_10\(\|E\|\) \+ ' + floatNum, ''.join(lines[:20]))
    if len(gainExp) > 0:
        gainOffset = float(gainExp[0])
    else:
        gainOffset = 0.0
    dprint(2,"frequency",freq,"gain offset",gainOffset);

    # Extract variables
    ndata = data.shape[0];
    theta, phi = data[:, 0],data[:, 1];
    E_theta, E_phi = data[:, 2] + 1j * data[:, 3], data[:, 4] + 1j * data[:, 5]
    
    if grid:
      # Now, figure out how to put this into a cube. For now, assume that theta cycles faster,
      # and phi slower
      # check that phi is sorted
      if (phi[1:]-phi[:-1]).min() < 0:
        raise TypeError,"%s: phi axis not monotonically increasing"%filename;
      nphi = len(numpy.unique(phi));
      ntheta = ndata//nphi;
#      print filename,nphi,ntheta,nphi*ntheta,ndata;
      if nphi*ntheta != ndata:
        raise TypeError,"%s: uneven number of samples"%filename;
      # reshape arrays into nphi x ntheta array
      shape = (nphi,ntheta);
      dprint(2,"data shape is",shape);
      phi = phi.reshape(shape);
      theta = theta.reshape(shape);
      E_theta = E_theta.reshape(shape);
      E_phi   = E_phi.reshape(shape);
      # check that results are sensible: phi had better be constant along the second axis,
      # and theta along the first
      if not (phi.min(1)==phi.max(1)).all():
        raise TypeError,"%s: phi axis malformed"%filename;
      if not (theta.min(0)==theta.max(0)).all():
        raise TypeError,"%s: phi axis malformed"%filename;
      phi = phi[:,0];
      theta = theta[0,:];
      dprint(2,"phi grid is (deg)",phi);
      dprint(2,"theta grid is (deg)",theta);
      # OK, now we have a definite rectangular grid. Check that phi goes from 0 to 360 properly, and
      # if we don't have a copy of the first row (phi==0) at the end (phi==360), then add it,
      # so as to aid interpolation
      if phi[0] == 0 and phi[-1] == 360:
        # average the two rows
        for E in E_theta,E_phi:
          if not (E[0,:] == E[-1,:]).all():
            dprint(0,"warning: conflicting values at phi=0,360 -- replacing with average");
            avg = (E[0,:]+E[-1,:])/2;
            E[0,:] = E[-1,:] = avg;
            #raise ValueError,"%s: conflicting values at phi=0 and phi=360"%filename;
      else:
        phi = numpy.array(list(phi) + (phi[0]+360));
        E_theta = numpy.append(E_theta,E_theta[0,:]);
        E_phi = numpy.append(E_phi,E_phi[0,:]);
      # Convert field strengths from antenna (theta, phi) coords to Stokes (x, y) coords
      phi *= DEG;
      theta *= DEG;
      cos_theta = numpy.cos(theta)[numpy.newaxis,:] if proj_theta else 1;
      cos_phi   = numpy.cos(phi)[:,numpy.newaxis] if rotate_xy else 1;
      sin_phi   = numpy.sin(phi)[:,numpy.newaxis] if rotate_xy else 0;
      Ex =  E_theta * cos_theta * cos_phi - E_phi * sin_phi;
      Ey = -E_theta * cos_theta * sin_phi - E_phi * cos_phi;

      # replace the theta=0 point with the average across the entire array (otherwise we have a discontinuity)
      if theta[0] == 0:
        Ex[:,0] = Ex[:,0].mean();
        Ey[:,0] = Ey[:,0].mean();
      
    # else generate ungridded data
    else:
      # Convert field strengths from antenna (theta, phi) coords to Stokes (x, y) coords
      # (NB: must do this before averaging at theta=0!!)
      cos_theta = numpy.cos(theta*DEG);
      cos_phi   = numpy.cos(phi*DEG);
      sin_phi   = numpy.sin(phi*DEG);
      Ex =  E_theta * cos_theta * cos_phi - E_phi * sin_phi;
      Ey = -E_theta * cos_theta * sin_phi - E_phi * cos_phi;
      # The EMSS pattern contains multiple versions of the beam center, where theta=0
      # but phi differs. Replace them with a single averaged version.
      beamCenter = (theta == 0.0)
      Ex[beamCenter] = Ex[beamCenter].mean()
      Ey[beamCenter] = Ey[beamCenter].mean()
      # select theta>0, phi<360
      select = (theta > 0.0) & (theta <= 90.0) & (phi < 360.0);
      # add a single theta=0 point to the selection
      select[numpy.where(beamCenter)[0][0]] = True
      # apply selection
      phi   = phi[select];
      theta = theta[select];
      Ex    = Ex[select];
      Ey    = Ey[select];
      # Convert field strengths from antenna (theta, phi) coords to Stokes (x, y) coords
      phi *= DEG;
      theta *= DEG;
      
    return Ex,Ey,phi,theta,freq,gainOffset;

def loadPatternSet (filenames,y=False,rotate_xy=True,proj_theta=False):
  """Loads a set of per-frequency patterns. Returns a tuple of E,phi,theta,freq.
  E is an (nphi,ntheta,nfreq) cube of complex gains (Ex if y=False, else Ey), while
  phi, theta and freq are vectors of coordinates.""";
  freqs = [];
  for ifreq,filename in enumerate(filenames):
    Ex,Ey,phi,theta,freq,gainOffset = _loadPattern(filename,grid=True,rotate_xy=rotate_xy,proj_theta=proj_theta);
    # change order of axis, since FITS has first axis last
    beam = Ey if y else Ex;
    freqs.append(freq);
    # check if it matches previous image
    if not ifreq:
      baseshape = beam.shape;
      phi0,theta0 = phi,theta;
      # setup 3D beam cube
      beamcube = numpy.zeros((len(phi),len(theta),len(filenames)),complex);
    else:
      if baseshape != beam.shape:
        raise TypeError,"file %s has different dimensions"%filename;
      if not((phi==phi0).all() and (theta==theta0).all()):
        raise TypeError,"file %s has different grid"%filename;
    beamcube[:,:,ifreq] = beam[:,:];
  dprint(2,"beam array has shape",beamcube.shape);
  dprint(2,"frequencies are",freqs);
  return beamcube,phi0,theta0,numpy.array(freqs);

class EMSSVoltageBeamPS (InterpolatedVoltageBeam):
  """This class implements a complex voltage beam that is read from an EMSS pattern file.
  This uses map_coordinates to interpolate values in polar coordinates (i.e. l/m inputs
  are converted to phi/theta, and interpolated in that grid)."""
  def __init__ (self,filenames,y=True,hier_interpol=True,spline_order=3,theta_step=1,phi_step=1,rotate=0,
      rotate_xy=True,proj_theta=False,normalization_factor=1,verbose=0):
    InterpolatedVoltageBeam.__init__(self,spline_order=spline_order,hier_interpol=hier_interpol);
    self._theta_step = theta_step;
    self._phi_step = phi_step;
    self._rotate = rotate*DEG;
    self._freq_warning = False;
    _verbosity.set_verbose(verbose);
    if verbose:
      _verbosity.enable_timestamps(True);
    self.read(filenames,y=y,rotate_xy=rotate_xy,proj_theta=proj_theta,normalization_factor=normalization_factor);
    
  def read (self,filenames,y=True,rotate_xy=True,proj_theta=False,normalization_factor=1):
    """Reads beam patterns from EMSS files. If y=True, uses the second (Ey) column""";
    beamcube,phi0,theta0,freqs = loadPatternSet(filenames,y=y,rotate_xy=rotate_xy,proj_theta=proj_theta);
    if self._theta_step != 1 or self._phi_step != 1:
      beamcube = beamcube[::self._phi_step,::self._theta_step,:];
      phi0 = phi0[::self._phi_step];
      theta0 = theta0[::self._theta_step];
    if normalization_factor != 1:
      beamcube /= normalization_factor;
    self.setFreqGrid(freqs);
    if len(freqs) > 1:
      self._freqmap = interpolate.interp1d(freqs,range(len(freqs)),'linear',bounds_error=False,fill_value=numpy.nan);
    else:
      self._freqmap = None;
    self._phi0 = phi0;
    self._theta0 = theta0;
    # note that coordinates out of the given phi/theta range will be converted to NANs
    self._phimap = interpolate.interp1d(phi0,range(len(phi0)),'linear',bounds_error=False,fill_value=numpy.nan);
    self._thetamap = interpolate.interp1d(theta0,range(len(theta0)),'linear',bounds_error=False,fill_value=numpy.nan);
    self.setBeamCube(beamcube);
  
  def hasFrequencyAxis (self):
    return len(self.freqGrid()) > 1;
  
  def freqToBeam (self,freq):
    if self._freqmap is None:
      raise RuntimeError,"attempting to interpolate in frequency, but frequency axis is not set. This is a bug!";
    # interpolate from grid back to linear 0...n-1 range
    freqcoord = self._freqmap(freq);
    # out-of-bounds values filled by NANs -- get a mask of these, issue warning, and reset to 0 and N-1
    oob = numpy.isnan(freqcoord);
    if oob.any():
      freqgrid = self.freqGrid();
      if not self._freq_warning:
        dprint(0,"WARNING: requested frequencies (%s) are outside the supplied beam frequency range (%f to %f MHz)"%
          (",".join(["%f"%(x*1e-6) for x in freq[oob]]),freqgrid[0]*1e-6,freqgrid[-1]*1e-6));
        dprint(0,"Doing frequency extrapolation instead");
        self._freq_warning = True;
      freqcoord[freq<=freqgrid[0]]  = - freq/freqgrid[0]     # lm multiplied by negative of this 
      freqcoord[freq>=freqgrid[-1]] = - freq/freqgrid[-1]   
    return freqcoord;

  @staticmethod
  def _normalizePhi (phi):
    # make sure phi is in 0-360 range
    phi = numpy.fmod(phi,2*math.pi);
    phi[phi<0] += 2*math.pi;
    return phi;

  def lmToBeam (self,l,m,rotate=None):  
    # m is direction sine, m1 is cosine
#    dprint(4,"lmToPolar",l/DEG,m/DEG);
    dprint(3,"lmToPolar [0]",l.ravel()[0]/DEG,m.ravel()[0]/DEG);
    l1 = numpy.sqrt(1-l**2);
    m1 = numpy.sqrt(1-m**2);
    phi = -numpy.arctan2(l,m) + self._rotate;  # flip phi, so it goes clockwise from N=0 to W=90
    if rotate is not None:
      phi += rotate;
    phi = self._normalizePhi(phi);
    theta = numpy.arccos(l1*m1);
    dprint(3,"phi,theta [0]",phi.ravel()[0]/DEG,theta.ravel()[0]/DEG);
#    dprint(4,"phi,theta are",phi/DEG,theta/DEG);
    return self._phimap(phi),self._thetamap(theta);
    
  def thetaPhiToBeam (self,theta,phi,rotate=None):
    # apply rotations to phi
    if self._rotate:
      phi = phi + self._rotate;
    if rotate is not None:
      phi = phi + rotate;
    # make sure phi is in 0-360 range
    phi = self._normalizePhi(phi);
    dprint(3,"phi,theta [0]",phi.ravel()[0]/DEG,theta.ravel()[0]/DEG);
    # interpolate
    return self._phimap(phi),self._thetamap(theta);
    
class EMSSVoltageBeamGridder (object):
  """This class implements a complex voltage beam that is read from an EMSS pattern file and gridded
  onto a regular l/m grid using matplotlib.mlab.griddata.
  """ 
  def __init__ (self,filenames,y=True,hier_interpol=True,spline_order=3,theta_step=1,phi_step=1,rotate=0,verbose=0):
    self._spline_order = spline_order;
    self._ampl_interpol = True;
    self._theta_step = theta_step;
    self._rotate = rotate*DEG;
    _verbosity.set_verbose(verbose);
    self.read(filenames,y=y);
    
  def read (self,filenames,y=True):
    """Reads beam patterns from EMSS files. If y=True, uses the second (Ey) column""";
    self._beams = [];
    for ifreq,filename in enumerate(filenames):
      Ex,Ey,phi,theta,freq,gainOffset = _loadPattern(filename,grid=False);
      beamcube = Ey if y else Ex;
      # convert coordinates and make arrays
      sin_theta = numpy.sin(theta);
      lcoord = sin_theta*numpy.sin(phi);
      mcoord = sin_theta*numpy.cos(phi);
      beam = Ey if y else Ex;
      beam_ampl = abs(beam) if self._ampl_interpol else None;
      self._beams.append((beam,beam_ampl,lcoord,mcoord,freq));
    self.freq = [ x[4] for x in self._beams ];
    
  def freqGrid (self):
    return self.freq;

  def grid (self,lgrid,mgrid,freqgrid):
    from matplotlib.mlab import griddata
    l = numpy.sin(lgrid);
    m = numpy.sin(mgrid);
    gridshape = (len(l),len(m),len(self.freq));
    beamgrid = numpy.zeros(gridshape,complex);
    # figure out useful range of l/m coordinates
    lmin,lmax = l.min(),l.max();
    mmin,mmax = m.min(),m.max();
    dl = (lmax - lmin)/10;
    dm = (mmax - mmin)/10;
    lmin,lmax = lmin-dl,lmax+dl;
    mmin,mmax = mmin-dm,mmax+dm;
    dprint(3,"regridding in area from %f,%f to %f,%f deg"%(lmin/DEG,mmin/DEG,lmax/DEG,mmax/DEG));
    # loop over all per-frequency patterns
    for ifreq,(beam,beam_ampl,lcoord,mcoord,freq) in enumerate(self._beams):
      # filter out values outside the range
      select = (lcoord < lmax)&(lcoord > lmin)&(mcoord < mmax)&(mcoord > mmin); 
      dprint(3,"selection reduces beam from %d to %d points"%(len(lcoord),select.sum()));
      # regrid onto each frequency plane
      for ifreq in range(len(self.freq)):
        gbeam = beamgrid[:,:,ifreq];
        gbeam.real = griddata(lcoord[select],mcoord[select],beam.real[select],l,m);
        gbeam.imag = griddata(lcoord[select],mcoord[select],beam.imag[select],l,m);
        if beam_ampl is not None:
          ampl = griddata(lcoord[select],mcoord[select],beam_ampl[select],l,m);
          ampl1 = abs(gbeam);
          wh = ampl1 != 0;
          gbeam[wh] = (gbeam[wh]/ampl1[wh])*ampl[wh];
      # replace nans with zeroes
      beamgrid[numpy.isnan(beamgrid)] = 0;
    # interpolate onto frequency grid
    dprint(4,"interpolated data range",beamgrid.min(),beamgrid.max());
    dprint(3,"interpolating frequencies");
    l1,m1 = numpy.meshgrid(numpy.arange(len(l)),numpy.arange(len(m)));
    if len(self.freq) > 1:
      f1 = interpolate.interp1d(self.freq,range(len(self.freq)),'linear')(freqgrid);
    else:
      f1 = numpy.array(0.);
    l1,f1 = unite_shapes(l1,f1);
    m1,f1 = unite_shapes(m1,f1);
    coords = numpy.vstack((l1.ravel(),m1.ravel(),f1.ravel()));
    output_shape = (len(l),len(m),len(freqgrid));
    output = numpy.zeros(output_shape,complex);
    output.real = interpolation.map_coordinates(beamgrid.real,coords,order=self._spline_order).reshape(output_shape);
    output.imag = interpolation.map_coordinates(beamgrid.imag,coords,order=self._spline_order).reshape(output_shape);
    dprint(3,"done");
    return output;

if __name__ == "__main__":
  import sys
  _verbosity.set_verbose(5);

  vb = EMSSVoltageBeam(spline_order=3);
  vb.read(('Gregorian_h_1000_full.pat','Gregorian_h_1900_full.pat'));
  
  l0 = numpy.array([-2,0,2])*DEG;
  l = numpy.vstack([l0]*len(l0));

  a = vb.interpolate(l,l.T,freq=[1.456e+9],freqaxis=2);
  b = vb.interpolate(l,l.T,freq=[1.456e+9,1.457e+9,1.458e+9],freqaxis=2);
  c = vb.interpolate(l,l.T,freq=[1.455e+9,1.457e+9,1.458e+9,1.46e+9],freqaxis=2);

  print "C",c.shape,c;
  print "B",b.shape,b;
  print "A",a.shape,a;
  
  freq0 = 1e+9;
  dfreq = 0.1e+9;
  nfreq = 8;
  dp = 0.01;
  radius = 500;
  
  l0 = numpy.arange(-radius,radius+1)*dp*DEG;
  l = numpy.vstack([l0]*len(l0));
  freq = numpy.arange(nfreq)*dfreq + freq0;
  c = vb.interpolate(l,l.T,freq=freq,freqaxis=2);
  
  # create FITS file
  import pyfits;
  x = abs(c.transpose());
  hdu = pyfits.PrimaryHDU(x);
  hdr = hdu.header;
  hdr.update('CTYPE1','L',"");
  hdr.update('CRPIX1',radius+1,"");
  hdr.update('CRVAL1',0,"");
  hdr.update('CDELT1',dp,"");
  hdr.update('CUNIT1','DEG',"");
  hdr.update('CTYPE2','M',"");
  hdr.update('CRPIX2',radius+1,"");
  hdr.update('CRVAL2',0,"");
  hdr.update('CDELT2',dp,"");
  hdr.update('CUNIT2','DEG',"");
  hdr.update('CTYPE3','FREQ',"");
  hdr.update('CRPIX3',1,"");
  hdr.update('CRVAL3',freq0,"");
  hdr.update('CDELT3',dfreq,"");
  hdr.update('CUNIT3','HZ',"");
  
  hdu.writeto('tmp.fits',clobber=True);
  
  
