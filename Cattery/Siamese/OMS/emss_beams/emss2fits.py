#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division
if __name__ == "__main__":

  # setup some standard command-line option parsing
  #
  import sys
  import numpy
  import math
  import os.path
  DEG = math.pi/180;
  from optparse import OptionParser,OptionGroup

  
  parser = OptionParser(usage="""%prog: [options] <*.pat file> <FITS file>""",
      description=
"""Converts EMSS ASCII beam patterns to FITS files. Multiple frequency patterns and griddings may be given by embedding
the strings @freq and @label in the pattern filename, and specifying the --freqs and --labels options.
For a compound pattern built from a  sum of components (i.e. feed, main reflector, subreflector), 
separate the per-component groups with '+'.""");

  parser.add_option("-x",action="store_true",
                    help="generate the X component (default)");
  parser.add_option("-y",action="store_true",
                    help="generate the Y beam component");
  parser.add_option("--freqs",metavar="FREQ1,FREQ2,...",type="string",default="",
                    help="frequency labels, to be substituted for @freq in the filename");
  parser.add_option("--labels",metavar="LABEL1,LABEL2,...",type="string",default="",
                    help="component labels, to be substituted for @label in the filename");
  parser.add_option("--rotate",metavar="DEG",type="float",default=0,
                    help="rotate the pattern by the specified number of degrees.");
  parser.add_option("--no-rotate-xy",action="store_false", dest="rotate_xy", default=True,
                    help="do not rotate data from theta/phi to Stokes coordinate system; default is to do the rotation");
  parser.add_option("--ampl",action="store_true",
                    help="generate the amplitude pattern.");
  parser.add_option("--real",action="store_true",
                    help="generate the real pattern.");
  parser.add_option("--imag",action="store_true",
                    help="generate the imaginary pattern. Default generates all three.");
#  parser.add_option("-i","--interpolator",choices=("ps","grid"),default="ps",
#                    help="select interpolation algorithm: 'ps' (polar spline), 'grid' (gridder). Default is %default.");
#  parser.add_option("-A","--no-ampl-interpol",action="store_true",
#                    help="disable amplitude interpolation");
  parser.add_option("-s","--spline-order",metavar="ORDER",type='int',default=3,
                    help="spline order, 1 to 3. Default is %default.");
  parser.add_option("--scale",type="float",default=1.,
                    help="rescale the beam (in size) by the given factor, >1 makes the beam broader.");
  parser.add_option("--theta-stepping",metavar="N",type="int",default=1,
                    help="sparsify the theta grid by a factor of N (i.e. only use every Nth sample for interpolation).");
  parser.add_option("--phi-stepping",metavar="N",type="int",default=1,
                    help="sparsify the phi grid by a factor of N.");
  parser.add_option("-t","--template",metavar="FILENAME",type="string",
                    help="generate the beam using the coordinate system (in l,m and frequency) of "+
                    "the given FITS file. The FITS file must have axes of CTYPE RA-*, DEC-* and FREQ. Note that "+
                    "this option overrides --radius, --resolution, --frequency, --num-freq and --delta-freq.");
  parser.add_option("-r","--radius",metavar="PIXELS",type="int",default=512,
                    help="radius of FITS file, in pixels. Resulting NAXISn will be 2*radius+1.");
  parser.add_option("-R","--resolution",metavar="DEG",type="float",default=0.01,
                    help="resolution of FITS file, in degrees. Default is %default.");
  parser.add_option("-f","--frequency",metavar="MHz",type="float",default=0,
                    help="first frequency plane, MHz.");
  parser.add_option("-n","--num-freq",type="int",default=1,
                    help="number of frequency planes. Default is %default.");
  parser.add_option("-d","--delta-freq",metavar="MHz",type="float",default=1,
                    help="frequency stepping, in MHz. Default is %default.");
  parser.add_option("--taper",metavar="FROM_DEG,TO_DEG",type="string",
                    help="apply a smooth taper to the beam. The taper function goes from 1 at radius FROM_DEG to 0 at radius TO_DEG");
  parser.add_option("--incremental",metavar="N",type="int",default=0,
                    help="compute the patterns in chunks of N grid rows at a time. Use this with very large "+
                    "patterns, if the full computation makes you run out of memory.");
  parser.add_option("-v","--verbose",metavar="LEVEL",type="int",default=0,
                    help="verbosity level, higher for more messages.");
  parser.add_option("-T","--timestamps",action="store_true",
                    help="include timestamps with messages.");

  (options,args) = parser.parse_args();
  
  if len(args) != 2:
    parser.error("One input *.pat file and one output .fits file must be given");
  
  patfile,fitsfile = args;
  
  if options.taper:
    try:
      taper_r0,taper_r1 = numpy.array(list(map(float,options.taper.split(","))))*DEG;
    except:
      parser.error("Invalid --taper %s setting"%options.taper);
      
  
  # get basename of FITS file, if all three need to be generated
  nout = sum(map(bool,(options.ampl,options.real,options.imag)));
  if nout > 1:
    parser.error("Only one of --ampl/--real/--imag may be specified");
  # if nout=0, generate three FITS files, so take the basename  
  elif nout == 0:
    if fitsfile.upper().endswith('.FITS'):
      fitsfile = os.path.splitext(fitsfile)[0];
  # else just 1 file, make sure it has fits extension
  else:
    if not fitsfile.upper().endswith('.FITS'):
      fitsfile += ".fits";
    
  from .EMSSVoltageBeam import _verbosity,EMSSVoltageBeamPS,EMSSVoltageBeamGridder,dprint,dprintf
  _verbosity.set_verbose(options.verbose);
  _verbosity.enable_timestamps(options.timestamps);
  interpolator = EMSSVoltageBeamPS;
  
#  if options.interpolator == "ps":
#    interpolator = EMSSVoltageBeamPS;
#  else:
#    interpolator = EMSSVoltageBeamGridder;
#  print "Using interpolator %s"%(interpolator.__name__);
#    "w/o amplitude interpolation" if options.no_ampl_interpol else "with amplitude interpolation");
  
  # how many compound beams do we have
  dprint(1,"reading pattern files");
  beam_components = [];
  labels = options.labels.split(",") or [""];
  freqs = options.freqs.split(",") or [""];
  for label in labels:
    patfile1 = patfile.replace("@label",label);
    patfile2 = [ patfile1.replace("@freq",freq) for freq in freqs ]; 
    vb = interpolator(patfile2,y=options.y,
                      spline_order=options.spline_order,
                      rotate_xy=options.rotate_xy,
                      theta_step=options.theta_stepping,phi_step=options.phi_stepping,rotate=options.rotate,
                      verbose=options.verbose);
    print("Creating %s voltage beam from %s"%("Y" if options.y else "X"," ".join(patfile2)));
    beam_components.append(vb);
  
  if not beam_components:
    parser.error("At least one input *.pat file must be given");
  
  dprint(1,"interpolating beam patterns");

  # if template FITS file is provided, get l/m/frequency grid from there
  l0 = m0 = freq = None;
  
  import pyfits
  if options.template:
    ff = pyfits.open(options.template);
    hdr = ff[0].header;
    # get axis info
    ndim = hdr['NAXIS'];
    for iaxis in range(ndim):
      axs = str(iaxis+1);
      npix = hdr['NAXIS'+axs];
      ctype,crval,cdelt,crpix = hdr.get('CTYPE'+axs,""),hdr.get('CRVAL'+axs,0),hdr.get('CDELT'+axs,0),hdr.get('CRPIX'+axs,0);
      name = ctype.strip().upper();
      # get the coordinate (relative to CRVALx) of each image pixel along this axis
      pixcoord = (numpy.arange(1,npix+1) - crpix)*cdelt;
      # have we found the coordinate axes?
      if name.startswith("RA"):
        l0 = numpy.sin(pixcoord*DEG);
        axis_ra = ctype,npix,crval,cdelt,crpix;
      elif name.startswith("DEC"):
        m0 = numpy.sin(pixcoord*DEG);
        axis_dec = ctype,npix,crval,cdelt,crpix;
      elif name == "FREQ":
        freq = crval + pixcoord;
        axis_freq = ctype,npix,crval,cdelt,crpix;
    # have we found l/m axes?
    if l0 is not None and m0 is not None and freq is not None:
      print("Using %dx%d grid from template file %s"%(len(l0),len(m0),options.template));
      print("Using frequencies of %sMHz"%(",".join(map(str,freq*1e-6))));
      l0,m0 = numpy.meshgrid(l0,m0);
    else:
      print("Template FITS file %s does not contain RA/DEC/FREQ axes"%options.template);
      sys.exit(1);
      
  # else generate an lm grid based on the command line
  else:
    # create coordinates for interpolation
    # here we make a regularly-gridded image in degrees
    l0 = numpy.arange(-options.radius,options.radius+1)*options.resolution*DEG/options.scale;
    l0,m0 = numpy.meshgrid(l0[::-1],l0);
## old code for when we used to use the grid interpolator too
#  elif options.interpolator == "grid":
#    l0 = numpy.arange(-options.radius,options.radius+1)*options.resolution*DEG/options.scale;
#    images = [ vb.grid(l0,l0,freq) for vb in beam_components ];

  # note that the interpolator classes expect direction cosines (i.e. sines w.r.t. centre), so take the sines now
  l0 = numpy.sin(l0);
  m0 = numpy.sin(m0);
  
  if freq is None:
    # if frequency was not loaded from template, pick the first one 
    freq0 = options.frequency;
    if options.frequency:
      freq0 = options.frequency*1e+6;
    else:
      freq0 = beam_components[0].freqGrid()[0];
    freq = numpy.arange(options.num_freq)*options.delta_freq*1e+6 + freq0;
    print("Using frequency grid of %sMHz"%(",".join(map(str,freq))));
  
    
  # INTERPOLATE!!!
  if options.incremental:
    parts = [];
    for row0 in range(0,len(l0),options.incremental):
      row1 = min(row0+options.incremental,len(l0));
      print("Generating pattern for rows %d:%d"%(row0,row1));
      images = [ vb.interpolate(l0[row0:row1,...],m0[row0:row1,...],freq=freq,freqaxis=2) 
                 for vb in beam_components ];
      parts.append(sum(images));
    print("Concatenating patterns");
    image = numpy.concatenate(parts);
    parts = None;
  else:
    dprint(1,"interpolating components");
    images = [ vb.interpolate(l0,m0,freq=freq,freqaxis=2) for vb in beam_components ];
    dprint(1,"summing components");
    image = sum(images);
  dprint(1,"generating FITS images");
  
  # apply taper
  if options.taper:
    r0 = numpy.sqrt(l0*l0+m0*m0).reshape(image.shape);
    taper = numpy.cos((r0-taper_r0)/(taper_r1-taper_r0)*math.pi/2);
    taper[r0<=taper_r0] = 1;
    taper[r0>=taper_r1] = 0;
    image *= taper;
  
  # create FITS file
#  image = image.transpose();

# reshape image into FORTRAN order (for FITS)
  # if freq axis missing, add a dummy one
  if image.ndim < 3:
    shape = [1] + list(image.shape);
  # elsem transpose so freq axis comes first (i.e. last in header)
  else:
    image = image.transpose((2,0,1))
    shape = image.shape;
  image = image.reshape(shape,order='F');
  
  
  import os
  import os.path
  
  if options.template:
    print("Generating FITS images using template %s"%options.template)
    
    def make_image (filename,data,what):
      hdu = pyfits.PrimaryHDU(data);
      hdr = hdu.header;
      hdr.update('CTYPE1',axis_ra[0],"");
      hdr.update('CRVAL1',axis_ra[2],"");
      hdr.update('CDELT1',axis_ra[3],"");
      hdr.update('CRPIX1',axis_ra[4],"");
      hdr.update('CUNIT1','DEG',"");
      
      hdr.update('CTYPE2',axis_dec[0],"");
      hdr.update('CRVAL2',axis_dec[2],"");
      hdr.update('CDELT2',axis_dec[3],"");
      hdr.update('CRPIX2',axis_dec[4],"");
      hdr.update('CUNIT2','DEG',"");
      
      hdr.update('CTYPE3',axis_freq[0],"");
      hdr.update('CRVAL3',axis_freq[2],"");
      hdr.update('CDELT3',axis_freq[3],"");
      hdr.update('CRPIX3',axis_freq[4],"");
      hdr.update('CUNIT3','HZ',"");
      
      if os.path.exists(filename):
        os.unlink(filename);
      hdu.writeto(filename,clobber=True);
      print("Wrote %s to %s"%(what,filename));

  else:
    npix = options.radius*2 + 1;
    print("Generating FITS images of %dx%d pixels, at a resolution of %g deg/pix"%(npix,npix,options.resolution));
    print("Full image size will be %fx%f deg"%(options.resolution*npix,options.resolution*npix));
    print("Frequency axis: %s MHz"%",".join(map(str,freq*1e-6)));
    
    def make_image (filename,data,what):
      hdu = pyfits.PrimaryHDU(data);
      hdr = hdu.header;
      hdr.update('CTYPE1','L',"");
      hdr.update('CRPIX1',options.radius+1,"");
      hdr.update('CRVAL1',0,"");
      hdr.update('CDELT1',-options.resolution,"");
      hdr.update('CUNIT1','DEG',"");
      hdr.update('CTYPE2','M',"");
      hdr.update('CRPIX2',options.radius+1,"");
      hdr.update('CRVAL2',0,"");
      hdr.update('CDELT2',options.resolution,"");
      hdr.update('CUNIT2','DEG',"");
      hdr.update('CTYPE3','FREQ',"");
      hdr.update('CRPIX3',1,"");
      hdr.update('CRVAL3',freq0,"");
      hdr.update('CDELT3',options.delta_freq*1e+6,"");
      hdr.update('CUNIT3','HZ',"");
      hdr.update('EMSSROT',options.rotate);
      if os.path.exists(filename):
        os.unlink(filename);
      hdu.writeto(filename,clobber=True);
      print("Wrote %s to %s"%(what,filename));
  
  if not nout or options.ampl:
    make_image(fitsfile if options.ampl else fitsfile+"_ampl.fits",abs(image),"amplitude pattern");
  if not nout or options.real:
    make_image(fitsfile if options.real else fitsfile+"_real.fits",image.real,"real pattern");
  if not nout or options.imag:
    make_image(fitsfile if options.imag else fitsfile+"_imag.fits",image.imag,"imaginary pattern");
