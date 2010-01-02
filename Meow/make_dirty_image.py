#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import os
import os.path

args = list(sys.argv);
args[0] = 'lwimager';
# insert cachesize option
for arg in args:
  if arg.startswith("cachesize="):
    break;
else:
  args.append("cachesize=4096");

viewer = 'kvis';
remove_img = False;
# get image_viewer and remove_image options; remove them from arguments as they don't get
# passed on to lwimager
for i,arg in enumerate(args):
  if arg.startswith("image_viewer="):
    viewer = arg.split('=')[1];
    del args[i];
  elif arg=="remove_image":
    remove_img = True;
    del args[i];


# run the imager
retcode = os.spawnvp(os.P_WAIT,args[0],args);

if not retcode:
  # get names of files
  files = {};
  for arg in args:
    for argname in ('image','fits','model','residual','restored'):
      if arg.startswith(argname+"="):
        files[argname] = arg.split('=')[1];
  
  # convert img files produced by CLEANing to FITS
  for filetype in ('model','residual'):  # restored image is already output via the 'fits' argument
    filename = files.get(filetype);
    if filename:
      fitsname = os.path.splitext(filename)[0]+".fits";
      if os.path.exists(fitsname):
        try:
          os.unlink(fitsname);
        except:
          pass;
      os.system("image2fits in=\\\\%s out=%s"%(filename,fitsname));

  # remove img files, with the exception of the CLEAN model
  if remove_img:
    for filetype in ('image','residual','restored'):
      filename = files.get(filetype);
      if filename:
        os.system("rm -fr %s"%filename);

  # run visualizer if a fits file is available
  filename = files.get('fits');
  if filename:
    if viewer == "none":
      print "Image %s ready. No image viewer selected (see TDL Exec menu)."%filename;
    else:
      if filename and os.path.exists(filename):
        try:
          os.execvp(viewer,[viewer,filename]);
        except OSError:
          print "WARNING: failed to start image viewer '%s'. Perhaps it is not installed?"%viewer;
  else:
    print "ERROR: no FITS image filename supplied to the make_dirty_image.py script. Please report this as a bug.";
