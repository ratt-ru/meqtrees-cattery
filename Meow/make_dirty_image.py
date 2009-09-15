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
# get image_viewer option
for i,arg in enumerate(args):
  if arg.startswith("image_viewer="):
    viewer = arg.split('=')[1];
    del args[i];
    break;
retcode = os.spawnvp(os.P_WAIT,args[0],args);

if not retcode:
  # find FITS file 
  filename = None;
  for arg in args:
    if arg.startswith("fits="):
      filename = arg.split('=')[1];
      break;
  # run visualizer
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
    print "ERROR: no image filename supplied to the make_dirty_image.py script. Please report this as a bug.";
    