#!/usr/bin/python
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
  args.append("cachesize=512");

retcode = os.spawnvp(os.P_WAIT,args[0],args);

if not retcode:
  # find FITS file and run visualizer
  filename = None;
  for arg in args:
    if arg.startswith("fits="):
      filename = arg.split('=')[1];
      break;
  if filename and os.path.exists(filename):
    os.execvp("kvis",["kvis",filename]);