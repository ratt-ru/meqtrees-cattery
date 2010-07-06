# -*- coding: utf-8 -*-
import os
import os.path
import sys
import traceback


Locations = [ os.path.expanduser("~/Tigger"),"/usr/lib/meqtrees/Tigger" ];

def importTigger (verbose=0):
  # tries to find and import Tigger, returns SkyModel class
  try:
    import Tigger.SiameseInterface
    if verbose:
      print "Found Tigger interface in Python path";
    return Tigger.SiameseInterface.TiggerSkyModel;
  except ImportError:
    pass;
  # try to find "tigger" binary, and append its directory to search path
  for path in os.getenv('PATH').split(':') + Locations:
    path = os.path.normpath(path);  # clean off trailing "/", etc.
    if verbose>1:
      print "Looking for Tigger in",path;
    if os.access(os.path.join(path,"tigger"),os.X_OK):
      oldpath = sys.path;
      path = os.path.abspath(path);
      sys.path.insert(0,os.path.dirname(path));
      try:
        import Tigger.SiameseInterface
        if verbose:
          print "Found Tigger interface in",path;
        del sys.path[0];
        return Tigger.SiameseInterface.TiggerSkyModel;
      except ImportError:
        if verbose>1:
          traceback.print_exc();
        del sys.path[0];
        pass;
  raise ImportError,"Cannot locate Tigger package";

def TiggerSkyModel (verbose=0,**kw):
  classobj = importTigger(verbose=verbose);
  return classobj(**kw);
