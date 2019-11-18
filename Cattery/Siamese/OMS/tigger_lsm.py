# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

import os
import os.path
import sys
import traceback


Locations = [ os.path.expanduser("~/Tigger"),"/usr/lib/meqtrees/Tigger" ];

def importTigger (verbose=0):
  # tries to find and import Tigger, returns SkyModel class
  try:
    import Tigger
    import Tigger.SiameseInterface
    print("Using LSM module from Tigger (%s) %s at %s (in path)"%(Tigger.release_string,
        Tigger.svn_revision_string,os.path.dirname(Tigger.__file__)));
    return Tigger.SiameseInterface.TiggerSkyModel;
  except ImportError:
    traceback.print_exc();
  # try to find "tigger" binary, and append its directory to search path
  for path in os.getenv('PATH').split(':') + Locations:
    path = os.path.normpath(path);  # clean off trailing "/", etc.
    if verbose>1:
      print("Looking for Tigger in",path);
    if os.access(os.path.join(path,"tigger"),os.X_OK):
      oldpath = sys.path;
      path = os.path.abspath(path);
      sys.path.insert(0,os.path.dirname(path));
      try:
        import Tigger
        import Tigger.SiameseInterface
        print("Using LSM module from Tigger (%s) %s at %s"%(Tigger.release_string,
            Tigger.svn_revision_string,os.path.dirname(Tigger.__file__)));
        del sys.path[0];
        return Tigger.SiameseInterface.TiggerSkyModel;
      except ImportError:
        traceback.print_exc();
  raise ImportError("Cannot locate Tigger package");

def TiggerSkyModel (verbose=0,**kw):
  classobj = importTigger(verbose=verbose);
  return classobj(**kw);
