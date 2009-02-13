Version = "0.1";

import os.path

# init debug printing
import Timba.utils
_verbosity = Timba.utils.verbosity(name="skypuss");
dprint = _verbosity.dprint;
dprintf = _verbosity.dprintf;

import Timba.GUI.pixmaps 
pixmaps = Timba.GUI.pixmaps.PixmapCache("purr");


import Timba.Apps.config
Config = Timba.Apps.config.section("SkyPuss");


# if GUI is enabled, this will be overwritten by a class that sets a busy cursor 
# on instantiation, and restores the cursor when deleted
class BusyIndicator (object):
  pass;

def progressMessage (msg,sub=False):
  """shows a progress message. If a GUI is available, this will be redefined by the GUI.
  If sub is true, message is a sub-message of previous one.""";
  pass;
