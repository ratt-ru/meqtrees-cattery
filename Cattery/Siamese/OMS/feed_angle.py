# -*- coding: utf-8 -*-
"""<P>This module implements a feed rotation Jones term. The term is a simple rotation matrix:</P>

<P align="center"><I>P = Rot(&rho;)</I></P>

<P>...where the rotation angle &rho; = &eta; + &rho;<sub>q</sub> is a sum of the (optional) parallactic angle
&eta;, and a per-station feed orientation &rho;<sub>q</sub>.</P>

<P>A feed orientation of 0 corresponds to the conventional Measurement Equation definition of X pointing North on the sky, and Y pointing West.</P>

<P align="right">Author: O. Smirnov &lt;<tt>smirnov@astron.nl</tt>&gt;</P>""";

__default_label__ = "P";
__default_name__  = "parallactic angle and feed orientation";

from Timba.TDL import *
from Meow import Context,MSUtils

import math
import cmath
import traceback
import os.path

DEG = math.pi/180;

def _read_ms_feed (enabled):
  """Called when the read_ms option is changed. Reads the MS/FEED subtable if enabled=True (and an MS is already set up.)"""
  if enabled and Context.mssel and Context.mssel.msname and MSUtils.TABLE:
    try:
      # read angles and antenna IDs from table
      feed = MSUtils.TABLE(os.path.join(Context.mssel.msname,"FEED"));
      angle_col = feed.getcol('RECEPTOR_ANGLE');
      ant_col = list(feed.getcol('ANTENNA_ID'));
      feed = None;
      # form up list of angles
      angles = [0]*len(Context.mssel.ms_antenna_names);
      for iant in range(len(angles)):
        # find first row with this antenna ID. We only look at the first entry per antenna. In princple the FEED table can contain a number of entires per antenna
        # (i.e. for different times, spectral windows, etc.), but I can't be bother to worry about such exotic configurations for now. Instead, we just get the first feed
        # angle for every antenna, or use 0 if none is found
        try:
          irow = ant_col.index(iant);
        except ValueError:
          continue;
        # get the X angle, as the first number from the column. We ignore the second number for now (tacitly assuming Y is orthogonal)
        angles[iant] = list(angle_col[irow])[0];
      # now make an angle string. If all angles are the same, use a single value
      if len(set(angles)) == 1:
        angle_opt.set_value("%g"%(angles[0]/DEG));
      # else join values
      else:
        angle_opt.set_value(" ".join(["%g"%(p/DEG) for p in angles]));
    except:
      traceback.print_exc();
      print "Error reading %s/FEED subtable. Feed angles not filled."%Context.mssel.msname;

def _select_ms (msname):
  """This is called when a new MS is selected. Reads the FEED table if that option is enabled""";
  if read_ms:
    _read_ms_feed(True);

read_ms = False;
if Context.mssel:
  Context.mssel.when_changed(_select_ms);

#
# now define our compile-time options
#
msopt = TDLCompileOption("read_ms","Read orientation from MS FEED subtable",False,doc="""<P>Select this to automatically
  fill in feed angles from the MS FEED subtable. Some MS's do not have this filled correctly, so
  we provide options for manually setting the feed orientation here.
  </P>""");
TDLCompileOption("enable_pa","Add parallactic angle (for alt-az mounts)",False,doc="""Include the (time-variable) parallactic angle &eta;. Use this option for alt-az mounts. Disable this
for equatorial mounts,""");
angle_opt = TDLCompileOption("feed_angle","Feed orientation, degrees",["0"],more=str,doc="""<P>This is the per-station feed orientation angle &rho;<sub>q</sub>. Enter a list
  of angles (one per station), separated by spaces. Enter a single value to use it for all stations (you can also enable the "Use same Jones for all stations" checkbox below.)</P>

  <P>Note that this list is absolute and not subject to antenna subset selection. If you have one feed angle for the first ten antennas, and a different angle for the next four (i.e. for
  WSRT "cross-pol" observations), then enter the first angle ten times, and the second angle four times, regardless of what antenna subset you have selected elsewhere.</P>
  """);

msopt.when_changed(_read_ms_feed);

# set a validator, so that only numeric input is accepted for the feed angles
import re
re_whitespace = re.compile("\s+");

def _validate_angle (ang):
  # this fill throw an exception if any of the angle components cannot be converted to a float
  return map(float,re_whitespace.split(ang));
angle_opt.set_validator(_validate_angle);

# helper function, convert string to float, return 0 on error
def _str_to_float (a):
  try:
    return float(a);
  except ValueError:
    return 0;


def compute_jones (Jones,stations=None,**kw):
  stations = stations or Context.array.stations();
  if enable_pa:
    radec = Context.observation.radec0();
    xyz = Context.array.xyz();

  # get the feed angles
  angles = map(_str_to_float,re_whitespace.split(feed_angle));

  for p in stations:
    # get feed angle for this antenna number. If list contains fewer angles than stations, the last angle is reused
    angle = angles[min(Context.array.number_of_station(p),len(angles)-1)]*DEG;
    Jj = Jones(p);
    ns = Jj.Subscope();

    # if parallactic angle is enabled, make nodes to add it to the feed angle
    if enable_pa:
      ns.fa0 << angle;
      ns.pa << Meq.ParAngle(radec=radec,xyz=xyz(p));
      ns.fa << ns.fa0 + ns.pa;
      if not Context.observation.circular():
        cos = ns.cos_fa << Meq.Cos(ns.fa);
        sin = ns.sin_fa << Meq.Sin(ns.fa);
      else:
        # pexp = ns.exp_pfa << Meq.Polar(1,ns.fa);
        # nexp = ns.exp_nfa << Meq.Polar(1,-ns.fa);
        pexp = ns.exp_pfa << Meq.Polar(1,-ns.fa);
        nexp = ns.exp_nfa << Meq.Polar(1,ns.fa);
    # no p.a., work out sines and cosines directly (as constants)
    else:
      if not Context.observation.circular():
        cos = math.cos(angle);
        sin = math.sin(angle);
      else:
        pexp = cmath.exp(1j*angle);
        nexp = cmath.exp(-1j*angle);
        
    # now make the rotation matrix. 'cos' and 'sin' may be nodes or constants at this point, it doesn't matter.
    if not Context.observation.circular():
      Jj << Meq.Matrix22(cos,-sin,sin,cos);
#      print "feed_angle: linear pol"
    else:
      Jj << Meq.Matrix22(pexp,0,0,nexp);
#      print "feed_angle: circular pol"

  return Jones;
