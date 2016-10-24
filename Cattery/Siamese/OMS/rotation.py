# -*- coding: utf-8 -*-
"""<P>This module implements a rotation Jones term. The term is a simple rotation matrix through -&rho;, corresponding to a feed
that has been rotated by &rho;:</P>

<TABLE>
<TR><TD> |</TD><TD>cos &rho;</TD><TD>sin &rho;</TD><TD>|</TD></TR>
<TR><TD> |</TD><TD>-sin &rho;</TD><TD>cos &rho;</TD><TD>|</TD></TR>
</TABLE>

<P>...where the rotation angle &rho; = &rho;<sub>q</sub>+&eta; is a sum of the (optional) parallactic angle &eta;, and a per-station feed orientation &rho;<sub>q</sub>.</P>

<P>Note that parallactic angle is normally interpreted in the standard IAU system (X=North, Y=East) as &eta;>0 corresponding to a rotation of local zenith from North towards East by &eta;. The rotation of the EMF vector (in the local frame) is then -&eta;. The feed angle is interpreted in the same way.</P> 

<P>An angle of 0 corresponds to a N-S X feed, and an E-W Y feed in the sky coordinate frame.
</P>

<P align="right">Author: O. Smirnov &lt;<tt>smirnov@astron.nl</tt>&gt;</P>""";

__default_label__ = "P";
__default_name__  = "parallactic angle and feed rotation";

from Timba.TDL import *
from Meow import Context,MSUtils

import math
import traceback
import os.path
import re

DEG = math.pi/180;


# setup validator for feed angle option
import re
re_whitespace = re.compile("\s+");

def _validate_angle (angle_str):
  # this fill throw an exception if any of the angle components cannot be converted to a float
  angles = [];
  for a in re_whitespace.split(angle_str):
    try: angles.append(float(a));
    except: pass;
  return angles;
  angle_opt.set_validator(_validate_angle);

# helper function, convert string to float, return 0 on error
def _str_to_float (a):
  try:
    return float(a);
  except ValueError:
    return 0;

#
# now define our compile-time options
#
PA_NONE = None;
PA_NORMAL = "enabled";
PA_INVERTED = "inverted";

class Rotation (object):
  def __init__ (self,label=__default_label__,pa=True,read_ms=True,feed_angle=True):
    self.tdloption_namespace = "rot_"+label;
    self._options = [];
    # include options based on arguments
    if feed_angle and read_ms:
      self.ms_opt = TDLOption("read_ms","Read orientation from MS FEED subtable",False,doc="""<P>Select this 
        to automatically
        fill in feed angles from the MS FEED subtable. Some MS's do not have this filled correctly, so
        we provide options for manually setting the feed orientation here.
        </P>""",namespace=self);
      self._options.append(self.ms_opt);
    else:
      self.read_ms = None;
    if pa:
      self.pa_opt = TDLOption("enable_pa","Parallactic angle rotation (for alt-az mounts)",
        [PA_NONE,PA_NORMAL,PA_INVERTED],doc="""<P>Include the (time-variable) parallactic angle &eta;. 
        Disable this for equatorial mounts (or when a sky derotator is present), 
        enable for alt-az mounts.</P>
        <P>To reverse the direction of rotation, use the "%s" setting. The feed rotation angle is then
        given by &rho; = &rho;<sub>q</sub>-&eta;. This is meant as an experimental feature only, since it doesn't match standard coordinate systems.</P>"""%(PA_INVERTED),namespace=self);
      self._options.append(self.pa_opt);
    else:
      self.enable_pa = None;
    if feed_angle:
      self.angle_opt = TDLCompileOption("feed_angle","Feed orientation, degrees",["0"],more=str,doc="""<P>
        This is the
        per-station feed orientation angle &rho;<sub>q</sub>. Enter a list
        of angles (one per station), separated by spaces. Enter a single value to use it for all stations (you can also enable the "Use same Jones for all stations" checkbox below.)</P>

        <P>Note that this list is absolute and not subject to antenna subset selection. If you have one feed angle for the first ten antennas, and a different angle for the next four (i.e. for
        WSRT "cross-pol" observations), then enter the first angle ten times, and the second angle four times, regardless of what antenna subset you have selected elsewhere.</P>
        """,namespace=self);
      self._options.append(self.angle_opt);
    else:
      self.feed_angle = None;
    # read FEED table when MS has changed
    if feed_angle and read_ms:
      self.ms_opt.when_changed(self._read_ms_feed);
      if Context.mssel:
        Context.mssel.when_changed(self._select_ms);
      
  def compile_options (self):
    return self._options;
    
  def _select_ms (self,msname):
    """This is called when a new MS is selected. Reads the FEED table if that option is enabled""";
    if self.read_ms:
      self._read_ms_feed(True);

  def _read_ms_feed (self,enabled):
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
          self.angle_opt.set_value("%g"%(angles[0]/DEG));
        # else join values
        else:
          self.angle_opt.set_value(" ".join(["%g"%(p/DEG) for p in angles]));
      except:
        traceback.print_exc();
        print "Error reading %s/FEED subtable. Feed angles not filled."%Context.mssel.msname;

  def compute_jones (self,Jones,sources=None,stations=None,**kw):
    """Returns Jones matrices. Note that this can also be used as a sky-Jones,
    (sources != None).
    """;
    stations = stations or Context.array.stations();

    # get the feed angles
    angles = map(_str_to_float,re_whitespace.split(self.feed_angle)) if self.feed_angle is not None else [0];

    for p in stations:
      # get feed angle for this antenna number. If list contains fewer angles than stations, the last angle is reused
      angle_deg = angles[min(Context.array.number_of_station(p),len(angles)-1)];
      angle = angle_deg*DEG;
      Jj = Jones(p) if sources is None else Jones(sources[0],p);
      ns = Jj.Subscope();

      # if parallactic angle is enabled, make nodes to add it to the feed angle
      if self.enable_pa != PA_NONE:
        pa = Context.observation.phase_centre.pa(Context.array.xyz(p));
        ns.fa0 << angle;
        if self.enable_pa == PA_NORMAL:
          ns.fa << ns.fa0 + pa;
        else:
          ns.fa << ns.fa0 - pa;
        if Context.observation.circular():
          pexp = ns.exp_pfa << Meq.Polar(1,ns.fa);
          nexp = ns.exp_nfa << Meq.Polar(1,-ns.fa);
        else:
          cos = ns.cos_fa << Meq.Cos(ns.fa);
          sin = ns.sin_fa << Meq.Sin(ns.fa);
      # no p.a., work out sines and cosines directly (as constants)
      else:
        if Context.observation.circular():
          pexp = cmath.rect(1,angle)
          pexp = cmath.rect(1,-angle)
        else:
          # treat these cases directly, to avoid rounding errors causing unseemly near-0 terms
          if angle_deg == 90:
            cos,sin = 0,1;
          elif angle_deg == 180:
            cos,sin = -1,0;
          elif angle_deg == 270:
            cos,sin = 0,-1;
          else:
            cos,sin = math.cos(angle),math.sin(angle);

      # now make the rotation matrix. 'cos' and 'sin' may be nodes or constants at this point, it doesn't matter.
      if Context.observation.circular():
        Jj << Meq.Matrix22(pexp,0,0,nexp);
      else:
        Jj << Meq.Matrix22(cos,sin,-sin,cos);
    
    # If used as a sky-Jones, return somethig that can be qualified with a source name
    if sources:
      for src in sources[1:]:
        for p in stations:
           Jones(src,p) << Meq.Identity(Jones(sources[0],p));
        
    return Jones;
