# -*- coding: utf-8 -*-
"""<P>This module implements a direction-independent leakage (D-)Jones term. The term takes the form of</P>

<TABLE>
<TR><TD> |</TD><TD>1</TD><TD>d</TD><TD>|</TD></TR>
<TR><TD> |</TD><TD>--d</TD><TD>1</TD><TD>|</TD></TR>
</TABLE>

<P align="right">Author: O. Smirnov &lt;<tt>smirnov@astron.nl</tt>&gt;</P>""";

__default_label__ = "D";
__default_name__  = "direction-independent leakage";

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

class Leakage (object):
  def __init__ (self,label=__default_label__):
    self.tdloption_namespace = "leakage_"+label
    self._options = []
    # include options based on arguments
    self.d_opt = TDLOption("d","Leakage value (d)",[0,0.01,0.04],more=float,doc="""<P>The value
        of the off-diagonal leakage term</P>""",namespace=self);
    self._options.append(self.d_opt);
      
  def compile_options (self):
    return self._options;

  def compute_jones (self,Jones,stations=None,**kw):
    """Returns Jones matrices. Note that this can also be used as a sky-Jones,
    (sources != None).
    """;
    stations = stations or Context.array.stations();

    for p in stations:
      Jones(p) << Meq.Matrix22(1,self.d,-self.d,1);
        
    return Jones;
