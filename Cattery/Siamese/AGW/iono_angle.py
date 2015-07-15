# -*- coding: utf-8 -*-
__default_label__ = "F";
__default_name__  = "ionospheric Faraday correction";

from Timba.TDL import *
from Meow import Context,MSUtils

import math
import cmath
import traceback
import os.path

from Timba.TDL import *
from Timba.Meq import meq
from Timba import pynode
from Timba import mequtils
from Timba.array import *
import os
import random

import Meow

#
# This routine uses the ASCII output file produced by the ALBUS ionsphere
# analysis software to either predict or correct for ionospheric Faraday
# rotation

# Note: for obscure reasons that only Oleg understands the file name
# provided in response to the following question must contain only lower
# case characters - e.g albus_filename.iono and not ALBUS_filename.iono
TDLCompileOption("iono_file","File with RMs",TDLFileSelect("*.iono"))

def compute_jones (Jones,sources,stations=None,**kw):
  # set up rotation matrix to predict Faraday rotation effects on an
  # incoming radio wave travelling through the ionosphere to a ground-based
  # radio telescope
  Jj = Jones();
  ns = Jj.Subscope();

  # get Faraday rotation angle
  Ia = ns.rm_angle << Meq.PyNode(class_name="PyGetRMAngle",module_name="PYRMAngle", albus_filename=iono_file)
  cosia = ns << Meq.Cos(Ia); 
  sinia = ns << Meq.Sin(Ia);

  # Make the ionosphere Prediction rotation matrix.
  # When correcting obseved visibilities, the correction matrix
  # is the inverse of the following (ie the angle is the opposite of that
  # given above).
  if not Context.observation.circular():
    Jj << Meq.Matrix22(cosia,Meq.Negate(sinia),sinia,cosia);
  else:
    Jj << Meq.Matrix22(cosia + sinia*1j,0.0,0.0,cosia - sinia*1j);

  for p in stations:
    for src in sources:
      Jj(src,p) << Meq.Identity(Jj)

  return Jones;
