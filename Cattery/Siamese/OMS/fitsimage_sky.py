# -*- coding: utf-8 -*-
#
#% $Id$ 
#
#
# Copyright (C) 2002-2007
# The MeqTree Foundation & 
# ASTRON (Netherlands Foundation for Research in Astronomy)
# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>,
# or write to the Free Software Foundation, Inc., 
# 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

from Timba.TDL import *
import Meow
import math
from Meow import FITSImageComponent
from Meow import Context

DEG = math.pi/180.;
ARCMIN = DEG/60;

# NB: use lm0=1e-20 to avoid our nasty bug when there's only a single source
# at the phase centre
def source_list (ns,basename="S",l0=None,m0=None):
  """Creates and returns selected model""";
  img = FITSImageComponent(ns,basename+'IMG',filename=image_filename,
                           direction=Context.observation.phase_centre);
  img.set_options(fft_pad_factor=pad_factor);
  return [ img ];


TDLCompileOption("image_filename","FITS image file",TDLFileSelect("*.fits"));
TDLCompileOption("pad_factor","Padding factor",[1,1.1,1.2,1.5,2],more=float,default=1.2,
  doc="""<P>During prediction from FITS images, flux scale may become severely distorted
  at the edges of the image. If your image is zero at the edges, this is not a problem.
  If your image contains non-zero flux close to its edges, you need to enable padding by setting the
  pad factor to >1. The image will then be padded out with zeroes out to
  a multiple of the original size (e.g. with a pad factor of 1.5, a 100x100 image
  will be padded out to a size of 150x150.)</P>

  <P>Padding increases your effective image size and therefore uses more RAM. With large images,
  don't enable it unless you have to.</P>""");

