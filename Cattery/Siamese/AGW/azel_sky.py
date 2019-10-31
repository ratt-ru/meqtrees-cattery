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
"""<P>This implements a single source that sits at a fixed azimuth and elevation (rather than
a fixed position on the sky). This is typical of RFI sources.</P> 

<P>Author: T. WIllis</P>""";
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division
from Timba.TDL import *
import Meow
import math

DEG = math.pi/180.;
ARCMIN = DEG/60.;
ARCSEC = ARCMIN/60.;

def make_source (ns,name,az,el,I=1):
  """Makes a source with a fixed azimuth/elevation direction """
  srcdir = Meow.AzElDirection(ns,name,az,el);
  if source_spi is not None:
    freq0 = source_freq0*1e+6 if source_freq0 else Meow.Context.observation.freq0();
  else:
    freq0 = None;
  if source_pol:
    Q = I*source_qi;
    U = I*source_ui;
    V = I*source_vi;
  else:
    Q = U = V = None;
  if source_type == "point":
    return Meow.PointSource(ns,name,srcdir,I=I,Q=Q,U=U,V=V,spi=source_spi,freq0=freq0);
  elif source_type == "gaussian":
    s1,s2 = gauss_smaj*ARCSEC,gauss_smin*ARCSEC;
    pa = gauss_pa*DEG;
    return Meow.GaussianSource(ns,name,srcdir,I=I,Q=Q,U=U,V=V,spi=source_spi,freq0=freq0,
        lproj=s1*math.sin(pa),mproj=s1*math.cos(pa),ratio=s2/s1);
  # else fall through and return none
  return None;

DefaultFlux = "default";

# NB: use lm0=1e-20 to avoid our nasty bug when there's only a single source
# at the phase centre
def source_list (ns,basename="S",l0=None,m0=None):
  """Creates and returns selected model""";
  global az_pos, el_pos
  sources = [make_source(ns,"S",az_pos*DEG, el_pos*DEG, source_flux)]
  return [x for x in sources if x];

TDLCompileOption("source_flux","Default source I flux, Jy",
      [1e-6,1e-3,1],more=float);

srctype_opt = TDLCompileOption("source_type","Source type",["point","gaussian"]);
gauss_menu = TDLCompileMenu("Gaussian shape",
    TDLCompileOption("gauss_smaj","Major axis extent, arcsec",10,more=float),
    TDLCompileOption("gauss_smin","Minor axis extent, arcsec",10,more=float),
    TDLCompileOption("gauss_pa","Position angle, deg",0,more=float));
srctype_opt.when_changed(lambda stype:gauss_menu.show(stype=="gaussian"));

TDLCompileMenu("Polarization",
  TDLOption("source_qi","Q/I ratio",
        [0,.1],more=float),
  TDLOption("source_ui","U/I ratio",
        [0,.1],more=float),
  TDLOption("source_vi","V/I ratio",
        [0,.1],more=float),toggle='source_pol');
TDLCompileOption("az_pos","Position of interfering source (az), in degrees",
      [0],more=float);
TDLCompileOption("el_pos","Position of interfering source (el), in degrees",
      [0],more=float);
TDLCompileOption("source_spi","Spectral index",[None],more=float);
TDLCompileOption("source_freq0","Reference frequency, MHz",[None],more=float);
