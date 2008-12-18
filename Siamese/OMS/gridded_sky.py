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

DEG = math.pi/180.;
ARCMIN = DEG/60;

def estimate_image_size (**kw):
  """Returns current image size, based on grid size and step""";
  return grid_size*grid_step;

def point_source (ns,name,l,m,I=1):
  """shortcut for making a pointsource with a direction. Returns None for sources out of the "sky"
  (l^2+m^2>1)""";
  l = math.sin(l);
  m = math.sin(m);
  if l*l + m*m <= 1:
    srcdir = Meow.LMDirection(ns,name,l,m);
    if source_spi is not None:
      freq0 = source_freq0 or Meow.Context.observation.freq0();
    else:
      freq0 = None;
    if source_pol:
      Q = I*source_qi;
      U = I*source_ui;
      V = I*source_vi;
    else:
      Q = U = V = None;
    return Meow.PointSource(ns,name,srcdir,I=I,Q=Q,U=U,V=V,spi=source_spi,freq0=freq0);
  else:
    return None;

def cross_model (ns,basename,l0,m0,dl,dm,nsrc,I):
  """Returns sources arranged in a cross""";
  model = [point_source(ns,basename+"+0+0",l0,m0,I)];
  dy = 0;
  for dx in range(-nsrc,nsrc+1):
    if dx:
      name = "%s%+d%+d" % (basename,dx,dy);
      model.append(point_source(ns,name,l0+dl*dx,m0+dm*dy,I));
  dx = 0;
  for dy in range(-nsrc,nsrc+1):
    if dy:
      name = "%s%+d%+d" % (basename,dx,dy);
      model.append(point_source(ns,name,l0+dl*dx,m0+dm*dy,I));
  return model;

def mbar_model (ns,basename,l0,m0,dl,dm,nsrc,I):
  """Returns sources arranged in a line along the m axis""";
  model = [];
  model.append(point_source(ns,basename+"+0",l0,m0,I));
  for dy in range(-nsrc,nsrc+1):
    if dy:
      name = "%s%+d" % (basename,dy);
      model.append(point_source(ns,name,l0,m0+dm*dy,I));
  return model;

def lbar_model (ns,basename,l0,m0,dl,dm,nsrc,I):
  """Returns sources arranged in a line along the m axis""";
  model = [];
  model.append(point_source(ns,basename+"+0",l0,m0,I));
  for dx in range(-nsrc,nsrc+1):
    if dx:
      name = "%s%+d" % (basename,dx);
      model.append(point_source(ns,name,l0+dl*dx,m0,I));
  return model;
  
def star8_model (ns,basename,l0,m0,dl,dm,nsrc,I):
  """Returns sources arranged in an 8-armed star""";
  model = [ point_source(ns,basename+"+0+0",l0,m0,I) ];
  for n in range(1,nsrc+1):
    for dx in (-n,0,n):
      for dy in (-n,0,n):
        if dx or dy:
          name = "%s%+d%+d" % (basename,dx,dy);
          model.append(point_source(ns,name,l0+dl*dx,m0+dm*dy,I));
  return model;

def grid_model (ns,basename,l0,m0,dl,dm,nsrc,I):
  """Returns sources arranged in a grid""";
  model = [ point_source(ns,basename+"+0+0",l0,m0,I) ];
  for dx in range(-nsrc,nsrc+1):
    for dy in range(-nsrc,nsrc+1):
      if dx or dy:
        name = "%s%+d%+d" % (basename,dx,dy);
        model.append(point_source(ns,name,l0+dl*dx,m0+dm*dy,I));
  return model;

def circ_grid_model (ns,basename,l0,m0,dl,dm,nsrc,I):
  """Returns sources arranged in a circular grid""";
  # start with a cross model
  model = cross_model(ns,basename,l0,m0,dl,dm,nsrc,I);
  # fill in diagonals
  dl /= math.sqrt(2);
  dm /= math.sqrt(2);
  for n in range(1,nsrc+1):
    for dx in (-n,0,n):
      for dy in (-n,0,n):
        if dx and dy:
          name = "%s%+d%+d" % (basename,dx,dy);
          model.append(point_source(ns,name,l0+dl*dx,m0+dm*dy,I));
  return model;

# NB: use lm0=1e-20 to avoid our nasty bug when there's only a single source
# at the phase centre
def source_list (ns,basename="S",l0=None,m0=None):
  """Creates and returns selected model""";
  l0 = l0 or grid_l0*ARCMIN;
  m0 = m0 or grid_m0*ARCMIN;
  if grid_size == 1 and not l0 and not m0:
    l0 = m0 = 1e-20;
  sources = model_func(ns,basename,l0,m0,
                       grid_step*ARCMIN,grid_step*ARCMIN,
                       (grid_size-1)/2,source_flux);
  return filter(lambda x:x,sources);

# model options
model_option = TDLCompileOption("model_func","Sky model type",
      [cross_model,grid_model,circ_grid_model,star8_model,lbar_model,mbar_model]);

TDLCompileOption("grid_size","Number of sources in each direction",
      [1,3,5,7],more=int);
TDLCompileOption("grid_step","Grid step, in arcmin",
      [.1,.5,1,2,5,10,15,20,30],more=float);
TDLCompileOption("source_flux","Source I flux, Jy",
      [1e-6,1e-3,1],more=float);
TDLCompileMenu("Polarization",
  TDLOption("source_qi","Q/I ratio",
        [0,.1],more=float),
  TDLOption("source_ui","U/I ratio",
        [0,.1],more=float),
  TDLOption("source_vi","V/I ratio",
        [0,.1],more=float),toggle='source_pol');
TDLCompileOption("grid_l0","Offset w.r.t. phase center (l), in arcmin",
      [0],more=float);
TDLCompileOption("grid_m0","Offset w.r.t. phase center (m), in arcmin",
      [0],more=float);
TDLCompileOption("source_spi","Spectral index",[None],more=float);
TDLCompileOption("source_freq0","Reference frequency, MHz",[None],more=float);

