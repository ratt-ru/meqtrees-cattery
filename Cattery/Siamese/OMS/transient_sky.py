from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from Timba.TDL import *
import Meow
import math

DEG = math.pi/180.;
ARCMIN = DEG/60;

def transient_source (ns,name,l,m,tburst,duration,Ipeak=1):
  """shortcut for making a transient pointsource with a direction. Returns None for sources out of the "sky"
  (l^2+m^2>1)""";
  l = math.sin(l);
  m = math.sin(m);
  t_rel = ns.t_rel ** (Meq.Time() - (ns.time0 << 0))
  lc = ns.lc << Meq.Exp((-Meq.Pow(t_rel-tburst,2))/(2*duration**2))
  Ilc = ns.Ilc << Ipeak * lc
  if l*l + m*m <= 1:
    srcdir = Meow.LMDirection(ns,name,l,m);
    return Meow.PointSource(ns,name,srcdir,I=Ilc);
  else:
    return None;

def transient_model (ns,basename,l0,m0,dl,dm,nsrc,tburst,duration,I):
  model = [transient_source(ns,basename,l0+dl,m0+dm,tburst,duration,I)]
  return model

def cross_model (ns,basename,l0,m0,dl,dm,nsrc,tburst,duration,I):
  """Returns sources arranged in a cross""";
  model = [transient_source(ns,basename+"+0+0",l0,m0,tburst,duration,I)];
  dy = 0;
  nsrc = int(nsrc)
  for dx in range(-nsrc,nsrc+1):
    if dx:
      name = "%s%+d%+d" % (basename,dx,dy);
      model.append(transient_source(ns,name,l0+dl*dx,m0+dm*dy,tburst,duration,I));
  dx = 0;
  for dy in range(-nsrc,nsrc+1):
    if dy:
      name = "%s%+d%+d" % (basename,dx,dy);
      model.append(transient_source(ns,name,l0+dl*dx,m0+dm*dy,tburst,duration,I));
  return model;

# NB: use lm0=1e-20 to avoid our nasty bug when there's only a single source
# at the phase centre
def source_list (ns,basename="S",l0=0,m0=0):
  """Creates and returns selected model""";
  if grid_size == 1 and not l0 and not m0:
    l0 = m0 = 1e-20;
  sources = model_func(ns,basename,l0,m0,
                       grid_step*ARCMIN,grid_step*ARCMIN,
                       (grid_size-1)/2,tburst,duration,source_flux);
  return [x for x in sources if x];

# model options
model_option = TDLCompileOption("model_func","Sky model type",
      [cross_model,transient_model]);

TDLCompileOption("tburst","Relative time of burst in seconds",
      [3600,7200,14400,28800,57600],more=int);
TDLCompileOption("duration","Duration of burst in seconds",
      [3600,7200,14400],more=int);
TDLCompileOption("grid_size","Number of sources in each direction",
      [1,3,5,7],more=int);
TDLCompileOption("grid_step","Grid step or offset, in arcmin",
      [.1,.5,1,2,5,10,15,20,30],more=float);
TDLCompileOption("source_flux","Source flux, Jy",
      [1e-6,1e-3,1],more=float);

