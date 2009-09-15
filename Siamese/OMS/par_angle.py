from Timba.TDL import *
from Meow import Context

import math

def compute_jones (Jones,stations=None,**kw):
  stations = stations or Context.array.stations();
  radec = Context.observation.radec0();
  xyz = Context.array.xyz();
  
  for p in stations:
    Jj = Jones(p);
    ns = Jj.Subscope();
    ns.pa << Meq.ParAngle(radec=radec,xyz=xyz(p));
    ns.cos_pa << Meq.Cos(ns.pa);
    ns.sin_pa << Meq.Sin(ns.pa);
    Jj << Meq.Matrix22(ns.cos_pa,-ns.sin_pa,ns.sin_pa,ns.cos_pa);
    
  return Jones;
