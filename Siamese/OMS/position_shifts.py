# -*- coding: utf-8 -*-
"""<P>This module implements a variety of per-source position shifts. You can choose to give sources
a fixed dl,dm offset, you can rotate the field (relative to l=m=0), or you can simulate differential
refraction, byt offseting each source's elevation by h*(el-el0), where el0 is elevation of the phase
centre.</P>

<P align="right">Author: O. Smirnov &lt;<tt>smirnov@astron.nl</tt>&gt;</P>""";

__default_label__ = "R";
__default_name__  = "position shifts";

from Timba.TDL import *
import math
import Meow
from Meow import Context
from Meow import StdTrees,ParmGroup

DEG = math.pi/180;
ARCMIN = DEG/60;
ARCSEC = ARCMIN/60;

LM_OFFSET = "lm offset";
DEL_OFFSET = "diff. elevation offset";
FIELD_ROT = "field rotation";

TDLCompileMenu("Type of shift",
  TDLMenu(LM_OFFSET,
    TDLOption("dl_arcsec","l, offset, arcsec",0,more=float),
    TDLOption("dm_arcsec","m, offset, arcsec",0,more=float),toggle='shift_lm_offset',
  ),
  TDLMenu(DEL_OFFSET,
    TDLOption("del_rate","Diffraction ratio",0,more=float),toggle='shift_del_offset',
  ),
  TDLMenu(FIELD_ROT,
    TDLOption("rot_deg","Rotation angle, deg",0,more=float),toggle='shift_field_rot',
  ),
  exclusive='shift_type'
);

TDLCompileOption('vis_phase_slope',"Visualize phase slope",True);

ANTX_dict = {
    '0':   0.0,        '1':  143.98881006, '2':  287.98154006, '3':   431.97187006,
    '4': 575.96470004, '5':  719.95646011, '6':  863.94757006, '7':  1007.93746007,
    '8':1151.92894011, '9': 1295.9213701 , 'A': 1331.92456019, 'B':  1403.91958018,
    'C':2627.84690046, 'D': 2699.84118052
};

def compute_jones (Jones,sources,stations=None,inspectors=[],meqmaker=None,label='R',**kw):
  """Creates the Z Jones for ionospheric phase, given TECs (per source, 
  per station).""";
  stations = stations or Context.array.stations;
  ns = Jones.Subscope();
  dir0 = Context.observation.phase_centre;
  if shift_del_offset:
    el0 = dir0.el();
    zyx0 = Context.array.xyz0();
  elif shift_field_rot:
    cos,sin = math.cos(rot_deg*DEG),math.sin(rot_deg*DEG);
    ns.RM << Meq.Matrix22(cos,-sin,sin,cos);

  uvw = Context.array.uvw();

  # now loop over sources
  for isrc,src in enumerate(sources):
    # work out dl,dm offset and dlmn node
    print shift_type;
    if shift_lm_offset:
      dl = ns.dl(src) << dl_arcsec*ARCSEC;
      dm = ns.dm(src) << dm_arcsec*ARCSEC;
      dlmn = ns.dlmn(src) << Meq.Composer(dl,dm,0);
    elif shift_del_offset:
      el = src.direction.el()
      ns.dEl(src) << el - el0;
      ddel = ns.ddEl(src) << ns.dEl(src)*del_rate;
      el1 = ns.el1(src) << el + ddel;
      ns.azel1(src) << Meq.Composer(src.direction.az(),el1);
      ns.radec1(src) << Meq.RADec(ns.azel1(src),Context.array.xyz0());
      ns.lmn1(src) << Meq.LMN(Context.observation.radec0(),ns.radec1(src));
      dlmn = ns.dlmn(src) << ns.lmn1(src) - src.lmn();
    elif shift_field_rot:
      lm1 = ns.lm1(src) << Meq.MatrixMultiply(ns.RM,src.direction.lm());
      dlmn = ns.dlmn(src) << Meq.Composer(lm1-src.direction.lm(),0);
    # work out corresponding phase shift
    for p in stations:
      Jones(src,p) << Meq.VisPhaseShift(lmn=dlmn,uvw=uvw(p));
    # compute phase slope
    if vis_phase_slope:
      Jones(src,'slope') << (Meq.Arg(Jones(src,stations[-1])) - Meq.Arg(Jones(src,stations[0])))/(ANTX_dict[stations[-1]]*DEG);
 
    
  # make bookmarks
  srcnames = [ src.name for src in sources ];
  meqmaker.make_bookmark_set(Jones,[ (src,p) for src in srcnames for p in stations ],
      "%s: inspector plot"%label,"%s: by source-station"%label,freqmean=True);
  meqmaker.make_bookmark_set(Jones,[ (src,stations[-1]) for src in srcnames ],
      "%s: inspector plot (last station)"%label,"%s: last station"%label,inspector_node=Jones('max'),freqmean=True);
  if vis_phase_slope:
    meqmaker.make_bookmark_set(Jones,[ (src,'slope') for src in srcnames ],
        "%s: phase slopes"%label,"%s: phase slopes"%label,inspector_node=Jones('slope'),freqmean=True);

  return Jones;

