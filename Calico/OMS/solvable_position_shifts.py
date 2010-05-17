# -*- coding: utf-8 -*-
"""<P>This module implements solvable positional offsets.</P>

<P align="right">Author: O. Smirnov &lt;<tt>smirnov@astron.nl</tt>&gt;</P>""";

__default_label__ = "R";
__default_name__  = "solvable position shifts";

from Timba.TDL import *
import math
import Meow
from Meow import Context
from Meow import StdTrees,ParmGroup

def compute_jones (Jones,sources,stations=None,inspectors=[],meqmaker=None,label='R',**kw):
  """Creates the Z Jones for ionospheric phase, given TECs (per source, 
  per station).""";
  stations = stations or Context.array.stations;
  ns = Jones.Subscope();
  
  parmdef = Meq.Parm(0,tags="pos_offset");
  parms = [];
  uvw = Context.array.uvw();
  # now loop over sources
  for isrc,src in enumerate(sources):
    parms += [ns.dl(src) << parmdef,ns.dm(src) << parmdef];
    dlmn = ns.dlmn(src) << Meq.Composer(ns.dl(src),ns.dm(src),0);
    for p in stations:
      Jones(src,p) << Meq.VisPhaseShift(lmn=dlmn,uvw=uvw(p));
  # make bookmarks
  srcnames = [ src.name for src in sources ];
  meqmaker.make_bookmark_set(Jones,[ (src,p) for src in srcnames for p in stations ],
      "%s: inspector plot"%label,"%s: by source-station"%label,freqmean=True);
  inspectors.append(ns.inspector(label,'dlmn') << \
      StdTrees.define_inspector(ns.dlmn,srcnames,label=label));

  # make parmgroups and solvejobs
  global pg;
  pg  = ParmGroup.ParmGroup(label,parms,table_name="%s.fmep"%label,bookmark=False);

  # make solvejobs
  ParmGroup.SolveJob("cal_"+label,"Calibrate %s (position shifts)"%label,pg);

  return Jones;

