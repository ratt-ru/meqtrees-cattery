# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

"""<P>This module implements a solvable differential refraction Jones matrix. Differential refraction is
computed by "compressing" the field-of-view in elevation by a solvable factor, relative to the phase
center or a designated central source.</P>

<P>NB: the current version ignores the <I>w(n-1)</I> term when computing the phase offsets due to differential
refraction. This is probably OK, as the position offsets (and thus changes in <I>n</I>) are small.</P>

<P align="right">Author: O. Smirnov &lt;<tt>smirnov@astron.nl</tt>&gt;</P>""";

__default_label__ = "Z";
__default_name__  = "simple gradient phase screen MIM";

from Timba.TDL import *
import math
import Meow
from Meow import Context
from Meow import StdTrees,ParmGroup

COMMON_GRAD = "common_gradient";
LOCAL_GRAD = "local_gradient";

TDLCompileOption("mim_type","MIM type",(
  (COMMON_GRAD,"global gradient"),
  (LOCAL_GRAD,"per-station gradients"),
));

def compose_ab (a,b,dims):
  if dims[-1] == 2:
    return Meq.Composer(a,b,dims=dims);
  else:
    return Meq.Composer(a,b,0,dims=dims);

def make_mim_parms (ns,stations,dims=[1,2]):
  """Makes gradient parameters within the given nodescope. Returns list
  of parms, and initializes ns.ab(p) for all p in stations.
  ns.ab(p) is 2/3 vector (of dimensions dims) of gradients in l and m.
  If this is a 3-vector, there's a 0 in the last position -- this
  simplifies multiplication by lmn in tensor-mode operations.""";
  if mim_type is COMMON_GRAD:
    # make common gradient parms
    ns.ab << compose_ab(
      ns.alpha << Meq.Parm(0,tags="mim"),
      ns.beta << Meq.Parm(0,tags="mim"),
      dims=dims);
    parmlist = [ns.alpha,ns.beta];
    # make per-station gradients
    xyz0 = Context.array.xyz0();
    for p in stations:
      ns.dx(p) << Meq.Norm(Context.array.xyz(p)-xyz0);
      ns.ab(p) << ns.dx(p)*ns.ab;
  else:
    # make local gradient parms
    parmlist = [];
    for p in stations:
      ns.ab(p) << compose_ab(
        ns.alpha(p) << Meq.Parm(0,tags="mim"),
        ns.beta(p) << Meq.Parm(0,tags="mim"),
        dims=dims);
      parmlist += [ns.alpha(p),ns.beta(p)];
  return parmlist;


def compute_jones (Jones,sources,stations=None,inspectors=[],meqmaker=None,label='Z',**kw):
  """Creates the Z Jones for ionospheric phase, given a TEC gradient.""";
  stations = stations or Context.array.stations;
  ns = Jones.Subscope();

  parmlist = make_mim_parms(ns,stations,dims=[1,2]);

  # choose frequency scaling so that a value of "1" produces 1 radian of phase
  # delay at 1 GHz.
  ns.freqscale = 1e+9/Meq.Freq();
  for p in stations:
    for isrc,src in enumerate(sources):
      dp = ns.dphase(src,p) << Meq.MatrixMultiply(ns.ab(p),src.direction.lm());
      Jones(src,p) << Meq.Polar(1,dp*ns.freqscale);

  # make bookmarks
  meqmaker.make_per_station_bookmarks(ns.ab,"MIM gradients",freqmean=False);
  meqmaker.make_per_source_per_station_bookmarks(ns.dphase,"MIM phase screen",
    sources,freqmean=False);

  # make solvejobs
  global pg;
  pg  = ParmGroup.ParmGroup(label,parmlist,table_name="%s.fmep"%label,bookmark=True);
  ParmGroup.SolveJob("cal_"+label,"Calibrate %s (gradient MIM)"%label,pg);

  return Jones;

def compute_jones_tensor (Jones,sources,stations=None,label="Z",lmn=None,meqmaker=None,inspectors=[],**kw):
  stations = stations or Context.array.stations;
  ns = Jones.Subscope();

  # if lmn tensor is not set for us, create a composer
  if lmn is None:
    lmn = ns.lmnT << Meq.Composer(dims=[0],*[ src.direction.lmn() for src in sources ]);

  parmlist = make_mim_parms(ns,stations,dims=[3]);

  # choose frequency scaling so that a value of "1" produces 1 radian of phase
  # delay at 1 GHz.
  ns.freqscale = 1e+9/Meq.Freq();
  xyz = Context.array.xyz;
  xyz0 = Context.array.xyz0();

  for p in stations:
    dp = ns.dphase(p) << Meq.MatrixMultiply(lmn,ns.ab(p));
    Jones(p) << Meq.Polar(1,dp*ns.freqscale);

  # make solvejobs
  global pg;
  pg  = ParmGroup.ParmGroup(label,parmlist,table_name="%s.fmep"%label,bookmark=True);
  ParmGroup.SolveJob("cal_"+label,"Calibrate %s (gradient MIM)"%label,pg);

  # make an inspector for the a and b parms
  meqmaker.make_per_station_bookmarks(ns.ab,"MIM gradients",freqmean=False);

  return Jones;
