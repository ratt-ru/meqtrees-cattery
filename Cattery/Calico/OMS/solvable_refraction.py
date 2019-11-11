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

__default_label__ = "R";
__default_name__  = "solvable differential refraction";

from Timba.TDL import *
import math
import Meow
from Meow import Context
from Meow import StdTrees,ParmGroup

# import measures tool
try:
  import pyrap_measures
  dm = pyrap_measures.measures();
  print("Imported measures tool from pyrap_measures");
except:
  try:
    import pyrap.measures
    dm = pyrap.measures.measures();
    print("Imported measures tool from pyrap.measures");
  except:
    dm = None;
    print("Cannot find pyrap.measures or pyrap_measures, measures tool not available");

TDLCompileOption('ref_source',"Reference direction for diff refraction",[None],more=str,
  doc="""<P>Specify a direction relative to which the field is refracted. This can be a source 
        name, a source number (0-based), a direction string (e.g. "J2000 30deg 45deg"), or None
        to use the phase centre. Typically, the reference direction corresponds to the brightest 
        source in your field -- selfcal will incorporate refraction towards that source into the 
        phase solutions, so that
        <I>differential</I> refraction manifests itself relative to that direction.</P>""");
TDLCompileOption('do_extinction',"Include differential extinction estimate",False,
  doc="""<P>If enabled, differential extinction will be included in the Jones matrix. 
        This will follow a 1/sin <I>h</i> law, where <i>h</i> is the refracted elevation angle.
        Extinction is an amplitude-only effect, which will be added
        to the phase-only effect of differential refraction.</P>""");
TDLCompileOption('coord_approx',"Use 2D approximation for coodinates",False,
  doc="""<P>If enabled, coordinates will be converted using a quick-and-dirty method.</P>""");

def compute_jones (Jones,sources,stations=None,inspectors=[],meqmaker=None,label='R',**kw):
  """Creates the Z Jones for ionospheric phase, given TECs (per source, 
  per station).""";
  stations = stations or Context.array.stations;
  ns = Jones.Subscope();
  
  # get reference source
  if ref_source:
    # treat as index first
    dir0 = None;
    try:
      dir0 = sources[int(ref_source)].direction; 
    except:
      pass;
    # else treat as name, find in list
    if not dir0:
      for src0 in sources:
        if src0.name == ref_source:
          dir0 = src0.direction;
          break;
    # else treat as direction string
    if not dir0:
      ff = list(ref_source.split());
      if len(ff) < 2 or len(ff) > 3:
        raise RuntimeError("invalid reference dir '%s' specified for %s-Jones"%(ref_source,label));
      global dm;
      if not dm:
        raise RuntimeError("pyrap measures module not available, cannot use direction strings for %s-Jones"%label);
      if len(ff) == 2:
        ff =  [ 'J2000' ] +ff;
      # treat as direction measure
      try:
        dmdir = dm.direction(*ff);
      except:
        raise RuntimeError("invalid reference dir '%s' specified for %s-Jones"%(ref_source,label));
      # convert to J2000 and make direction object
      dmdir = dm.measure(dmdir,'J2000');
      ra,dec = dm.getvalue(dmdir)[0].get_value(),dm.getvalue(dmdir)[1].get_value();
      dir0 = Meow.Direction(ns,"refdir",ra,dec,static=True);
  else:
    dir0 = Context.observation.phase_centre;
 
  # make refraction scale node
  scale = ns.scale(0) << Meq.Parm(0,tags="refraction");

  xyz0 = Context.array.xyz0();
  if coord_approx:
    # get PA, and assume it's the same over the whole field
    pa = ns.pa0 << Meq.ParAngle(dir0.radec(),xyz0);
    # second column of the Rot(-PA) matrix. Multiply this by del to get a rotation of (0,del) into the lm plane.
    # The third component (0) is for convenience, as it immediately gives us dl,dm,dn, since we assume dn~0
    rot_pa = ns.rotpa0 << Meq.Composer(Meq.Sin(pa),Meq.Cos(pa),0);

  # el0: elevation of field centre
  el0 = dir0.el();
  if do_extinction:
    ns.inv_ext0 << Meq.Sin(el0);  # inverse of extinction towards el0
  # station UVWs
  uvw = Context.array.uvw();
  # now loop over sources
  for isrc,src in enumerate(sources):
    # reference direction: no refraction at all
    if src.direction is dir0:
      for p in stations:
        Jones(src,p) << 1;
      continue;
    # dEl is source elevation minus el0
    # ddEl = scale*dEl: amount by which source refracts (negative means field is compressed)
    el = src.direction.el()
    ns.dEl(src) << el - el0;
    ddel = ns.ddEl(src) << ns.dEl(src)*scale;
    # get el1: refracted elevation angle
    if not coord_approx or do_extinction:
      el1 = ns.el1(src) << el + ddel;
    # compute extinction component
    if do_extinction:
      # compute inverse of extinction towards the refracted direction el1
      iext = ns.inv_ext(src) << Meq.Sin(el1);  # 
      # and differential extinction is then ext1/ext0
      ext = ns.dext(src) << ns.inv_ext0/iext;
    # Compute dlmn offset in lm plane.
    if coord_approx:
      # Approximate mode: ddel is added to elevation, so to get the lm offset, we need
      # to apply Rot(PA) to the column vector (0,ddel), and then take the sine of the result.
      dlmn = ns.dlmn(src) << Meq.Sin(ddel*rot_pa);
    else:
      ns.azel1(src) << Meq.Composer(src.direction.az(),el1);
      ns.radec1(src) << Meq.RADec(ns.azel1(src),xyz0);
      ns.lmn1(src) << Meq.LMN(Context.observation.radec0(),ns.radec1(src));
      dlmn = ns.dlmn(src) << ns.lmn1(src) - src.lmn();
    # get per-station phases
    for p in stations:
      if do_extinction:
        Jones(src,p) << ext*(ns.phase(src,p) << Meq.VisPhaseShift(lmn=dlmn,uvw=uvw(p)));
      else:
        Jones(src,p) << Meq.VisPhaseShift(lmn=dlmn,uvw=uvw(p));
  # make bookmarks
  srcnames = [ src.name for src in sources ];
  meqmaker.make_bookmark_set(Jones,[ (src,p) for src in srcnames for p in stations ],
      "%s: inspector plot"%label,"%s: by source-station"%label,freqmean=True);
  inspectors.append(ns.inspector(label,'scale') << \
      StdTrees.define_inspector(ns.scale,[0],label=label));
  inspectors.append(ns.inspector(label,'delta-el') << \
      StdTrees.define_inspector(ns.ddEl,srcnames,label=label));
  inspectors.append(ns.inspector(label,'delta-el') << \
      StdTrees.define_inspector(ns.ddEl,srcnames,label=label));
  inspectors.append(ns.inspector(label,'dlmn') << \
      StdTrees.define_inspector(ns.dlmn,srcnames,label=label));
  if do_extinction:
    inspectors.append(ns.inspector(label,'inv-ext') << \
        StdTrees.define_inspector(ns.inv_ext,srcnames,label=label));
    inspectors.append(ns.inspector(label,'diff-ext') << \
        StdTrees.define_inspector(ns.dext,srcnames,label=label));

  # make parmgroups and solvejobs
  global pg;
  pg  = ParmGroup.ParmGroup(label,[scale],table_name="%s.fmep"%label,bookmark=False);

  # make solvejobs
  ParmGroup.SolveJob("cal_"+label,"Calibrate %s (differential refraction)"%label,pg);

  return Jones;

