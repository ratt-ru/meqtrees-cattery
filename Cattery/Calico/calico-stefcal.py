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

 # standard preamble
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from Timba.TDL import *
from Timba.Meq import meq
import math
import os
import os.path

import Meow
from Meow import ParmGroup,Bookmarks,StdTrees,Context

# This defines some ifr subsets that are commonly used for WSRT data,
# to be offered as defaults in the GUI wherever ifrs are selected.
# MSUtils already has this for WSRT, so reuse it here.
STD_IFR_SUBSETS = Meow.MSUtils.STD_IFR_SUBSETS["WSRT"];

TDLCompileOptionSeparator("MS selection");
# MS options first
mssel = Meow.Context.mssel = Meow.MSUtils.MSSelector(has_input=True,tile_sizes=[1,8,16,32],max_tiles=[1],
                  read_flags=True,write_flags=True,
                  hanning=True,invert_phases=True);
# MS compile-time options
TDLCompileOptions(*mssel.compile_options());
TDLCompileOption("run_purr","Start Purr on this MS",False);
# MS run-time options
TDLRuntimeMenu("Data selection & flag handling",*mssel.runtime_options());


# now load optional modules for the ME maker
from Meow import TensorMeqMaker
meqmaker = TensorMeqMaker.TensorMeqMaker(solvable=True,
                            use_jones_inspectors=True,
                            use_skyjones_visualizers=False,
                            use_decomposition=False
);

## specify available sky models
## these will show up in the menu automatically
# from Calico.OMS import central_point_source
# from Siamese.OMS import fitsimage_sky,gridded_sky
# models = [central_point_source,fitsimage_sky,gridded_sky]

models = [];
try:
  from Siamese.OMS.tigger_lsm import TiggerSkyModel
  models.insert(0,TiggerSkyModel());
except:
  traceback.print_exc();
  pass;

meqmaker.add_sky_models(models);

read_ms_model_opt = TDLCompileOption("read_ms_model","Read additional uv-model visibilities from MS",False,doc="""
  <P>If enabled, then an extra set of <i>model</i> visibilities will be read from a second column
  of the MS (in addition to the input data.) These can either be added to whatever is predicted by the
  sky model <i>in the uv-plane</i> (i.e. subject to uv-Jones but not sky-Jones corruptions), or directly
  subtracted from the input data. See also the "Equation type" option below.

  <P>If you are repeatedly running a large sky model, and are not solving for any image-plane effects, then this
  feature lets you compute your sky model (or a part of it) just once, save the result to the uv-model column,
  and reuse it in subsequent steps. This can save a lot of processing time.</P>
  """);


from Siamese.OMS.rotation import Rotation
from Siamese.OMS import oms_dipole_projection
meqmaker.add_sky_jones('L','parallactic angle or dipole rotation',[Rotation('L',feed_angle=False),oms_dipole_projection])

# E - beam
# add a fixed primary beam first
from Calico.OMS import wsrt_beams  #,wsrt_beams_zernike
from Siamese.OMS import pybeams_fits
from Siamese.OMS.emss_beams import emss_polar_beams

meqmaker.add_sky_jones('E','primary beam',[wsrt_beams,pybeams_fits,emss_polar_beams]); # ,wsrt_beams_zernike]);

# P - feed angle
from Siamese.OMS import feed_angle
meqmaker.add_uv_jones('P','feed orientation',[feed_angle]);

# very important -- insert meqmaker's options properly
TDLCompileOptions(*meqmaker.compile_options());


TDLCompileOptionSeparator("Stefcal options");

# if table access is available, add baseline selection options
if Meow.MSUtils.TABLE:
  calib_ifrs_opt = TDLCompileOption('calibrate_ifrs',"Interferometer subset for calibration",
    ["all"]+STD_IFR_SUBSETS,more=str,doc="""<P>
      You can restrict calibration to a subset of interferometers. Note that this selection
      applies on top of (not instead of) the global interferometer selection specified above.
      """+Meow.IfrArray.ifr_spec_syntax);
  mssel.when_changed(lambda msname:calib_ifrs_opt.set_option_list(mssel.ms_ifr_subsets));

CORRECTED_DATA = "CORR_DATA";
CORRECTED_DATA_SUB = "CORR_DATA_SUB";
RESIDUALS = "RES";
CORRECTED_RESIDUALS = "CORR_RES";
CORRUPTED_MODEL = "PREDICT";
CORRUPTED_MODEL_ADD = "DATA+PREDICT";

output_option = TDLCompileOption('do_output',"Output visibilities",
                                 [  (CORRECTED_DATA,"corrected data"),
                                    (CORRECTED_DATA_SUB,"corrected data minus DDS"),
                                    (CORRECTED_RESIDUALS,"corrected residuals") ],
  doc="""<P>This selects what sort of visibilities get written to the output column.</P>
  <ul>

  <li><B>Corrected data</B>: input data corrected for the direction-independent (DI) Jones terms;</li>

  <li><B>Corrected residuals</B>: corrected data minus complete sky model;</li>

  <li><B>Corrected data minus DDS</B>: corrected data minus only those sky model components that have a dE term.</li>

  </ul>

  """);


diffgain_tag = 'dE';
diffgain_group = 'cluster';

from Calico.OMS.StefCal.GainOpts import GainOpts,MODE_SOLVE_SAVE,MODE_SOLVE_NOSAVE,MODE_SOLVE_APPLY

gopts = GainOpts("direction-independent gain","gain","G","stefcal");
TDLCompileOptions(*gopts.tdl_options);
bopts = GainOpts("direction-independent gain","gain1","B","stefcal");
TDLCompileOptions(*bopts.tdl_options);
dgsel = TensorMeqMaker.SourceSubsetSelector("Apply diffgains to selected sources",
          tdloption_namespace='de_subset',annotate=False);
deopts = GainOpts("differential gain","diffgain","dE","stefcal",pre_opts=dgsel.options);
TDLCompileOptions(*deopts.tdl_options);

DIAGONLY,ALLFOUR = "diag","full";
TDLCompileMenu("Use interferometer errors",
  TDLOption("stefcal_ifr_gain_mode","Solution mode",
    {MODE_SOLVE_SAVE:"solve and save",MODE_SOLVE_NOSAVE:"solve, do not save",MODE_SOLVE_APPLY:"load and apply"}),
  TDLOption("stefcal_ifr_gain_reset","Ignore previously saved solutions",False),
  TDLOption("stefcal_diagonal_ifr_gains","Polarizations",
            {DIAGONLY:"parallel-hand only",ALLFOUR:"full 2x2"}),
  TDLOption("stefcal_per_chan_ifr_gains","Solve on a per-channel basis",False),
  TDLOption("stefcal_ifr_gain_table","Filename for solutions",["ifrgains.ma"],more=str),
  toggle="stefcal_ifr_gains",
);
TDLCompileOption("stefcal_nmajor","Number of major loops",[1,2,3,5],more=int,default=2);
TDLCompileOption("stefcal_rescale","Rescale data to model before solving",["no","scalar","per slot"]);
TDLCompileOption("stefcal_noise_per_chan","Use per-channel noise estimates",True);
stefcal_downsample = False;
#TDLCompileMenu("Use on-the-fly downsampling",
#  TDLCompileOption("stefcal_downsample_timeint","Downsampling interval, time axis (1 for full resolution)",[1],more=int,default=1),
# TDLCompileOption("stefcal_downsample_freqint","Downsampling interval, freq axis (1 for full resolution)",[1],more=int,default=1),
#  toggle="stefcal_downsample");
TDLCompileOption("critical_flag_threshold","Critical flag threshold",[10,20,50,100],more=int,default=20,
  doc=
  """If percentage of flagged data exceeds this threshold, stop with an error message. Set to 100 to disable. This
  is a sanity check meant to catch situations when incorrect StefCal settings cause all (or too much) data to be flagged.
  """
  );
TDLCompileMenu("Enable visualizers",
  TDLCompileOption("visualize_flag_unity","Flag zero/unity solutions in visualizers",True),
  TDLCompileOption("visualize_norm_offdiag","Normalize off-diagonal terms by diagonals",True),
  toggle="stefcal_visualize");
TDLCompileOption("stefcal_verbose","Stefcal verbosity level",[0,1,2,3],more=int);

import Purr.Pipe

def _define_forest(ns,parent=None,**kw):
  if not mssel.msname:
    raise RuntimeError("MS not set");
  if run_purr:
    Timba.TDL.GUI.purr(mssel.msname+".purrlog",[mssel.msname,'.']);
  # create Purr pipe
  global purrpipe;
  purrpipe = Purr.Pipe.Pipe(mssel.msname);

  # setup contexts from MS
  mssel.setup_observation_context(ns);
  array = Meow.Context.array;

  # make spigot nodes for data
  mssel.enable_input_column(True);
  spigots = array.spigots(corr=mssel.get_corr_index());
  meqmaker.make_per_ifr_bookmarks(spigots,"Input visibilities");

  # data tensor
  ns.DT << Meq.Composer(dims=[0],mt_polling=True,*[ spigots(p,q) for p,q in array.ifrs() ]);
  # list of model tensors
  models = [];

  if read_ms_model:
    mssel.enable_model_column(True);
    model_spigots = array.spigots(column="PREDICT",corr=mssel.get_corr_index());
    meqmaker.make_per_ifr_bookmarks(model_spigots,"UV-model visibilities");
    mtuv = ns.MT("uv") << Meq.Composer(dims=[0],mt_polling=True,*[ model_spigots(p,q) for p,q in array.ifrs() ]);

  # predict tree using the MeqMaker
  all_sources = meqmaker.get_source_list(ns);
  dg_sources = deopts.enabled and dgsel.filter(all_sources);
  if dg_sources:
    # first group all sources without a diffgain on them
    groups = [ [ src for src in all_sources if not src in dg_sources ] ];
    # group diffgain-enabled sources by grouping tag
    clusters = set([src.get_attr(diffgain_group,None) for src in dg_sources]);
    dg_groups = [ (name,[ src for src in dg_sources if src.get_attr(diffgain_group) == name ]) for name in clusters if name ];
    # add sources without a grouping tag individually, as single-source groups
    dg_groups += [ (src.name,[src]) for src in dg_sources if not src.get_attr(diffgain_group,None) ];
    # now sort by brightness
    flux_dgg = [ (sum([src.get_attr('Iapp',0) or src.get_attr('I') for src in dgg[1]]),dgg) for dgg in dg_groups ];
    from past.builtins import cmp
    from functools import cmp_to_key
    flux_dgg.sort(key=cmp_to_key(lambda a,b:cmp(b[0],a[0])));
    diffgain_labels = [ dgg[0] for flux,dgg in flux_dgg ];
    groups += [ dgg[1] for flux,dgg in flux_dgg ];
    num_diffgains = len(flux_dgg);
    # now make predict trees
    for i,group in enumerate(groups):
      MT = ns.MT(group[0].name if i else "all");
      # first tensor is "MT", rest qualified by source names
      predict = meqmaker.make_predict_tree(MT.Subscope(),sources=group);
      MT << Meq.Composer(dims=[0],mt_polling=True,*[ predict(p,q) for p,q in array.ifrs() ]);
      # if reading an extra uv-model, add to first term
      if not i and read_ms_model:
        MT = ns.MT << Meq.Add(MT,mtuv)
      models.append(MT);
    print("Number of diffgain predict groups:",len(groups));
  else:
    diffgain_labels = [];
    num_diffgains = 0;
    predict = meqmaker.make_predict_tree(ns);
    MT = ns.MT("all") << Meq.Composer(dims=[0],mt_polling=True,*[ predict(p,q) for p,q in array.ifrs() ]);
    if read_ms_model:
      MT = ns.MT << Meq.Add(MT,mtuv)
    models.append(MT)
    
  solve_ifrs  = array.subset(calibrate_ifrs,strict=False).ifrs();
  downsample_subtiling = [ stefcal_downsample_timeint,stefcal_downsample_freqint ] if stefcal_downsample else [1,1];

  import Calico.OMS.StefCal.StefCal
  kwopts = {}
  gopts.set_stefcal_node_options(kwopts,visualize=stefcal_visualize);
  bopts.set_stefcal_node_options(kwopts,visualize=stefcal_visualize);
  deopts.set_stefcal_node_options(kwopts,visualize=stefcal_visualize);
  ns.stefcal << Meq.PyNode(class_name="StefCalNode",module_name=Calico.OMS.StefCal.StefCal.__file__,
                           ifrs=[ "%s:%s"%(p,q) for p,q in array.ifrs() ],
                           baselines=[ array.baseline(ip,iq) for (ip,p),(iq,q) in array.ifr_index() ],
                           solve_ifrs=[ "%s:%s"%(p,q) for p,q in solve_ifrs ],
                           noise_per_chan=stefcal_noise_per_chan,
                           downsample_subtiling=downsample_subtiling,
                           num_major_loops=stefcal_nmajor,
                           regularization_factor=1e-6,#
                           rescale=stefcal_rescale,
                           init_from_previous=False,
                           critical_flag_threshold=critical_flag_threshold,
                           diffgain_labels=diffgain_labels,
                           # flagging options
                           output_flag_bit=Meow.MSUtils.FLAGMASK_OUTPUT,
                           # IFR gain solution options
                           apply_ifr_gains=stefcal_ifr_gains,
                           solve_ifr_gains=(stefcal_ifr_gain_mode != MODE_SOLVE_APPLY),
                           reset_ifr_gains=stefcal_ifr_gain_reset,
                           save_ifr_gains=(stefcal_ifr_gain_mode == MODE_SOLVE_SAVE),
                           ifr_gain_table=stefcal_ifr_gain_table,
                           per_chan_ifr_gains=stefcal_per_chan_ifr_gains,
                           diag_ifr_gains=(stefcal_diagonal_ifr_gains == DIAGONLY),
                           residuals=(do_output == CORRECTED_RESIDUALS),
                           subtract_dgsrc=(do_output == CORRECTED_DATA_SUB),
                           verbose=stefcal_verbose,
                           children=[ns.DT]+models,**kwopts);
                           
  inspectors = meqmaker.get_inspectors() or [];
  # make output bookmarks
  nv = 0;
  for p,q in array.ifrs():
    sel = ns.output_sel(p,q) << Meq.Selector(ns.stefcal,index=list(range(nv,nv+4)),multi=True);
    ns.output(p,q) << Meq.Composer(sel,dims=[2,2]);
    nv += 4;
  meqmaker.make_per_ifr_bookmarks(ns.output,"Output visibilities");
  
  Bookmarks.Page("StefCal outputs").add(ns.stefcal,viewer="Record Browser");

  if gopts.enabled and gopts.visualize and stefcal_visualize:
    ns.stefcal_vis_G << Meq.PyNode(class_name="StefCalVisualizer",module_name=Calico.OMS.StefCal.StefCal.__file__,
      label="G",flag_unity=visualize_flag_unity,norm_offdiag=visualize_norm_offdiag,
      vells_label=Context.correlations);
    ns.stefcal_vis_G_avg << Meq.PyNode(class_name="StefCalVisualizer",module_name=Calico.OMS.StefCal.StefCal.__file__,
      label="G",freq_average=True,flag_unity=visualize_flag_unity,norm_offdiag=visualize_norm_offdiag,
      vells_label=Context.correlations);
    Bookmarks.Page("StefCal G plotter").add(ns.stefcal_vis_G,viewer="Result Plotter");
    Bookmarks.Page("StefCal G inspector").add(ns.stefcal_vis_G_avg,viewer="Collections Plotter");
    inspectors += [ ns.stefcal_vis_G,ns.stefcal_vis_G_avg ];
  if bopts.enabled and bopts.visualize and stefcal_visualize:
    ns.stefcal_vis_B << Meq.PyNode(class_name="StefCalVisualizer",module_name=Calico.OMS.StefCal.StefCal.__file__,
      label="B",flag_unity=visualize_flag_unity,norm_offdiag=visualize_norm_offdiag,
      vells_label=Context.correlations);
    ns.stefcal_vis_B_avg << Meq.PyNode(class_name="StefCalVisualizer",module_name=Calico.OMS.StefCal.StefCal.__file__,
      label="B",freq_average=True,flag_unity=visualize_flag_unity,norm_offdiag=visualize_norm_offdiag,
      vells_label=Context.correlations);
    Bookmarks.Page("StefCal B plotter").add(ns.stefcal_vis_B,viewer="Result Plotter");
    Bookmarks.Page("StefCal B inspector").add(ns.stefcal_vis_B_avg,viewer="Collections Plotter");
    inspectors += [ ns.stefcal_vis_B,ns.stefcal_vis_B_avg ];
  if deopts.enabled and deopts.visualize and stefcal_visualize:
    for i,label in enumerate(diffgain_labels):
      vde = ns.stefcal_vis_dE(label) << Meq.PyNode(class_name="StefCalVisualizer",module_name=Calico.OMS.StefCal.StefCal.__file__,
        label="dE:%s"%label,flag_unity=visualize_flag_unity,norm_offdiag=visualize_norm_offdiag,
        vells_label=Context.correlations);
      vde_avg = ns.stefcal_vis_dE_avg(label) << Meq.PyNode(class_name="StefCalVisualizer",module_name=Calico.OMS.StefCal.StefCal.__file__,
                                    label="dE:%s"%label,freq_average=True,flag_unity=visualize_flag_unity,norm_offdiag=visualize_norm_offdiag,
                                    vells_label=Context.correlations);
      Bookmarks.Page("StefCal dE:%s plotter"%label).add(vde,viewer="Result Plotter");
      Bookmarks.Page("StefCal dE:%s inspector"%label).add(vde_avg,viewer="Collections Plotter");
      inspectors += [ vde,vde_avg ];

  # make sinks
  StdTrees.make_sinks(ns,ns.output,spigots=spigots,post=inspectors,
      corr_index=mssel.get_corr_index());
  # this should come last, since runtime options may be built up during compilation.
  TDLRuntimeOptions(*meqmaker.runtime_options(nest=False));

  # finally, setup imaging options
  imsel = mssel.imaging_selector(npix=512,arcmin=meqmaker.estimate_image_size());
  TDLRuntimeMenu("Make an image from this MS",*imsel.option_list());

  # and close meqmaker -- this exports annotations, etc
  meqmaker.close();
  
  # add options to clear all solutions 
  from Calico.OMS.StefCal import StefCal
  TDLRuntimeOption("stefcal_reset_all","Remove all existing solutions",False);
  for opt in gopts,bopts,deopts:
    if opt.enabled:
      TDLRuntimeOption("reset","Remove existing %s solutions (%s)"%(opt.label,os.path.basename(opt.table)),False,namespace=opt);
  if stefcal_ifr_gains:
    TDLRuntimeOption("stefcal_reset_ifr_gains","Remove existing interferometer errors (%s)"%(
        os.path.basename(stefcal_ifr_gain_table)),False);
  TDLRuntimeJob(_run_stefcal,"Run StefCal",job_id="stefcal");

def _run_stefcal (mqs,parent,wait=False):
  # collect which tables to remove
  to_remove = [];
  for opt in gopts,bopts,deopts:
    if stefcal_reset_all or (opt.enabled and opt.reset):
      to_remove.append(opt.table);
  if stefcal_reset_all or (stefcal_ifr_gains and stefcal_reset_ifr_gains):
    to_remove.append(stefcal_ifr_gain_table);
  # remove tables
  for fname in to_remove:
    if os.path.exists(fname):
      print("Removing %s as requested"%fname);
      try:
        os.unlink(fname);
      except:
        traceback.print_exc();
        print("Error removing %s"%fname);
    else:
      print("%s does not exist, so not trying to remove"%fname);
  mqs.clearcache('VisDataMux');
  mqs.execute('VisDataMux',mssel.create_io_request(),wait=wait);
