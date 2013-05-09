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

# specify available sky models
# these will show up in the menu automatically
#from Calico.OMS import central_point_source
#from Siamese.OMS import fitsimage_sky,gridded_sky
#models = [central_point_source,fitsimage_sky,gridded_sky]

models = [];
try:
  from Siamese.OMS.tigger_lsm import TiggerSkyModel
  models.insert(0,TiggerSkyModel());
except:
  traceback.print_exc();
  pass;

meqmaker.add_sky_models(models);

# E - beam
# add a fixed primary beam first
from Calico.OMS import wsrt_beams  #,wsrt_beams_zernike

meqmaker.add_sky_jones('E','primary beam',[wsrt_beams]); # ,wsrt_beams_zernike]);

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
RESIDUALS = "RES";
CORRECTED_RESIDUALS = "CORR_RES";
CORRUPTED_MODEL = "PREDICT";
CORRUPTED_MODEL_ADD = "DATA+PREDICT";

output_option = TDLCompileOption('do_output',"Output visibilities",
                                 [  (CORRECTED_DATA,"corrected data"),
                                    (RESIDUALS,"uncorrected residuals"),
                                    (CORRECTED_RESIDUALS,"corrected residuals"),
                                    (CORRUPTED_MODEL,"predict"),
                                    (CORRUPTED_MODEL_ADD,"data+predict")],
  doc="""<P>This selects what sort of visibilities get written to the output column:</P>
  <ul>

  <li><B>Predict</B> refers to the visibilities given by the sky model (plus an optional uv-model column),
  corrupted by the current instrumental model using the Measurement Equation specified below.</li>

  <li><B>Corrected data</B> is the input data corrected for the instrumental model (by applying the inverse of the
  M.E.)</li>

  <li><B>Uncorrected residuals</B> refer to input data minus predict. This corresponds to whatever signal is
  left in your data that is <b>not</b> represented by the model, and still subject to instrumental corruptions.</li>

  <li><B>Corrected residuals</B> are residuals corrected for the instrumental model. This is what you usually
  want to see during calibration.</li>

  <li><B>Data+predict</B> is a special mode where the predict is <i>added</I> to the input data. This is used
  for injecting synthetic sources into your data, or for accumulating a uv-model in several steps. (In
  the latter case your input column needs to be set to the uv-model column.)</li>
  </ul>

  </P>If calibration is enabled above, then a calibration step is executed prior to generating output data. This
  will update the instrumental and/or sky models. If calibration is not enabled, then the current models
  may still be determined by the results of past calibration, since these are stored in persistent <i>MEP
  tables.</i></P>

  """);


diffgain_tag = 'dE';
diffgain_group = 'cluster';

MODE_SOLVE_SAVE = "solve-save";
MODE_SOLVE_NOSAVE = "solve-nosave"
MODE_SOLVE_APPLY = "apply"

TDLCompileOption("stefcal_implementation","Stefcal implementation",
  ["GainDiag","Gain2x2","Gain2x2-r8768.Gain2x2" ] );
TDLCompileOption("stefcal_gain_mode","Solution mode",
  {MODE_SOLVE_SAVE:"solve and save",MODE_SOLVE_NOSAVE:"solve, do not save",MODE_SOLVE_APPLY:"load and apply"}); 
TDLCompileOption("stefcal_gain_reset","Ignore previously saved G/dE solutions",False);
TDLCompileMenu("Filenames for solution tables",
  TDLOption("stefcal_gain_table","Gains",["gains.cp"],more=str),
  TDLOption("stefcal_diff_gain_table","Diferential gains",["diffgains.cp"],more=str),
);
dgsel = TensorMeqMaker.SourceSubsetSelector("Enable source-based differential gains (dEs)",
          tdloption_namespace='de_subset',annotate=False);
TDLCompileOptions(*dgsel.options);
TDLCompileMenu("Convergence settings (G)",
  TDLOption("stefcal_weigh_g","Apply noise-based weights",True),
  TDLOption("stefcal_niter_g","Max iterations, first major loop",[20,50,100],more=int,default=50),
  TDLOption("stefcal_niter1_g","Max iterations, middle loops",[20,50,100],more=int,default=20),
  TDLOption("stefcal_niter2_g","Max iterations, last loop",[20,50,100],more=int,default=50),
  TDLOption("stefcal_epsilon_g","Convergence epsilon",[1e-4,1e-5,1e-6],more=float,default=1e-6),
  TDLOption("stefcal_quota_g","Convergence epsilon quota",[.95,.99,1],more=float),
  TDLOption("stefcal_delta","Convergence delta, first major loop",[1e-5,1e-6,1e-7],more=float,default=1e-6),
  TDLOption("stefcal_delta_1","Convergence delta, middle loops",[1e-5,1e-6,1e-7],more=float,default=1e-5),
  TDLOption("stefcal_delta_2","Convergence delta, last loop",[1e-5,1e-6,1e-7],more=float,default=1e-5),
  TDLOption("stefcal_omega","Averaging weight (omega), first major loop",[0.5,1.8],more=float),
  TDLOption("stefcal_average","Averaging mode, first major loop",[0,1,2],default=2),
  TDLOption("stefcal_ff","Enable feed-forward averaging, first major loop",False),
  TDLOption("stefcal_omega_1","Averaging weight (omega), rest of major loop",[0.5,1.8],more=float),
  TDLOption("stefcal_average_1","Averaging mode, rest of major loop",[0,1,2],default=2),
  TDLOption("stefcal_ff_1","Enable feed-forward averaging, rest of major loop",False),
  TDLOption("stefcal_g_timeint","Solution interval, time axis (0 for full axis)",[0,1],more=int,default=1),
  TDLOption("stefcal_g_freqint","Solution interval, freq axis (0 for full axis)",[0,1],more=int,default=1),
  TDLOption("stefcal_g_timesmooth","Smoothing kernel, time axis",[0],more=int),
  TDLOption("stefcal_g_freqsmooth","Smoothing kernel, freq axis",[0],more=int),
);
TDLCompileMenu("Convergence settings (dE)",
  TDLOption("stefcal_weigh_de","Apply noise-based weights",True),
  TDLOption("stefcal_niter_de","Max iterations",[20,50,100],more=int),
  TDLOption("stefcal_epsilon_de","Convergence epsilon",[1e-4,1e-5,1e-6],more=float),
  TDLOption("stefcal_quota_de","Convergence epsilon quota",[.95,.99,1],more=float,default=.99),
  TDLOption("stefcal_delta_de","Convergence delta",[1e-5,1e-6,1e-7],more=float,default=1e-6),
  TDLOption("stefcal_omega_de","Averaging weight (omega)",[0.5,1.8],more=float),
  TDLOption("stefcal_average_de","Averaging mode",[0,1,2],default=2),
  TDLOption("stefcal_ff_de","Enable feed-forward averaging",False),
  TDLOption("stefcal_de_timeint","Solution interval, time axis (0 for full axis)",[0,30,60,120],more=int,default=60),
  TDLOption("stefcal_de_freqint","Solution interval, freq axis (0 for full axis)",[0,30,60,120],more=int,default=0),
  TDLOption("stefcal_de_timesmooth","Smoothing kernel, time axis",[0],more=int),
  TDLOption("stefcal_de_freqsmooth","Smoothing kernel, freq axis",[0],more=int),
);
TDLCompileOption("stefcal_nmajor","Number of major loops",[1,2,3,5],more=int,default=2);
TDLCompileOption("stefcal_rescale","Rescale data to model before solving",True);
DIAGONLY,ALLFOUR = "diag","full";
TDLCompileMenu("Enable IFR-based gains",
  TDLOption("stefcal_ifr_gain_mode","Solution mode",
    {MODE_SOLVE_SAVE:"solve and save",MODE_SOLVE_NOSAVE:"solve, do not save",MODE_SOLVE_APPLY:"load and apply"}),
  TDLOption("stefcal_ifr_gain_reset","Ignore previously saved solutions",False),
  TDLOption("stefcal_diagonal_ifr_gains","Polarizations",
            {DIAGONLY:"parallel-hand only",ALLFOUR:"full 2x2"}),
  TDLOption("stefcal_per_chan_ifr_gains","Solve on a per-channel basis",False),
  TDLOption("stefcal_ifr_gain_table","Filename for solutions",["ifrgains.cp"],more=str),
  toggle="stefcal_ifr_gains",
);

EQTYPE_MODEL = "model";
EQTYPE_DATA  = "data";
TDLCompileOption("stefcal_eqtype","Equation type",{EQTYPE_MODEL:"GMG*->D",EQTYPE_DATA:"GDG*->M"});
TDLCompileMenu("Include visualizers...",
  TDLCompileOption("visualize_G","For G solutions",False),
  TDLCompileOption("visualize_dE","For dE solutions",False),
  TDLCompileOption("visualize_flag_unity","Flag zero/unity solutions in visualizers",True),
  TDLCompileOption("visualize_norm_offdiag","Normalize off-diagonal terms by diagonals",True),
  toggle="stefcal_visualize");
TDLCompileOption("stefcal_verbose","Verbosity level",[0,1,2,3],more=int);

import Purr.Pipe

def _define_forest(ns,parent=None,**kw):
  if not mssel.msname:
    raise RuntimeError,"MS not set";
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

  # predict tree using the MeqMaker
  all_sources = meqmaker.get_source_list(ns);
  dg_sources = dgsel.subset_enabled and dgsel.filter(all_sources);
  if dg_sources:
    # first group all sources without a diffgain on them
    groups = [ [ src for src in all_sources if not src in dg_sources ] ];
    # group diffgain-enabled sources by grouping tag
    clusters = set([src.get_attr(diffgain_group,None) for src in dg_sources]);
    dg_groups = [ [ src for src in dg_sources if src.get_attr(diffgain_group) == name ] for name in clusters if name ];
    # add sources without a grouping tag individually, as single-source groups
    dg_groups += [ [ src ] for src in dg_sources if not src.get_attr(diffgain_group,None) ];
    # now sort by brightness
    flux_dgg = [ (sum([src.get_attr('Iapp',0) or src.get_attr('I') for src in dgg]),dgg) for dgg in dg_groups ];
    flux_dgg.sort(lambda a,b:cmp(b[0],a[0]));
    groups += [ dgg for flux,dgg in flux_dgg ];
    num_diffgains = len(flux_dgg);
    # now make predict trees
    models = [];
    for i,group in enumerate(groups):
      MT = ns.MT(i);
      predict = meqmaker.make_predict_tree(MT.Subscope(),sources=group);
      ns.MT(i) << Meq.Composer(dims=[0],mt_polling=True,*[ predict(p,q) for p,q in array.ifrs() ]);
      models.append(ns.MT(i));
  else:
    predict = meqmaker.make_predict_tree(ns);
    ns.MT << Meq.Composer(dims=[0],mt_polling=True,*[ predict(p,q) for p,q in array.ifrs() ]);
    models = [ ns.MT ];
    global visualize_dE;
    visualize_dE = False;

  solve_ifrs  = array.subset(calibrate_ifrs,strict=False).ifrs();
  gain_subtiling = [ stefcal_g_timeint,stefcal_g_freqint ];
  gain_smoothing = [ stefcal_g_timesmooth,stefcal_g_freqsmooth ];
  for i in range(2):
    if gain_smoothing[i]:
      gain_subtiling[i] = 1;
  if all([not x for x in gain_smoothing]):
    gain_smoothing = [];

  diffgain_subtiling = [ stefcal_de_timeint,stefcal_de_freqint ];
  diffgain_smoothing = [ stefcal_de_timesmooth,stefcal_de_freqsmooth ];
  for i in range(2):
    if diffgain_smoothing[i]:
      diffgain_subtiling[i] = 1;
  if all([not x for x in diffgain_smoothing]):
    diffgain_smoothing = [];

  import Calico.OMS.StefCal.StefCal
  global visualize_G,visualize_dE;
  visualize_G = stefcal_visualize and visualize_G;
  visualize_dE = stefcal_visualize and visualize_dE;
  ns.stefcal << Meq.PyNode(class_name="StefCalNode",module_name=Calico.OMS.StefCal.StefCal.__file__,
                           ifrs=[ "%s:%s"%(p,q) for p,q in array.ifrs() ],
                           baselines=[ array.baseline(ip,iq) for (ip,p),(iq,q) in array.ifr_index() ],
                           solve_ifrs=[ "%s:%s"%(p,q) for p,q in solve_ifrs ],
                           max_iter=stefcal_niter_g,max_iter_1=stefcal_niter1_g,max_iter_2=stefcal_niter2_g,
                           delta=stefcal_delta,delta_1=stefcal_delta_1,delta_2=stefcal_delta_2,
                           epsilon=stefcal_epsilon_g,
                           weigh_gains=stefcal_weigh_g,
                           diffgain_max_iter=stefcal_niter_de,
                           diffgain_delta=stefcal_delta_de,
                           diffgain_epsilon=stefcal_epsilon_de,
                           weigh_diffgains=stefcal_weigh_de,
                           max_major=stefcal_nmajor,
                           omega=stefcal_omega,omega_1=stefcal_omega_1,omega_de=stefcal_omega_de,
                           average=stefcal_average,average_1=stefcal_average_1,average_de=stefcal_average_de,
                           feed_forward=stefcal_ff,feed_forward_1=stefcal_ff_1,feed_forward_de=stefcal_ff_de,
                           convergence_quota=stefcal_quota_g,diffgain_convergence_quota=stefcal_quota_de,
                           gain_subtiling=gain_subtiling,diffgain_subtiling=diffgain_subtiling,
                           gain_smoothing=gain_smoothing,diffgain_smoothing=diffgain_smoothing,
                           implementation=stefcal_implementation,
                           regularization_factor=1e-6,#
                           rescale=stefcal_rescale,
                           init_from_previous=False,
                           # gain solution options
                           solve_gains=(stefcal_gain_mode != MODE_SOLVE_APPLY),
                           save_gains=(stefcal_gain_mode == MODE_SOLVE_SAVE),
                           reset_gains=stefcal_gain_reset,
                           gain_table=stefcal_gain_table,
                           diff_gain_table=stefcal_diff_gain_table,
                           # IFR gain solution options
                           apply_ifr_gains=stefcal_ifr_gains,
                           solve_ifr_gains=(stefcal_ifr_gain_mode != MODE_SOLVE_APPLY),
                           reset_ifr_gains=stefcal_ifr_gain_reset,
                           save_ifr_gains=(stefcal_ifr_gain_mode == MODE_SOLVE_SAVE),
                           ifr_gain_table=stefcal_ifr_gain_table,
                           per_chan_ifr_gains=stefcal_per_chan_ifr_gains,
                           diag_ifr_gains=(stefcal_diagonal_ifr_gains == DIAGONLY),
#                          gains_on_data=True,
#                           gain_bounds=[.01,1000],
                           gains_on_data=(stefcal_eqtype == EQTYPE_DATA),
#                           gain_bounds=[1e-5,10],
                           visualize_gains=visualize_G,visualize_diffgains=visualize_dE,
#                           dump_domain=0,dump_diffgain=-1,
                           correct=(do_output in (CORRECTED_DATA,CORRECTED_RESIDUALS)),
                           residuals=(do_output in (RESIDUALS,CORRECTED_RESIDUALS)),
                           verbose=stefcal_verbose,
                           children=[ns.DT]+models);

  if visualize_G:
    ns.stefcal_vis_G << Meq.PyNode(class_name="StefCalVisualizer",module_name=Calico.OMS.StefCal.StefCal.__file__,
      label="G",flag_unity=visualize_flag_unity,norm_offdiag=visualize_norm_offdiag,
      vells_label=Context.correlations);
    ns.stefcal_vis_G_avg << Meq.PyNode(class_name="StefCalVisualizer",module_name=Calico.OMS.StefCal.StefCal.__file__,
      label="G",freq_average=True,flag_unity=visualize_flag_unity,norm_offdiag=visualize_norm_offdiag,
      vells_label=Context.correlations);
  if visualize_dE:
    for i in range(num_diffgains):
      ns.stefcal_vis_dE(i) << Meq.PyNode(class_name="StefCalVisualizer",module_name=Calico.OMS.StefCal.StefCal.__file__,
        label="dE:%d"%i,flag_unity=visualize_flag_unity,norm_offdiag=visualize_norm_offdiag,
        vells_label=Context.correlations);
      ns.stefcal_vis_dE_avg(i) << Meq.PyNode(class_name="StefCalVisualizer",module_name=Calico.OMS.StefCal.StefCal.__file__,
                                    label="dE:%d"%i,freq_average=True,flag_unity=visualize_flag_unity,norm_offdiag=visualize_norm_offdiag,
                                    vells_label=Context.correlations);

  nv = 0;
  for p,q in array.ifrs():
    sel = ns.output_sel(p,q) << Meq.Selector(ns.stefcal,index=range(nv,nv+4),multi=True);
    ns.output(p,q) << Meq.Composer(sel,dims=[2,2]);
    nv += 4;
  meqmaker.make_per_ifr_bookmarks(ns.output,"Output visibilities");

  inspectors = meqmaker.get_inspectors() or [];
  if visualize_G:
    Bookmarks.Page("StefCal G plotter").add(ns.stefcal_vis_G,viewer="Result Plotter");
    Bookmarks.Page("StefCal G inspector").add(ns.stefcal_vis_G_avg,viewer="Collections Plotter");
    inspectors += [ ns.stefcal_vis_G,ns.stefcal_vis_G_avg ];
  if visualize_dE:
    for i in range(num_diffgains):
      Bookmarks.Page("StefCal dE:%d plotter"%i).add(ns.stefcal_vis_dE(i),viewer="Result Plotter");
      Bookmarks.Page("StefCal dE:%d inspector"%i).add(ns.stefcal_vis_dE_avg(i),viewer="Collections Plotter");
      inspectors += [ ns.stefcal_vis_dE(i),ns.stefcal_vis_dE_avg(i) ];

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
  for table,varname,desc in [ 
        (stefcal_gain_table,'stefcal_reset_gains',"gain solutions"),
        (stefcal_diff_gain_table,'stefcal_reset_diff_gains',"differential gain solutions"),
        (stefcal_ifr_gain_table,'stefcal_reset_ifr_gains',"IFR-based gain solutions") ]:
    TDLRuntimeOption(varname,"Remove existing %s (%s)"%(desc,os.path.basename(table)),False);
  TDLRuntimeJob(_run_stefcal,"Run StefCal",job_id="stefcal");

def _run_stefcal (mqs,parent,wait=False):
  for table,varname,desc in [ 
        (stefcal_gain_table,'stefcal_reset_gains',"gain solutions"),
        (stefcal_diff_gain_table,'stefcal_reset_diff_gains',"differential gain solutions"),
        (stefcal_ifr_gain_table,'stefcal_reset_ifr_gains',"IFR-based gain solutions") ]:
    if stefcal_reset_all or globals().get(varname,False):
      if os.path.exists(table):
        print "Removing %s as requested"%table;
        try:
          os.unlink(table);
        except:
          traceback.print_exc();
          print "Error removing %s"%table;
      else:
        print "%s does not exist, so not trying to remove"%table;
  mqs.clearcache('VisDataMux');
  mqs.execute('VisDataMux',mssel.create_io_request(),wait=wait);
