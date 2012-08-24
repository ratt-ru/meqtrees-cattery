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
import fnmatch

import Meow
from Meow import ParmGroup,Bookmarks,StdTrees

# This defines some ifr subsets that are commonly used for WSRT data,
# to be offered as defaults in the GUI wherever ifrs are selected.
# MSUtils already has this for WSRT, so reuse it here.
STD_IFR_SUBSETS = Meow.MSUtils.STD_IFR_SUBSETS["WSRT"];

TDLCompileOptionSeparator("MS selection");
# MS options first
mssel = Meow.Context.mssel = Meow.MSUtils.MSSelector(has_input=True,tile_sizes=None,
                  read_flags=True,write_flags=True,
                  hanning=True,invert_phases=True);
# MS compile-time options
TDLCompileOptions(*mssel.compile_options());
TDLCompileOption("run_purr","Start Purr on this MS",False);
# MS run-time options
TDLRuntimeMenu("Data selection & flag handling",*mssel.runtime_options());

TDLCompileOptionSeparator("Processing options");

# setup calibration mode menu
# some string constants for the menu entries
CAL = record(VIS="VIS",AMPL="AMPL",LOGAMPL="LOGAMPL",PHASE="PHASE",DATA="DATA",DIFF="DIFF");

cal_type_opt = TDLOption('cal_type',"Equation type",
                         [  (CAL.DATA,"data vs. predict"),
                            (CAL.DIFF,"(data - uvmodel) vs. predict") ],
  doc="""<P>This determines how the calibration equations are structured. The first setting is
  for normal calibration (data is compared to a predict:)</P>

  <P align="center">corrupt(sky_model + uv_model) &rarr; data</P>

  <P>The second setting is used if your uv-model column is already corrupted:</P>

  <P align="center">corrupt(sky_model) &rarr; data - uv_model</P>
  """
);
cal_what_opt = TDLOption('cal_what',"Calibrate on",
                        [(CAL.VIS,"complex visibilities"),(CAL.AMPL,"amplitudes"),(CAL.LOGAMPL,"log(amplitudes)"),(CAL.PHASE,"phases")],doc="""
  <P>Select "complex visibilities" to directly fit a complex model to complex data, using the real and imaginary
parts. The other options (warning: highly experimental!) will fit only the complex amplitudes (or logs of amplitudes) or phases.</P>
""");

cal_options = [ cal_type_opt,cal_what_opt ];

# if table access is available, add baseline selection options
if Meow.MSUtils.TABLE:
  calib_ifrs_opt =  TDLOption('calibrate_ifrs',"...using interferometers",
    ["all"]+STD_IFR_SUBSETS,more=str,doc="""<P>
      You can restrict calibration to a subset of interferometers. Note that this selection
      applies on top of (not instead of) the global interferometer selection specified above.
      """+Meow.IfrArray.ifr_spec_syntax);
  cal_options.append(calib_ifrs_opt);
  mssel.when_changed(lambda msname:calib_ifrs_opt.set_option_list(mssel.ms_ifr_subsets));

read_ms_model_opt = TDLCompileOption("read_ms_model","Read additional uv-model visibilities from MS",False,doc="""
  <P>If enabled, then an extra set of <i>model</i> visibilities will be read from a second column
  of the MS (in addition to the input data.) These can either be added to whatever is predicted by the
  sky model <i>in the uv-plane</i> (i.e. subject to uv-Jones but not sky-Jones corruptions), or directly
  subtracted from the input data. See also the "Equation type" option below.

  <P>If you are repeatedly running a large sky model, and are not solving for any image-plane effects, then this
  feature lets you compute your sky model (or a part of it) just once, save the result to the uv-model column,
  and reuse it in subsequent steps. This can save a lot of processing time.</P>
  """);
read_ms_model_opt.when_changed(cal_type_opt.show);

cal_toggle = TDLCompileMenu("Calibrate (fit corrupted model to data)",
     toggle='do_solve',open=True,doc="""<P>Select this to include a calibration step in your tree. Calibration
  involves comparing predicted visibilities to input data, and iteratively adjusting the sky and/or
  instrumental model for the best fit.
  </P>
  """,*cal_options
  );

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

FIRST_SOURCE = "first source in model";
correct_sky_options = TDLCompileOption('do_correct_sky',"Include sky-Jones correction for which source",
  [None,FIRST_SOURCE],more=str,doc=
  """<P>If enabled, then any sky-Jones terms associated with the named source will be corrected
  for (note that individual Jones terms may be excluded from the correction via their
  Advanced options menu). The source name may include shell-style wildcards such as "*" and "?",
  in which case the first matching source will be used.</P>""");
output_option.when_changed(lambda output:correct_sky_options.show(output in [CORRECTED_RESIDUALS,CORRECTED_DATA]));

flag_jones_opt = TDLOption("flag_jones","Flag on out-of-bounds Jones terms",False,
    doc="""<P>Your tree can automatically flag visibility points where certain Jones terms go out of bounds.
    The relevant Jones terms will have individual flagging options in their "Advanced options" menus, which
    you must also enable for this option to have effect. This option enables and disables Jones-bassed
    flagging globally.</P>""");
flag_res_opt = TDLOption("flag_res","Flag on residual amplitudes >",[None],more=float,
    doc="""<P>If selected, your tree will flag all visibilities (per IFR/timeslot/channel) where the residual
    complex amplitude exceeds the given value.</P>
    """);
flag_meanres_opt = TDLOption("flag_mean_res","Flag on mean residual amplitudes (over all IFRs) >",[None],more=float,
    doc="""<P>If selected, your tree will flag entire timeslots (per channel) where the mean residual
    complex amplitude over all IFRs within the timeslots exceeds the given value.</P>
    """);

flag_menu = TDLCompileMenu("Flag output visibilities",flag_jones_opt,flag_res_opt,flag_meanres_opt,toggle="flag_enable");

def _select_output (output):
  if output in [RESIDUALS,CORRECTED_RESIDUALS,CORRECTED_DATA]:
    flag_menu.show(True);
    flag_jones_opt.show(output in [CORRECTED_RESIDUALS,CORRECTED_DATA]);
#    flag_res_opt.show(output in [RESIDUALS,CORRECTED_RESIDUALS]);
#    flag_meanres_opt.show(output in [RESIDUALS,CORRECTED_RESIDUALS]);
  else:
    flag_menu.show(False);

output_option.when_changed(_select_output);

# now load optional modules for the ME maker
from Meow import TensorMeqMaker
meqmaker = TensorMeqMaker.TensorMeqMaker(solvable=True,
                            use_correction=True,
                            use_jones_inspectors=True,
                            use_skyjones_visualizers=False,
);

flag_menu.when_changed(mssel.enable_write_flags);

# specify available sky models
# these will show up in the menu automatically
from Calico.OMS import central_point_source
from Siamese.OMS import fitsimage_sky,gridded_sky
models = [central_point_source,fitsimage_sky,gridded_sky]

try:
  from Siamese.OMS.tigger_lsm import TiggerSkyModel
  models.insert(0,TiggerSkyModel());
except:
  traceback.print_exc();
  pass;

meqmaker.add_sky_models(models);

# E - beam
# add a fixed primary beam first
from Calico.OMS.wsrt_cos3_beam import WSRTCos3Beam #,wsrt_beams_zernike
meqmaker.add_sky_jones('E','WSRT cos^3 beam',[WSRTCos3Beam("E",solvable=False)]);
meqmaker.add_sky_jones('Es','solvable WSRT cos^3 beam',[WSRTCos3Beam("Es",solvable=True)]);
## add solvable refraction
# from Calico.OMS import solvable_position_shifts
# meqmaker.add_sky_jones('R','position shifts',solvable_position_shifts);
# add differential gains
from Calico.OMS import solvable_sky_jones
meqmaker.add_sky_jones('dE','differential gains',
  [ solvable_sky_jones.DiagRealImag('dE'),
    solvable_sky_jones.FullRealImag('dE'),
    solvable_sky_jones.DiagAmplPhase('dE')
  ]);

# P - feed angle
from Siamese.OMS import feed_angle
meqmaker.add_uv_jones('P','feed orientation',[feed_angle]);

# B - bandpass, G - gain
from Calico.OMS import solvable_jones
meqmaker.add_uv_jones('B','bandpass',
  [ solvable_jones.DiagRealImag("B"),
    solvable_jones.FullRealImag("B"),
    solvable_jones.DiagAmplPhase("B") ],flaggable=True);
meqmaker.add_uv_jones('G','receiver gains/phases',
  [ solvable_jones.DiagRealImag("G"),
    solvable_jones.FullRealImag("G"),
    solvable_jones.DiagAmplPhase("G") ],flaggable=True);

from Calico.OMS import ifr_based_errors
meqmaker.add_vis_proc_module('IG','multiplicative interferometer-based gains (IG)',[ifr_based_errors.IfrGains()]);
meqmaker.add_vis_proc_module('IC','additive interferometer-based biases (IC)',[ifr_based_errors.IfrBiases()]);

TDLCompileOption("enable_inspectors","Enable extra visualization nodes (inspectors, etc.)",True);

# very important -- insert meqmaker's options properly
TDLCompileOptions(*meqmaker.compile_options());

# TDLCompileOption("include_inspectors","Include inspectors",True);

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
  mssel.setup_observation_context(ns,prefer_baseline_uvw=True);
  array = Meow.Context.array;

  # make spigot nodes for data
  if do_solve or do_output not in [CORRUPTED_MODEL]:
    mssel.enable_input_column(True);
    spigots = spigots0 = outputs = array.spigots(corr=mssel.get_corr_index());
    if enable_inspectors:
      meqmaker.make_per_ifr_bookmarks(spigots,"Input visibilities");
    # add IFR-based errors, if any
    spigots = meqmaker.apply_visibility_processing(ns,spigots);
  else:
    mssel.enable_input_column(False);
    spigots = spigots0 = None;

  # make spigot nodes for model
  corrupt_uvdata = model_spigots = None;
  if read_ms_model:
    mssel.enable_model_column(True);
    model_spigots = array.spigots(column="PREDICT",corr=mssel.get_corr_index());
    if enable_inspectors:
      meqmaker.make_per_ifr_bookmarks(model_spigots,"UV-model visibilities");
    # if calibrating on (input-corrupt model), make corrupt model
    if do_solve and cal_type == CAL.DIFF:
      corrupt_uvdata = meqmaker.corrupt_uv_data(ns,model_spigots);

  # if needed, then make a predict tree using the MeqMaker
  if do_solve or do_output != CORRECTED_DATA:
    if model_spigots and not corrupt_uvdata:
      uvdata = model_spigots;
    else:
      uvdata = None;
    predict = meqmaker.make_predict_tree(ns,uvdata=uvdata);
  else:
    predict = None;
  output_title = "Uncorrected residuals";

  # make nodes to compute residuals
  if do_output in [CORRECTED_RESIDUALS,RESIDUALS]:
    residuals = ns.residuals;
    for p,q in array.ifrs():
      if corrupt_uvdata:
        residuals(p,q) << Meq.Subtract(spigots(p,q),corrupt_uvdata(p,q),predict(p,q));
      else:
        residuals(p,q) << spigots(p,q) - predict(p,q);
    if enable_inspectors:
      meqmaker.make_per_ifr_bookmarks(residuals,"Uncorrected residuals");
    outputs = residuals;

  # and now we may need to correct the outputs
  if do_output in [CORRECTED_DATA,CORRECTED_RESIDUALS]:
    if do_correct_sky:
      srcs = meqmaker.get_source_list(ns);
      if do_correct_sky is FIRST_SOURCE:
        sky_correct = srcs and srcs[0];
      else:
        srcs = [ src for src in srcs if fnmatch.fnmatchcase(src.name,do_correct_sky) ];
        sky_correct = srcs and srcs[0];
    else:
      sky_correct = None;
    outputs = meqmaker.correct_uv_data(ns,outputs,sky_correct=sky_correct,
                                      flag_jones=flag_jones);
    output_title = "Corrected data" if do_output is CORRECTED_DATA else "Corrected residuals";
  elif do_output == CORRUPTED_MODEL:
    outputs = predict;
    output_title = "Predict";
  elif do_output == CORRUPTED_MODEL_ADD:
    outputs = ns.output;
    for p,q in array.ifrs():
      outputs(p,q) << spigots(p,q) + predict(p,q);
    output_title = "Data+predict";

  # make flaggers
  if flag_enable and do_output in [ CORRECTED_DATA,RESIDUALS,CORRECTED_RESIDUALS ]:
    flaggers = [];
    if flag_res is not None or flag_mean_res is not None:
      for p,q in array.ifrs():
        ns.absres(p,q) << Meq.Abs(outputs(p,q));
    # make flagger for residuals
    if flag_res is not None:
      for p,q in array.ifrs():
        ns.flagres(p,q) << Meq.ZeroFlagger(ns.absres(p,q)-flag_res,oper='gt',flag_bit=Meow.MSUtils.FLAGMASK_OUTPUT);
      flaggers.append(ns.flagres);
      # ...and an inspector for them
      if enable_inspectors:
        meqmaker.make_per_ifr_bookmarks(ns.flagres,"Residual amplitude flags");
    # make flagger for mean residuals
    if flag_mean_res is not None:
      ns.meanabsres << Meq.Mean(*[ns.absres(p,q) for p,q in array.ifrs()]);
      ns.flagmeanres << Meq.ZeroFlagger(ns.meanabsres-flag_mean_res,oper='gt',flag_bit=Meow.MSUtils.FLAGMASK_OUTPUT);
      Meow.Bookmarks.Page("Mean residual amplitude flags").add(ns.flagmeanres,viewer="Result Plotter");
      flaggers.append(lambda p,q:ns.flagmeanres);

    # merge flags into output
    if flaggers:
      if enable_inspectors:
        meqmaker.make_per_ifr_bookmarks(outputs,output_title+" (unflagged)");
      for p,q in array.ifrs():
        ns.flagged(p,q) << Meq.MergeFlags(outputs(p,q),*[f(p,q) for f in flaggers]);
      outputs = ns.flagged;

  if enable_inspectors:
    meqmaker.make_per_ifr_bookmarks(outputs,output_title);
    abs_outputs = outputs('abs');
    for p,q in array.ifrs():
      abs_outputs(p,q) << Meq.Abs(outputs(p,q));
    meqmaker.make_per_ifr_bookmarks(abs_outputs,output_title+" (mean amplitudes)");

  # make solve trees
  if do_solve:
    # parse ifr specification
    solve_ifrs  = array.subset(calibrate_ifrs,strict=False).ifrs();
    if not solve_ifrs:
      raise RuntimeError,"No interferometers selected for calibration. Check your ifr specification (under calibration options).";
    # inputs to the solver are based on calibration type
    if corrupt_uvdata:
      [ ns.diff(p,q) << spigots(p,q) - corrupt_uvdata(p,q) for p,q in solve_ifrs ];
      rhs = ns.diff;
    else:
      rhs = spigots;
    lhs = predict;
    weights = modulo = None;
    # if calibrating visibilities, feed them to condeq directly, else take ampl/phase
    if cal_what == CAL.VIS:
      pass;
    elif cal_what == CAL.AMPL:
      [ x('ampl',p,q) << Meq.Abs(x(p,q)) for p,q in ifrs for x in rhs,lhs ];
      lhs = lhs('ampl');
      rhs = rhs('ampl');
    elif cal_what == CAL.LOGAMPL:
      [ x('logampl',p,q) << Meq.Log(Meq.Abs(x(p,q))) for p,q in ifrs for x in rhs,lhs ];
      lhs = lhs('logampl');
      rhs = rhs('logampl');
    elif cal_what == CAL.PHASE:
      [ x('phase',p,q) << Meq.Arg(x(p,q)) for p,q in ifrs for x in rhs,lhs ];
      [ rhs('ampl',p,q) << Meq.Abs(rhs(p,q)) for p,q in ifrs  ];
      lhs = lhs('phase');
      rhs = rhs('phase');
      weights = rhs('ampl');
      modulo = 2*math.pi;
    else:
      raise ValueError,"unknown cal_what setting: "+str(cal_what);
    # make a solve tree
    solve_tree = StdTrees.SolveTree(ns,lhs,solve_ifrs=solve_ifrs,weights=weights,modulo=modulo);
    # the output of the sequencer is either the residuals or the spigots,
    # according to what has been set above
    outputs = solve_tree.sequencers(inputs=rhs,outputs=outputs);

  post = ( ( enable_inspectors and meqmaker.get_inspectors() ) or [] );
  StdTrees.make_sinks(ns,outputs,spigots=spigots0,post=post,corr_index=mssel.get_corr_index());

  if not do_solve:
    name = "Generate "+output_title.lower();
    comment = "Generated "+output_title.lower();
    if name:
      # make a TDL job to run the tree
      def run_tree (mqs,parent,wait=False,**kw):
        global tile_size;
        purrpipe.title("Calibrating").comment(comment);
        return mqs.execute(Meow.Context.vdm.name,mssel.create_io_request(tile_size),wait=wait);
      TDLRuntimeMenu(name,
        TDLOption('tile_size',"Tile size, in timeslots",[10,60,120,240],more=int,
                  doc="""Input data is sliced by time, and processed in chunks (tiles) of
                  the indicated size. Larger tiles are faster, but use more memory."""),
        TDLJob(run_tree,name,job_id='generate_visibilities')
      );

  # very important -- insert meqmaker's runtime options properly
  # this should come last, since runtime options may be built up during compilation.
  TDLRuntimeOptions(*meqmaker.runtime_options(nest=False));
  # insert solvejobs
  if do_solve:
    TDLRuntimeOptions(*ParmGroup.get_solvejob_options());
  # finally, setup imaging options
  imsel = mssel.imaging_selector(npix=512,arcmin=meqmaker.estimate_image_size());
  TDLRuntimeMenu("Make an image from this MS",*imsel.option_list());

  # and close meqmaker -- this exports annotations, etc
  meqmaker.close();
