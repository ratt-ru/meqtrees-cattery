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

import Meow
from Meow import ParmGroup,Bookmarks,StdTrees

# This defines some ifr subsets that are commonly used for WSRT data,
# to be offered as defaults in the GUI whereverifrs are selected.
STD_IFR_SUBSETS = [
    "-45 -56 -67",
    "-45 -56 -67 -9A -AB -CD",
    "-45 -46 -56 -67 -68 -9A -AB -CD",
    "FM",
    "FM -9A -9B",
];

TDLCompileOptionSeparator("MS selection");
# MS options first
mssel = Meow.Context.mssel = Meow.MSUtils.MSSelector(has_input=True,tile_sizes=None,
                  read_flags=True,write_flags=True,
                  hanning=True,invert_phases=True);
# MS compile-time options
TDLCompileOptions(*mssel.compile_options());
TDLCompileOption("run_purr","Start Purr on this MS",True);
# MS run-time options
TDLRuntimeMenu("Data selection & flag handling",*mssel.runtime_options());

TDLCompileOptionSeparator("Processing options");

# setup calibration mode menu
# some string constants for the menu entries
CAL = record(VIS="visibilities",AMPL="amplitudes",LOGAMPL="log-amplitudes",PHASE="phases",
  DATA = "data vs. M.E.",
  DIFF = "(data - uv model) vs. M.E."
);
cal_type_opt = TDLOption('cal_type',"Equation type",[CAL.DATA,CAL.DIFF]);
cal_what_opt = TDLOption('cal_what',"Calibrate on",[CAL.VIS,CAL.AMPL,CAL.LOGAMPL,CAL.PHASE]);
 
cal_options = [ cal_type_opt,cal_what_opt ];

# if table access is available, add baseline selection options
if Meow.MSUtils.TABLE:
  calib_ifrs_opt =  TDLOption('calibrate_ifrs',"...using interferometers",
    ["all"]+STD_IFR_SUBSETS,more=str,doc=Meow.IfrArray.ifr_spec_syntax);
  cal_options.append(calib_ifrs_opt);
  mssel.when_changed(lambda msname:calib_ifrs_opt.set_option_list(mssel.ms_ifr_subsets));

read_ms_model_opt = TDLCompileOption("read_ms_model","Read uv-model visibilities from MS",False);
read_ms_model_opt.when_changed(cal_type_opt.show);

cal_toggle = TDLCompileMenu("Calibrate (fit corrupted model to data)",
     toggle='do_solve',open=True,doc="""Select this to include calibration in your tree.""",
     *cal_options
  );

CORRECTED_DATA = "corrected data";
RESIDUALS = "uncorrected residuals";
CORRECTED_RESIDUALS = "corrected residuals";
MODEL = "sky model";
CORRUPTED_MODEL = "corrupted model";

output_option = TDLCompileOption('do_output',"Output visibilities",[CORRECTED_DATA,RESIDUALS,CORRECTED_RESIDUALS,CORRUPTED_MODEL]);
flag_options = TDLCompileMenu("Flag out-of-bounds Jones terms",toggle='do_correct_flag',doc="""If selected,
      your tree will flag visiblity points where the norm of a Jones term is above or below a threshold.""",
    *(
            TDLOption('correct_flag_jmax',"Flag if |J|>",[None,10,100],more=float),
            TDLOption('correct_flag_jmin',"Flag if |J|<",[None,.1,.01],more=float)
    ));

do_correct_sky = False;
# correct_sky_options = TDLCompileOption('do_correct_sky',"...include sky-Jones correction for first source in model",False));
# output_option.when_changed(lambda output:correct_sky_options.show(output in [CORRECTED_RESIDUALS,CORRECTED_DATA]));
output_option.when_changed(lambda output:flag_options.show(output in [CORRECTED_RESIDUALS,CORRECTED_DATA]));

# now load optional modules for the ME maker
from Meow import MeqMaker
meqmaker = MeqMaker.MeqMaker(solvable=True,
                            use_jones_inspectors=True,
                            use_skyjones_visualizers=False,
                            use_decomposition=False
);

flag_options.when_changed(mssel.enable_write_flags);

# specify available sky models
# these will show up in the menu automatically
from Calico.OMS import central_point_source
from Siamese.OMS import fitsimage_sky
import Meow.LSM
lsm = Meow.LSM.MeowLSM(include_options=False);

meqmaker.add_sky_models([lsm,central_point_source,fitsimage_sky]);

# now add optional Jones terms
# these will show up in the menu automatically

# E - beam
# add a fixed primary beam first
from Calico.OMS import wsrt_beams
from Calico.OMS import solvable_pointing_errors
meqmaker.add_sky_jones('E','primary beam',[wsrt_beams],
  pointing=solvable_pointing_errors);
# then add differential gains
from Calico.OMS import solvable_sky_jones
meqmaker.add_sky_jones('dE','differential gains',
  [ solvable_sky_jones.DiagRealImag('dE'),
    solvable_sky_jones.FullRealImag('dE') ]);

# P - feed angle
from Siamese.OMS import feed_angle
meqmaker.add_uv_jones('P','feed orientation',[feed_angle]);

# B - bandpass, G - gain
from Calico.OMS import solvable_jones
meqmaker.add_uv_jones('B','bandpass',
  [ solvable_jones.DiagRealImag("B"),
    solvable_jones.FullRealImag("B"),
    solvable_jones.DiagAmplPhase("B") ]);
meqmaker.add_uv_jones('G','receiver gains/phases',
  [ solvable_jones.DiagRealImag("G"),
    solvable_jones.FullRealImag("G"),
    solvable_jones.DiagAmplPhase("G") ]);

from Calico.OMS import ifr_based_errors
meqmaker.add_vis_proc_module('IG','multiplicative IFR errors',[ifr_based_errors.IfrGains()]);
meqmaker.add_vis_proc_module('IC','additive IFR errors',[ifr_based_errors.IfrBiases()]);

# very important -- insert meqmaker's options properly
TDLCompileOptions(*meqmaker.compile_options());

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

  # all inspector nodes will be added to this list
  inspectors = [];
  # make spigot nodes for data
  if do_solve or do_output not in [MODEL,CORRUPTED_MODEL]:
    mssel.enable_input_column(True);
    spigots = spigots0 = outputs = array.spigots(corr=mssel.get_corr_index());

    # ...and an inspector for them
    StdTrees.vis_inspector(ns.inspector('data'),spigots,
                                bookmark="Inspect data visibilities");
    inspectors += [ ns.inspector('data') ];
    Bookmarks.make_node_folder("Data visibilities by baseline",
      [ spigots(p,q) for p,q in array.ifrs() ],sorted=True,ncol=2,nrow=2);
  else:
    mssel.enable_input_column(False);
    spigots = spigots0 = None;

  # make spigot nodes for model
  corrupt_uvdata = model_spigots = None;
  if read_ms_model:
    mssel.enable_model_column(True);
    model_spigots = array.spigots(column="PREDICT",corr=mssel.get_corr_index());
    # ...and an inspector for them
    StdTrees.vis_inspector(ns.inspector('model'),model_spigots,
                                bookmark="Inspect model visibilities");
    inspectors += [ ns.inspector('model') ];
    Bookmarks.make_node_folder("Model visibilities by baseline",
      [ model_spigots(p,q) for p,q in array.ifrs() ],sorted=True,ncol=2,nrow=2);
    # if calibrating on (input-corrupt model), make corrupt model
    if do_solve and cal_type == CAL.DIFF:
      corrupt_uvdata = meqmaker.corrupt_uv_data(ns,model_spigots);

  inspect_ifrs = array.ifrs();
  # if needed, then make a predict tree using the MeqMaker
  if do_solve or do_output != CORRECTED_DATA:
    if model_spigots and not corrupt_uvdata:
      uvdata = model_spigots;
    else:
      uvdata = None;
    predict = meqmaker.make_predict_tree(ns,uvdata=uvdata);
    # make a ParmGroup and solve jobs for source parameters, if we have any
    if do_solve:
      parms = {};
      for src in meqmaker.get_source_list(ns):
        parms.update([(p.name,p) for p in src.get_solvables()]);
      if parms:
        pg_src = ParmGroup.ParmGroup("source",parms.values(),
                    table_name="sources.fmep",
                    individual=True,bookmark=True);
        # now make a solvejobs for the source
        ParmGroup.SolveJob("cal_source","Calibrate source model",pg_src);
  else:
    predict = None;

  # make nodes to compute residuals
  if do_output in [CORRECTED_RESIDUALS,RESIDUALS]:
    residuals = ns.residuals;
    for p,q in array.ifrs():
      if corrupt_uvdata:
        residuals(p,q) << Meq.Subtract(spigots(p,q),corrupt_uvdata(p,q),predict(p,q));
      else:
        residuals(p,q) << spigots(p,q) - predict(p,q);
    outputs = residuals;

  # and now we may need to correct the outputs
  if do_output in [CORRECTED_DATA,CORRECTED_RESIDUALS]:
    if do_correct_sky:
      srcs = meqmaker.get_source_list(ns);
      sky_correct = srcs and srcs[0];
    else:
      sky_correct = None;
    global do_correct_flag;
    if do_correct_flag and correct_flag_jmin is None and correct_flag_jmax is None:
      do_correct_flag = False;
    outputs = meqmaker.correct_uv_data(ns,outputs,sky_correct=sky_correct,
                                      flag_jones_minmax=do_correct_flag and (correct_flag_jmin,correct_flag_jmax),
                                      inspect_ifrs=inspect_ifrs);
  elif do_output == CORRUPTED_MODEL:
    outputs = predict;

  # make solve trees
  if do_solve:
    # parse ifr specification
    solve_ifrs  = array.subset(calibrate_ifrs,strict=False).ifrs();
    if not solve_ifrs:
      raise RuntimeError,"No interferometers selected for calibration. Check your ifr specification under calibration options.";
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

  # make sinks and vdm.
  # The list of inspectors must be supplied here
  inspectors += meqmaker.get_inspectors() or [];
  StdTrees.make_sinks(ns,outputs,spigots=spigots0,post=inspectors);
  Bookmarks.make_node_folder("Corrected/residual visibilities by baseline",
    [ outputs(p,q) for p,q in array.ifrs() ],sorted=True,ncol=2,nrow=2);

  if not do_solve:
    name = "Generate "+do_output;
    comment = "Generated "+do_output;
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
