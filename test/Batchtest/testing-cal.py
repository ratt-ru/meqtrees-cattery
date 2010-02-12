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

# MS options first
mssel = Meow.Context.mssel = Meow.MSUtils.MSSelector(has_input=True,tile_sizes=None,read_flags=True,
                  hanning=True,invert_phases=True);
# MS compile-time options
TDLCompileOptions(*mssel.compile_options());
TDLCompileOption("run_purr","Start Purr on this MS",True);
# MS run-time options
TDLRuntimeMenu("Data selection & flag handling",*mssel.runtime_options());
## also possible:

# setup output mode menu
CAL = record(VIS="visibilities",AMPL="amplitudes",LOGAMPL="log-amplitudes",PHASE="phases");
cal_options = [
     TDLOption('cal_type',"Calibrate on",[CAL.VIS,CAL.AMPL,CAL.LOGAMPL,CAL.PHASE]),
];
# if table access is available, add baseline selection options
if Meow.MSUtils.TABLE:
  cal_options += [
     TDLOption('min_baseline',"Ignore baselines shorter than or equal to (m)",[None,72,144],more=int),
     TDLOption('max_baseline',"Ignore baselines longer than or equal to  (m)",[None],more=int)   ];

cal_toggle = TDLMenu("Calibrate",
     toggle='do_solve',open=True,
     *cal_options
  );
TDLCompileMenu("What do we want to do",
  cal_toggle,
  TDLOption('do_subtract',"Subtract sky model and generate residuals",True),
  TDLOption('do_correct',"Correct the data or residuals",True),
);
do_correct_sky = False;
  #TDLOption('do_correct_sky',"...include sky-Jones correction for first source in model",True));

# now load optional modules for the ME maker
from Meow import MeqMaker
meqmaker = MeqMaker.MeqMaker(solvable=True,use_jones_inspectors=True);

# disable source decomposition when calibrating
def enable_calibration (enable):
  meqmaker.use_decomposition_opt.show(not enable);
  if enable:
    meqmaker.use_decomposition_opt.set(False);
cal_toggle.when_changed(enable_calibration);

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
meqmaker.add_sky_jones('dE','differential gains',[
    solvable_sky_jones.FullRealImag('dE')]);

# P - feed angle
from Siamese.OMS import feed_angle
meqmaker.add_uv_jones('P','feed orientation',[feed_angle]);

# B - bandpass, G - gain
from Calico.OMS import solvable_jones
meqmaker.add_uv_jones('B','bandpass',
  [ solvable_jones.DiagAmplPhase("B"),
    solvable_jones.FullRealImag("B") ]);
meqmaker.add_uv_jones('G','receiver gains/phases',
  [ solvable_jones.DiagAmplPhase("G"),
    solvable_jones.FullRealImag("G") ]);

from Calico.OMS import ifr_based_errors
meqmaker.add_vis_proc_module('IG','interferometer gains',[ifr_based_errors.IfrGains()]);
meqmaker.add_vis_proc_module('IC','interferometer biases',[ifr_based_errors.IfrBiases()]);

# very important -- insert meqmaker's options properly
TDLCompileOptions(*meqmaker.compile_options());

import Purr.Pipe

def _define_forest(ns,parent=None,**kw):
  if run_purr:
    Timba.TDL.GUI.purr(mssel.msname+".purrlog",[mssel.msname,'.']);
  # create Purr pipe
  global purrpipe;
  purrpipe = Purr.Pipe.Pipe(mssel.msname);
  
  # get antennas from MS
  ANTENNAS = mssel.get_antenna_set(range(1,15));
  array = Meow.IfrArray(ns,ANTENNAS,mirror_uvw=False);
  stas = array.stations();
  # get phase centre from MS, setup observation
  observation = Meow.Observation(ns,phase_centre=mssel.get_phase_dir(),
          linear=mssel.is_linear_pol(),
          circular=mssel.is_circular_pol());
  Meow.Context.set(array,observation);
  # get active correlations from MS
  Meow.Context.active_correlations = mssel.get_correlations();
  
  # make spigot nodes
  spigots = spigots0 = outputs = array.spigots(corr=mssel.get_corr_index());

  # ...and an inspector for them
  StdTrees.vis_inspector(ns.inspector('input'),spigots,
                              bookmark="Inspect input visibilities");
  inspectors = [ ns.inspector('input') ];
  Bookmarks.make_node_folder("Input visibilities by baseline",
    [ spigots(p,q) for p,q in array.ifrs() ],sorted=True,ncol=2,nrow=2);

  inspect_ifrs = array.ifrs();
  if do_solve:
    # filter solvable baselines by baseline length
    solve_ifrs = [];
    antpos = mssel.ms_antenna_positions;
    if (min_baseline or max_baseline) and antpos is not None:
      for (ip,p),(iq,q) in array.ifr_index():
        baseline = math.sqrt(((antpos[ip,:]-antpos[iq,:])**2).sum());
        if (not min_baseline or baseline > min_baseline) and \
           (not max_baseline or baseline < max_baseline):
          solve_ifrs.append((p,q));
    else:
      solve_ifrs = array.ifrs();
    inspect_ifrs = solve_ifrs;

  # make a predict tree using the MeqMaker
  if do_solve or do_subtract:
    predict = meqmaker.make_predict_tree(ns);
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

  # make nodes to compute residuals
  if do_subtract:
    residuals = ns.residuals;
    for p,q in array.ifrs():
      residuals(p,q) << spigots(p,q) - predict(p,q);
    outputs = residuals;

  # and now we may need to correct the outputs
  if do_correct:
    if do_correct_sky:
      srcs = meqmaker.get_source_list(ns);
      sky_correct = srcs and srcs[0];
    else:
      sky_correct = None;
    outputs = meqmaker.correct_uv_data(ns,outputs,sky_correct=sky_correct,inspect_ifrs=inspect_ifrs);

  # make solve trees
  if do_solve:
    # inputs to the solver are based on calibration type
    # if calibrating visibilities, feed them to condeq directly
    if cal_type == CAL.VIS:
      observed = spigots;
      model    = predict;
    # else take ampl/phase component
    else:
      model = ns.model;
      observed = ns.observed;
      if cal_type == CAL.AMPL:
        for p,q in array.ifrs():
          observed(p,q) << Meq.Abs(spigots(p,q));
          model(p,q)  << Meq.Abs(predict(p,q));
      elif cal_type == CAL.LOGAMPL:
        for p,q in array.ifrs():
          observed(p,q) << Meq.Log(Meq.Abs(spigots(p,q)));
          model(p,q)  << Meq.Log(Meq.Abs(predict(p,q)));
      elif cal_type == CAL.PHASE:
        for p,q in array.ifrs():
          observed(p,q) << 0;
          model(p,q)  << Meq.Abs(predict(p,q))*Meq.FMod(Meq.Arg(spigots(p,q))-Meq.Arg(predict(p,q)),2*math.pi);
      else:
        raise ValueError,"unknown cal_type setting: "+str(cal_type);
    # make a solve tree
    solve_tree = StdTrees.SolveTree(ns,model,solve_ifrs=solve_ifrs);
    # the output of the sequencer is either the residuals or the spigots,
    # according to what has been set above
    outputs = solve_tree.sequencers(inputs=observed,outputs=outputs);

  # make sinks and vdm.
  # The list of inspectors must be supplied here
  inspectors += meqmaker.get_inspectors() or [];
  StdTrees.make_sinks(ns,outputs,spigots=spigots0,post=inspectors);
  Bookmarks.make_node_folder("Corrected/residual visibilities by baseline",
    [ outputs(p,q) for p,q in array.ifrs() ],sorted=True,ncol=2,nrow=2);

  if not do_solve:
    if do_subtract:
      name = "Generate residuals";
      comment = "Generated residual visibilities.";
    elif do_correct:
      name = "Generate corrected data";
      comment = "Generated corrected visibilities.";
    else:
      name = None;
    if name:
      # make a TDL job to runsthe tree
      def run_tree (mqs,parent,**kw):
        global tile_size;
        purrpipe.title("Calibrating").comment(comment);
        mqs.execute(Meow.Context.vdm.name,mssel.create_io_request(tile_size),wait=False);
      TDLRuntimeMenu(name,
        TDLOption('tile_size',"Tile size, in timeslots",[10,60,120,240],more=int,
                  doc="""Input data is sliced by time, and processed in chunks (tiles) of
                  the indicated size. Larger tiles are faster, but use more memory."""),
        TDLRuntimeJob(run_tree,name)
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