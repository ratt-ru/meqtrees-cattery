
 # standard preamble
from Timba.TDL import *
from Timba.Meq import meq
import math

import Meow
import Meow.StdTrees
from Meow import Context,ParmGroup,Bookmarks
from Meow.Parameterization import resolve_parameter

# MS options first
mssel = Context.mssel = Meow.MSUtils.MSSelector(has_input=True,tile_sizes=None,flags=False,hanning=True);
# MS compile-time options
TDLCompileOptions(*mssel.compile_options());
# UVW options
TDLCompileOptions(*Meow.IfrArray.compile_options());
# MS run-time options
TDLRuntimeMenu("MS/data selection options",*mssel.runtime_options());
## also possible:

CORRELATIONS = [ "XX","XY","YX","YY" ];
CORR_INDICES = dict([(corr,index) for index,corr in enumerate(CORRELATIONS)]);
ALL_CORRS = "XX XY YX YY";
DIAG_CORRS = "XX YY";
CROSS_CORRS = "YX XY";
SINGLE_CORR = "XX";

# output mode menu
TDLCompileMenu("What do we want to do",
  TDLMenu("Calibrate",
     TDLOption('cal_vis',"Calibrate visibilities",True),
     TDLOption('cal_ampl',"Calibrate amplitudes",False),
     TDLOption('cal_log_ampl',"Calibrate log-amplitudes",False),
     TDLOption('cal_phase',"Calibrate phases",False),
     TDLOption('cal_corr',"Use correlations",
                [ALL_CORRS,DIAG_CORRS,CROSS_CORRS,SINGLE_CORR]),
     toggle='do_solve',open=True,exclusive='solve_type'
  ),
  TDLOption('do_subtract',"Subtract sky model and generate residuals",True),
  TDLOption('do_correct',"Correct the data or residuals",True),
  TDLOption('do_invert_phase',"Invert phases in input data",True,
     doc="""Inverts phases in the input data. Some e.g. WSRT MSs require this."""),
);
do_correct_sky = False;
  #TDLOption('do_correct_sky',"...include sky-Jones correction for first source in model",True));

# now load optional modules for the ME maker
from Meow import MeqMaker
meqmaker = MeqMaker.MeqMaker(solvable=True);

# specify available sky models
# these will show up in the menu automatically
import Meow.LSM
lsm = Meow.LSM.MeowLSM(include_options=False);

meqmaker.add_sky_models([lsm]);

# now add optional Jones terms
# these will show up in the menu automatically


# B - bandpass, G - gain
import solvable_jones
meqmaker.add_uv_jones('B','bandpass',
  [ solvable_jones.DiagAmplPhase(),
    solvable_jones.FullRealImag() ]);
meqmaker.add_uv_jones('G','receiver gains/phases',
  [ solvable_jones.DiagAmplPhase(),
    solvable_jones.FullRealImag() ]);


# very important -- insert meqmaker's options properly
TDLCompileOptions(*meqmaker.compile_options());

# resampling option
TDLCompileMenu("Enable resampling",
  TDLOption('resample_time',"Resampling factor in time",[3,5,10],more=int),
  TDLOption('resample_freq',"Resampling factor in freq",[3,5,10],more=int),
  toggle='do_resample');

TDLCompileOption('do_ifr_errors',"Use interferometer-based errors",False);


def _define_forest (ns):
  ANTENNAS = mssel.get_antenna_set(range(1,15));
  array = Meow.IfrArray(ns,ANTENNAS,resamplers=do_resample);
  observation = Meow.Observation(ns);
  Meow.Context.set(array,observation);
  stas = array.stations();

  # make spigot nodes
  spigots = spigots0 = outputs = array.spigots();
  # invert phases if necessary
  if do_invert_phase:
    for p,q in array.ifrs():
      ns.conj_spigot(p,q) << Meq.Conj(spigots(p,q));
    spigots = ns.conj_spigot;

  # ...and an inspector for them
  Meow.StdTrees.vis_inspector(ns.inspector('input'),spigots,
                              bookmark="Inspect input visibilities");
  inspectors = [ ns.inspector('input') ];
  Bookmarks.make_node_folder("Input visibilities by baseline",
    [ spigots(p,q) for p,q in array.ifrs() ],sorted=True,ncol=2,nrow=2);

  # make a ParmGroup and solve jobs for source fluxes
  srcs = meqmaker.get_source_list(ns);
  if srcs and do_solve:
    # make dictionary of source parameters. Use a dict since
    # sources may share paramaters, so we don't want them in the list twice
    parms = {};
    for src in srcs:
      for p in src.visibilities().search(tags="source solvable"):
        parms[p.name] = p;
    if parms:
      pg_src = ParmGroup.ParmGroup("source",parms.values(),
                  table_name="sources.mep",
                  individual=True,bookmark=True);
      # now make a solvejobs for the source
      ParmGroup.SolveJob("cal_source","Calibrate source model",pg_src);

  # make list of selected correlations
  selected_corrs = cal_corr.split(" ");
  # make a predict tree using the MeqMaker
  if do_solve or do_subtract:
    predict = meqmaker.make_tree(ns);
    # add ifr errors if necessary
    if do_ifr_errors:
      ifr_gain = ns.ifr_gain;
      ifr_bias = ns.ifr_bias;
      def1 = Meow.Parm(1);
      def0 = Meow.Parm(0);
      for p,q in array.ifrs():
        for corr in CORRELATIONS:
          if corr in selected_corrs:
            ifr_gain(p,q,corr) << Meq.ToComplex(
                resolve_parameter(corr,ifr_gain(p,q,corr,'r'),def1,tags="ifr gain real"),
                resolve_parameter(corr,ifr_gain(p,q,corr,'i'),def0,tags="ifr gain imag"));
            ifr_bias(p,q,corr) << Meq.ToComplex(
                resolve_parameter(corr,ifr_bias(p,q,corr,'r'),def0,tags="ifr bias real"),
                resolve_parameter(corr,ifr_bias(p,q,corr,'i'),def0,tags="ifr bias imag"));
          else:
            ifr_gain(p,q,corr) << 1;
            ifr_bias(p,q,corr) << 0;
        ifr_gain(p,q) << \
            Meq.Matrix22(*[ifr_gain(p,q,corr) for corr in CORRELATIONS]);
        ifr_bias(p,q) << \
            Meq.Matrix22(*[ifr_bias(p,q,corr) for corr in CORRELATIONS]);
        ns.predict_ifr(p,q) << predict(p,q)*ifr_gain(p,q) + ifr_bias(p,q);
      predict = ns.predict_ifr;
      # add parmgroups
      pg_ifr_ampl = ParmGroup.ParmGroup("ifr_gain",
                       ifr_gain.search(tags="solvable ifr gain"),
                       table_name="ifr_error_gain.mep",bookmark=False);
      pg_ifr_bias  = ParmGroup.ParmGroup("ifr_bias",
                       ifr_bias.search(tags="solvable ifr bias"),
                       table_name="ifr_error_bias.mep",bookmark=False);
      Bookmarks.make_node_folder("Interferometer-based gains",
        [ ifr_gain(p,q) for p,q in array.ifrs() ],sorted=True,nrow=2,ncol=2);
      Bookmarks.make_node_folder("Interferometer-based biases",
        [ ifr_bias(p,q) for p,q in array.ifrs() ],sorted=True,nrow=2,ncol=2);
      ParmGroup.SolveJob("cal_ifr_ampl",
                          "Calibrate interferometer-based gains",pg_ifr_ampl);
      ParmGroup.SolveJob("cal_ifr_bias",
                          "Calibrate interferometer-based biases",pg_ifr_bias);
    # add resampling if needed
    if do_resample:
      for p,q in array.ifrs():
        modres = ns.modres(p,q) << Meq.ModRes(predict(p,q),
                                              upsample=[resample_time,resample_freq]);
        ns.resampled(p,q) << Meq.Resampler(modres,mode=2);
      predict = ns.resampled;
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
    outputs = meqmaker.correct_uv_data(ns,outputs,sky_correct=sky_correct);
    if do_ifr_errors:
      for p,q in array.ifrs():
        ns.correct_ifr(p,q) << (outputs(p,q) - ns.ifr_bias(p,q))/ns.ifr_gain(p,q);
      outputs = ns.correct_ifr;

  # make solve trees
  if do_solve:
    # extract selected correlations
    if cal_corr != ALL_CORRS:
      index = [ CORR_INDICES[c] for c in selected_corrs ];
      for p,q in array.ifrs():
        ns.sel_predict(p,q) << Meq.Selector(predict(p,q),index=index,multi=True);
        ns.sel_spigot(p,q)  << Meq.Selector(spigots(p,q),index=index,multi=True);
      spigots = ns.sel_spigot;
      predict = ns.sel_predict;
    # inputs to the solver are based on calibration type
    # if calibrating visibilities, feed them to condeq directly
    if solve_type == 'cal_vis':
      observed = spigots;
      model    = predict;
    # else take ampl/phase component
    else:
      model = ns.model;
      observed = ns.observed;
      if solve_type == 'cal_ampl':
        for p,q in array.ifrs():
          observed(p,q) << Meq.Abs(spigots(p,q));
          model(p,q)  << Meq.Abs(predict(p,q));
      elif solve_type == 'cal_log_ampl':
        for p,q in array.ifrs():
          observed(p,q) << Meq.Log(Meq.Abs(spigots(p,q)));
          model(p,q)  << Meq.Log(Meq.Abs(predict(p,q)));
      elif solve_type == 'cal_phase':
        for p,q in array.ifrs():
          observed(p,q) << 0;
          model(p,q)  << Meq.Abs(predict(p,q))*Meq.FMod(Meq.Arg(spigots(p,q))-Meq.Arg(predict(p,q)),2*math.pi);
      else:
        raise ValueError,"unknown solve_type setting: "+str(solve_type);
    # make a solve tree
    solve_tree = Meow.StdTrees.SolveTree(ns,model);
    # the output of the sequencer is either the residuals or the spigots,
    # according to what has been set above
    outputs = solve_tree.sequencers(inputs=observed,outputs=outputs);

  # make sinks and vdm.
  # The list of inspectors must be supplied here
  inspectors += meqmaker.get_inspectors() or [];
  Meow.StdTrees.make_sinks(ns,outputs,spigots=spigots0,post=inspectors);
  Bookmarks.make_node_folder("Corrected/residual visibilities by baseline",
    [ outputs(p,q) for p,q in array.ifrs() ],sorted=True,ncol=2,nrow=2);

  if not do_solve:
    global _run_tree;
    if do_subtract:
      TDLRuntimeJob(_run_tree,"Generate residuals");
    elif do_correct:
      TDLRuntimeJob(_run_tree,"Generate corrected data");

  # very important -- insert meqmaker's runtime options properly
  # this should come last, since runtime options may be built up during compilation.
  TDLRuntimeOptions(*meqmaker.runtime_options(nest=False));
  # and insert all solvejobs
  TDLRuntimeOptions(*ParmGroup.get_solvejob_options());
  # finally, setup imaging options
  imsel = mssel.imaging_selector(npix=512,arcmin=meqmaker.estimate_image_size());
  TDLRuntimeMenu("Imaging options",*imsel.option_list());


# runs the tree
def _run_tree (mqs,parent,tiling=240,**kw):
  mqs.execute(Context.vdm.name,mssel.create_io_request(tiling),wait=False);


