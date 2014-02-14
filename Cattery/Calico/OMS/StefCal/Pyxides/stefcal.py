from Pyxis.ModSupport import *

import std,ms,lsm,mqt,imager

# register ourselves with Pyxis, and define what superglobals we use (these come from ms)
register_pyxis_module(superglobals="MS LSM DESTDIR OUTFILE STEP");

define("STEFCAL_SCRIPT","${mqt.CATTERY}/Calico/calico-stefcal.py","stefcal TDL script")
define("STEFCAL_SECTION","stefcal","default TDL config section")
define("STEFCAL_JOBNAME","stefcal","default TDL job name")
define("STEFCAL_TDLOPTS","","extra TDL options for stefcal")
define("STEFCAL_GAIN_Template","$MS/gain${SUFFIX}.cp","current file for gain solutions")
define("STEFCAL_GAIN1_Template","$MS/gain1${SUFFIX}.cp","current file for gain1 solutions")
define("STEFCAL_IFRGAIN_Template","$MS/ifrgain${SUFFIX}.cp","current file for IFR gain solutions")
define("STEFCAL_DIFFGAIN_Template","$MS/diffgain${SUFFIX}.cp","current file for diffgain solutions")
define("STEFCAL_DIFFGAIN1_Template","$MS/diffgain${SUFFIX}.cp","current file for diffgain solutions")
define("STEFCAL_GAIN_SAVE_Template","$OUTFILE.gain${SUFFIX}.cp","archive destination for gain solutions")
define("STEFCAL_GAIN1_SAVE_Template","$OUTFILE.gain1${SUFFIX}.cp","archive destination for gain1 solutions")
define("STEFCAL_DIFFGAIN_SAVE_Template","$OUTFILE.diffgain${SUFFIX}.cp","archive destination for diffgain solutions")
define("STEFCAL_IFRGAIN_SAVE_Template","$OUTFILE.ifrgain${SUFFIX}.cp","archive destination for IFR gain solutions")
define("STEFCAL_DIFFGAIN_SMOOTHING","","smoothing kernel (time,freq) for diffgains, overrides TDL config file")
define("STEFCAL_DIFFGAIN_INTERVALS","","solution intervals (time,freq) for diffgains, overrides TDL config file")

define("STEFCAL_STEP_INCR",1,"automatically increment v.STEP with each call to stefcal");

def stefcal ( msname="$MS",section="$STEFCAL_SECTION",
              diffgains=None,
              apply_only=False,
              gain_apply_only=False,
              diffgain_apply_only=False,
              ifrgain_apply_only=False,
              diffgain_intervals=None,diffgain_smoothing=None,
              flag_threshold=None,
              output="CORR_RES",
              plotvis="${ms.PLOTVIS}",
              dirty=True,restore=False,restore_lsm=True,
              label=None,
              args=[],options={},
              **kws):
  """Generic function to run a stefcal job.
  
  'section'         TDL config file section
  'label'           will be assigned to the global LABEL for purposes of file naming
  'apply_only'      if true, will only apply saved solutions
  'diffgains'       set to a source subset string to solve for diffgains. Set to True to use "=dE"
  'diffgain_mode'   'solve-save' to solve & save, 'solve-nosave' to not save, 'apply' to apply only
  'flag_threshold'  threshold flaging post-solutions. Give one threshold to flag with --above,
                    or T1,T2 for --above T1 --fm-above T2
  'output'          output visibilities ('CORR_DATA','CORR_RES', 'RES' are useful)
  'plotvis'         if not empty, specifies which output visibilities to plot using plot-ms (see plot.ms.py --help) 
  'dirty','restore' 
  'restore_lsm'     image output visibilities (passed to imager.make_image above as is)
  'args','options'  passed to the stefcal job as is (as a list of arguments and kw=value pairs), 
                    can be used to supply extra TDL options
  extra keywords:   passed to the stefcal job as kw=value, can be used to supply extra TDL options, thus
                    overriding settings in the TDL config file. Useful arguments of this kind are e.g.:
                    stefcal_reset_all=True to remove prior gains solutions.
  """
  msname,section,lsm,label,plotvis = interpolate_locals("msname section lsm label plotvis");
  
  makedir(v.DESTDIR);
  
  # increment step counter and assign global label
  
  if label is not None:
    v.LABEL = str(label);
  if type(v.STEP) is int and STEFCAL_STEP_INCR:
    v.STEP += STEFCAL_STEP_INCR;

  # setup stefcal options and run 
  info("Running stefcal ${step <STEP} ${(<LABEL>)}");
  # setup args
  args0 = [ """${ms.MS_TDL} ${ms.CHAN_TDL} ${lsm.LSM_TDL} ms_sel.ms_ifr_subset_str=${ms.IFRS} 
    stefcal_gain.enabled=1 stefcal_diffgain.enabled=%d %s"""%
    ((1 if diffgains else 0),STEFCAL_TDLOPTS) ];
  if diffgains:
    if diffgains is True:
      diffgains = "=dE";
    args0.append("de_subset.subset_enabled=1 de_subset.source_subset=$diffgains"); 
  opts = {
    'do_output': output,
    'stefcal_gain.mode': "apply" if apply_only or gain_apply_only else "solve-save",
    'stefcal_gain1.mode': "apply" if apply_only  or gain_apply_only else "solve-save",
    'stefcal_diffgain.mode': "apply" if apply_only or diffgain_apply_only else "solve-save",
    'stefcal_diffgain1.mode': "apply" if apply_only or diffgain_apply_only  else "solve-save",
    'stefcal_ifr_gain_mode': "apply" if apply_only or ifrgain_apply_only else "solve-save",
    'stefcal_gain.table': STEFCAL_GAIN,
    'stefcal_gain1.table': STEFCAL_GAIN1,
    'stefcal_diffgain.table': STEFCAL_DIFFGAIN,
    'stefcal_diffgain1.table': STEFCAL_DIFFGAIN1,
    'stefcal_ifr_gain_table': STEFCAL_IFRGAIN,
    'stefcal_visualize': False
  }
  timesmooth,freqsmooth = diffgain_smoothing or STEFCAL_DIFFGAIN_SMOOTHING or (0,0);
  timeint,freqint = diffgain_intervals or STEFCAL_DIFFGAIN_INTERVALS or (0,0);
  opts['stefcal_diffgain.timeint'] = 0 if timesmooth else timeint;
  opts['stefcal_diffgain.freqint'] = 0 if freqsmooth else freqint;
  opts['stefcal_diffgain.timesmooth'] = timesmooth;
  opts['stefcal_diffgain.freqsmooth'] = freqsmooth;

  # add user-defined args
  args0 += list(args);
  opts.update(options);
  opts.update(kws);
  # run the job
  mqt.run(STEFCAL_SCRIPT,STEFCAL_JOBNAME,section=section,args=args0,options=opts);
  
  # copy gains
  if not apply_only:
    if os.path.exists(STEFCAL_GAIN):
      std.copy(STEFCAL_GAIN,STEFCAL_GAIN_SAVE);
    if os.path.exists(STEFCAL_GAIN1):
      std.copy(STEFCAL_GAIN1,STEFCAL_GAIN1_SAVE);
    if os.path.exists(STEFCAL_DIFFGAIN):
      std.copy(STEFCAL_DIFFGAIN,STEFCAL_DIFFGAIN_SAVE);
    if os.path.exists(STEFCAL_IFRGAIN):
      std.copy(STEFCAL_IFRGAIN,STEFCAL_IFRGAIN_SAVE);
    
  # post-calibration flagging
  if flag_threshold:
    if isinstance(flag_threshold,(list,tuple)):
      t0,t1 = flag_threshold;
    else:
      t0,t1 = flag_threshold,None;
    ms.flagms("--above %g"%t0,"-f threshold -c");
    if t1:
      ms.flagms("--fm-above %g"%t1,"-f fmthreshold -c");

  # plot residuals
  if plotvis:
    info("Plotting visibilities ($plotvis)");
    ms.PLOTVIS = plotvis;
    ms.plotms("-o ${OUTFILE}_${output}${_s<STEP}${_<label}.png");
    
  # make images
  imager.make_image(msname,dirty=dirty,restore=restore,restore_lsm=restore_lsm);

# document global options for stefcal()
document_globals(stefcal,"MS LSM mqt.TDLCONFIG STEFCAL_* ms.DDID ms.CHANRANGE ms.IFRS ms.PLOTVIS STEP LABEL");
  