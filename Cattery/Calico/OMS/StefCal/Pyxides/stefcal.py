from Pyxis.ModSupport import *

import std,ms,lsm,mqt,imager

# register ourselves with Pyxis, and define what superglobals we use (these come from ms)
register_pyxis_module(superglobals="MS LSM DESTDIR OUTFILE STEP");

define("STEFCAL_SCRIPT","${mqt.CATTERY}/Calico/calico-stefcal.py","stefcal TDL script")
define("STEFCAL_SECTION","stefcal","default TDL config section")
define("STEFCAL_JOBNAME","stefcal","default TDL job name")
define("STEFCAL_TDLOPTS","","extra TDL options for stefcal")
define("STEFCAL_OUTPUT_COLUMN","CORRECTED_DATA","default output column")
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
define("STEFCAL_DIFFGAIN_PLOT_PREFIX","dE","automatically plot diffgain solutions. Set to empty string to disable.")

define("STEFCAL_STEP_INCR",1,"automatically increment v.STEP with each call to stefcal");

def stefcal ( msname="$MS",section="$STEFCAL_SECTION",
              diffgains=None,
              apply_only=False,
              reset=False,
              gain_apply_only=False,
              gain_reset=False,
              diffgain_apply_only=False,
              diffgain_reset=False,
              diffgain_plot_prefix="$STEFCAL_DIFFGAIN_PLOT_PREFIX",
              ifrgain_apply_only=False,
              ifrgain_reset=False,
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
  'apply_only'      if true, will only apply saved solutions rather than re-solve
  '{gain,diffgain,ifrgain}_apply_only'      
                    if true, will only apply saved gain/diffgain/IFR gain solutions rather than re-solve
  'reset'           if true, will reset all saved solutuions prior to starting
  '{gain,diffgain,ifrgain}_reset'      
                    if true, will reset saved gain/diffgain/IFR gain solutuions prior to starting
  'diffgains'       set to a source subset string to solve for diffgains. Set to True to use "=dE"
  'diffgain_mode'   'solve-save' to solve & save, 'solve-nosave' to not save, 'apply' to apply only
  'diffgain_plot'   automatically invoke make_diffgain_plots() if True
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
  msname,section,lsm,label,plotvis,diffgain_plot_prefix = \
    interpolate_locals("msname section lsm label plotvis diffgain_plot_prefix");
  
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
    ms_sel.output_column=$STEFCAL_OUTPUT_COLUMN
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
    'stefcal_gain.reset': int(reset or gain_reset),
    'stefcal_diffgain.mode': "apply" if apply_only or diffgain_apply_only else "solve-save",
    'stefcal_diffgain1.mode': "apply" if apply_only or diffgain_apply_only  else "solve-save",
    'stefcal_diffgain.reset': int(reset or diffgain_reset),
    'stefcal_diffgain1.reset': int(reset or diffgain_reset),
    'stefcal_ifr_gain_mode': "apply" if apply_only or ifrgain_apply_only else "solve-save",
    'stefcal_ifr_gain_reset': int(reset or ifrgain_reset),
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
      if diffgain_plot_prefix:
        make_diffgain_plots(STEFCAL_DIFFGAIN_SAVE,prefix=diffgain_plot_prefix);
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
  imager.make_image(msname,column=STEFCAL_OUTPUT_COLUMN,dirty=dirty,restore=restore,restore_lsm=restore_lsm);

# document global options for stefcal()
document_globals(stefcal,"MS LSM mqt.TDLCONFIG STEFCAL_* ms.DDID ms.CHANRANGE ms.IFRS ms.PLOTVIS STEP LABEL");

import cPickle

define("DIFFGAIN_PLOT_DIR_Template","$OUTFILE.diffgain${SUFFIX}.plots","directory for diffgain plots")
define("DIFFGAIN_PLOT_AMPL_YLIM",(0,2),"explicit Y axis limits for amplitude plot, as (y1,y2) tuple, or None for auto-scale")

def make_diffgain_plots (filename="$STEFCAL_DIFFGAIN_SAVE",prefix="dE",dir="$DIFFGAIN_PLOT_DIR",ylim=None,subset=None,ant=None):
  """Makes a set of diffgain plots from the specified saved file ('filename'). Plots are placed into directory
  'dir' and filenames are prefixed by 'prefix'. 
  'ylim' can be used to override DIFFGAIN_PLOT_AMPL_YLIM setting.
  'subset' can be set to a whitespace-separated list of parameter names
  'ant' can be set to a whitespace-separated list of antennas. Wildcard patterns are allowed in both cases"""
  import pylab
  
  filename,prefix,dir = interpolate_locals("filename prefix dir");
  if dir:
    makedir(dir);
  else:
    dir = ".";

  info("loading diffgain solutions from $filename");
  DG0 = cPickle.load(file(filename))

  DG = DG0['gains']
  srcnames = sorted(DG.keys())
  ncols = len(srcnames)
  antennas = sorted(DG[srcnames[0]]['solutions'].keys(),lambda a,b:cmp(int(a),int(b)));
  
  ylim = ylim or DIFFGAIN_PLOT_AMPL_YLIM;
  
  # apply subsets
  if subset:
    info("applying subset '$src' to parameters",*srcnames);
    src = subset.split();
    srcnames = [ x for x in srcnames if any([ fnmatch.fnmatch(x,patt) for patt in src ]) ];
  if ant:
    info("applying subset '$ant' to antennas",*antennas);
    ant = ant.split();
    antennas = [ x for x in antennas if any([ fnmatch.fnmatch(x,patt) for patt in ant ]) ];
    
  if not srcnames or not antennas:
    info("no parameters or antennas to plot");
    return;

  info("making diffgain plots for",*srcnames);
  info("and %d antennas"%len(antennas));

  for ant in antennas:
    filename = II("$dir/$prefix-ant-$ant.png");   
    info("making plot $filename");
    pylab.figure(figsize=(5*ncols,3*8))
    for i,src in enumerate(sorted(srcnames)):
      sols = DG[src]['solutions'][ant]
      for j,xx in enumerate(("RR","RL","LR","LL")):
        pylab.subplot(8,ncols,j*2*ncols+i+1)
        pylab.title("%s:%s:%s:ampl"%(src,ant,xx));
        if not hasattr(sols[j],'shape'):
          continue;
        ntime,nfreq = sols[j].shape;
        pylab.xticks([]);
        pylab.xlim(0,ntime-1);
        if ylim:
          pylab.ylim(*ylim);
        if nfreq>1:
          x = range(ntime);
          pylab.fill_between(x,abs(sols[j][:,0]),abs(sols[j][:,-1]),color='grey')
          pylab.plot(abs(sols[j][:,nfreq/2]))
        else:
          pylab.plot(abs(sols[j][:,0]))
        pylab.subplot(8,ncols,(j*2+1)*ncols+i+1)
        ph = numpy.angle(sols[j])*180/math.pi;
        if nfreq>1:
          pylab.plot(ph[:,0],'.',ms=0.5,mec='0.2')
          pylab.plot(ph[:,-1],'.',ms=0.5,mec='0.2')
          pylab.plot(ph[:,nfreq/2],'.b',ms=0.5)
        else:
          pylab.plot(ph[:,0],'.',ms=0.5)
        pylab.title("%s:%s:%s:phase (deg)"%(src,ant,xx));
        pylab.xticks([]);
        pylab.xlim(0,ntime-1)
    pylab.savefig(filename,dpi=150);
    pylab.close();

  ncols = len(srcnames)
  nrows = len(antennas)
  filename = II("$dir/$prefix-ampl-summary.png")
  info("making plot $filename");
  pylab.figure(figsize=(5*ncols,3*nrows));
  for row,ant in enumerate(antennas):
    for isrc,src in enumerate(sorted(srcnames)):
      sols = DG[src]['solutions'][ant][0]
      pylab.subplot(nrows,ncols,row*ncols+isrc+1)
      pylab.title("%s:%s:%s:ampl"%(src,ant,"RR"));
      if not hasattr(sols,'shape'):
        continue;
      ntime,nfreq = sols.shape;
      if nfreq>1:
        x = range(ntime);
        pylab.fill_between(x,abs(sols[:,0]),abs(sols[:,-1]),color='grey')
        pylab.plot(abs(sols[:,nfreq/2]))
      else:
        pylab.plot(abs(sols[:,0]))
      pylab.xticks([]);
      pylab.xlim(0,ntime-1)
      if ylim:
        pylab.ylim(*ylim);
  pylab.savefig(filename,dpi=150);
  pylab.close();
  
document_globals(make_diffgain_plots,"DIFFGAIN_PLOT* STEFCAL_DIFFGAIN_SAVE");