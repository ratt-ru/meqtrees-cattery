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
define("STEFCAL_GAIN_PLOT_PREFIX","G","automatically plot gain (G) solutions. Set to empty string to disable.")
define("STEFCAL_GAIN1_PLOT_PREFIX","B","automatically plot gain1 (B) solutions. Set to empty string to disable.")
define("STEFCAL_IFRGAIN_PLOT_PREFIX","IG","automatically IFR gain solutions. Set to empty string to disable.")
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
              gain_plot_prefix="$STEFCAL_GAIN_PLOT_PREFIX",
              gain1_plot_prefix="$STEFCAL_GAIN1_PLOT_PREFIX",
              ifrgain_plot_prefix="$STEFCAL_IFRGAIN_PLOT_PREFIX",
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
  msname,section,lsm,label,plotvis,gain_plot_prefix,gain1_plot_prefix,ifrgain_plot_prefix,diffgain_plot_prefix = \
    interpolate_locals("msname section lsm label plotvis gain_plot_prefix gain1_plot_prefix ifrgain_plot_prefix diffgain_plot_prefix");
  
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
      if gain_plot_prefix:
        make_gain_plots(STEFCAL_GAIN_SAVE,prefix=gain_plot_prefix);
    if os.path.exists(STEFCAL_GAIN1):
      std.copy(STEFCAL_GAIN1,STEFCAL_GAIN1_SAVE);
      if gain_plot_prefix:
        make_gain_plots(STEFCAL_GAIN1_SAVE,prefix=gain_plot_prefix);
    if os.path.exists(STEFCAL_DIFFGAIN):
      std.copy(STEFCAL_DIFFGAIN,STEFCAL_DIFFGAIN_SAVE);
      if diffgain_plot_prefix:
        make_diffgain_plots(STEFCAL_DIFFGAIN_SAVE,prefix=diffgain_plot_prefix);
    if os.path.exists(STEFCAL_IFRGAIN):
      std.copy(STEFCAL_IFRGAIN,STEFCAL_IFRGAIN_SAVE);
      if ifrgain_plot_prefix:
        make_ifrgain_plots(STEFCAL_IFRGAIN_SAVE,prefix=ifrgain_plot_prefix);
                  
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



###################### PLOTTING ROUTINES

def _cmp_antenna (sa,sb):
  """Helper function to sort antenna names. Try numeric compare first, fall back to text compare if failed""";
  try:
    return cmp(int(sa),int(sb));
  except:
    return cmp(sa,sb);

import cPickle

define("FEED_TYPE","rl","feed type for plot labels: rl or xy")

define("DIFFGAIN_PLOT_DIR_Template","$OUTFILE.diffgain${SUFFIX}.plots","directory for diffgain plots")
define("DIFFGAIN_PLOT_AMPL_YLIM",(0,2),"explicit Y axis limits for amplitude plot, as (y1,y2) tuple, or None for auto-scale")

def make_diffgain_plots (filename="$STEFCAL_DIFFGAIN_SAVE",prefix="dE",dir="$DIFFGAIN_PLOT_DIR",ylim=None,subset=None,ant=None):
  """Makes a set of diffgain plots from the specified saved file ('filename'). Plots are placed into directory
  'dir' and filenames are prefixed by 'prefix'. 
  'ylim' can be used to override DIFFGAIN_PLOT_AMPL_YLIM setting.
  'subset' can be set to a whitespace-separated list of parameter names
  'ant' can be set to a whitespace-separated list of antennas. Wildcard patterns are allowed in both cases"""
  import pylab
  from Timba.Meq import meq
  
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
  antennas = sorted(DG[srcnames[0]]['solutions'].keys(),_cmp_antenna);
  
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
  # feed labels
  feeds = "RR","RL","LR","LL" if FEED_TYPE.upper() == "RL" else "XX","XY","YX","YY";  

  for ant in antennas:
    filename = II("$dir/$prefix-ant-$ant.png");   
    info("making plot $filename");
    pylab.figure(figsize=(5*ncols,3*8))
    for i,src in enumerate(sorted(srcnames)):
      sols = DG[src]['solutions'][ant]
      for j,xx in enumerate(feeds):
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
  
document_globals(make_diffgain_plots,"DIFFGAIN_PLOT* STEFCAL_DIFFGAIN_SAVE FEED_TYPE");


_IFRGAIN_DIR = "."
_IFRGAIN_PREFIX = "IG"
_IFRGAIN_TYPE = "rrll"

define("IFRGAIN_PLOT_Template","$OUTFILE.${_IFRGAIN_PREFIX}-${_IFRGAIN_TYPE}${-<_IFRGAIN_COMP}${-<SUFFIX}.png","filename for ifrgain plots")

import cmath

def _normifrgain (rr):
  """Converts gains to offsets with std""";
  if type(rr) in (float,complex):
    return abs(rr-1),0;
  else:
    offset = abs(rr[rr!=1]-1);
    return float(offset.mean()),float(offset.std());

def _complexifrgain (rr):
  """Converts gains to complex offsets with std""";
  if type(rr) in (float,complex):
    return rr-1,0;
  else:
    vals = rr[rr!=1]-1;
    offset = float(abs(vals).mean());
    mean = vals.mean().ravel()[0];
    mean = cmath.rect(offset,cmath.phase(mean));
    return mean,float((vals-mean).std());

def _is_unity (rr,ll):
  return type(rr) in (float,complex) and rr == 1 and type(ll) in (float,complex) and ll == 1;

def make_ifrgain_plots (filename="$STEFCAL_DIFFGAIN_SAVE",prefix="IG",feed="$IFRGAIN_PLOT_FEED"):
  """Makes a set of ifrgain plots from the specified saved file.""" 
  import pylab   
  from Timba.Meq import meq
  
  filename,prefix = interpolate_locals("filename prefix");
  global __IFRGAIN_PREFIX,_IFRGAIN_TYPE,_IFRGAIN_COMP;
  _IFRGAIN_PREFIX = prefix;
  dir = os.path.dirname(IFRGAIN_PLOT);
  if dir:
    makedir(dir)

  # feed labels
  feeds = "RR","LL","RL","LR" if FEED_TYPE.upper() == "RL" else "XX","YY","XY","YX";  

  ig = cPickle.load(file(filename))         

  def plot_xy (content,title):
    """Plots x vs y""";
    pylab.figure(figsize=(16,12))
    # plot errors bars, if available
    pylab.errorbar(
      [x for l,(x,xe),(y,ye) in content],[y for l,(x,xe),(y,ye) in content],
      [ye for l,(x,xe),(y,ye) in content],[xe for l,(x,xe),(y,ye) in content],
      fmt=None,ecolor="lightgrey"
    );
    # plot labels
    for label,(x,xe),(y,ye) in content:
      pylab.text(x,y,label,horizontalalignment='center',verticalalignment='center',size=8)
    pylab.title(title)
    pylab.savefig(II("$IFRGAIN_PLOT"),dpi=100);
    info("generated plot $IFRGAIN_PLOT")

  def plot_complex (content,title):
    """Plots x vs y""";
    pylab.figure(figsize=(16,12))
    # plot errors bars, if available
    pylab.errorbar(
      [x.real for l1,l2,(x,xe),(y,ye) in content],[x.imag for l1,l2,(x,xe),(y,ye) in content],
      [xe for l1,l2,(x,xe),(y,ye) in content],[xe for l1,l2,(x,xe),(y,ye) in content],
      fmt=None,ecolor="lightgrey"
    );
    # plot labels
    for l1,l2,(x,xe),(y,ye) in content:
      pylab.text(x.real,x.imag,l1,horizontalalignment='center',verticalalignment='center',color='blue',size=8)
      pylab.text(y.real,y.imag,l2,horizontalalignment='center',verticalalignment='center',color='red',size=8)
    pylab.title(title)
    pylab.savefig(II("$IFRGAIN_PLOT"),dpi=100);
    info("generated plot $IFRGAIN_PLOT")

  def plot_ants (content,title):
    """Plots x vs y""";
    pylab.figure(figsize=(16,12))
    # plot errors bars, if available
    # pylab.errorbar(
    #   range(1,len(content)+1),[x for lx,x,xe,ly,y,ye in content],
    #   [xe for lx,x,xe,ly,y,ye in content],
    #   fmt=None,ecolor="lightgrey"
    # );
    # plot labels
    for i,(p,gainlist) in enumerate(content):
      for label,color,(value,std) in gainlist: 
        if value:
          pylab.plot(i,value,'w.')
          pylab.text(i,value,label,horizontalalignment='center',verticalalignment='center',size=8,color=color)
    pylab.xlim(-1,len(content));
    pylab.xticks(range(len(content)),[p for p,gainlist in content ]);
    pylab.title(title)
    pylab.savefig(II("$IFRGAIN_PLOT"),dpi=100);
    info("generated plot $IFRGAIN_PLOT")
  
  # plot RR vs LL offsets
  crl = [ ("%s-%s"%(p,q),_normifrgain(rr),_normifrgain(ll)) 
    for (p,q),(rr,rl,lr,ll) in ig.iteritems() if not _is_unity(rr,ll) ];
  have_diag = bool(crl);
  if have_diag:
    _IFRGAIN_TYPE = ("%s%s-offset"%feeds[:2]).lower();
    plot_xy(crl,"IFR offset amplitude (%s vs. %s)"%feeds[:2]);

  # plot RL vs LR offsets, if present
  crl = [ ("%s-%s"%(p,q),_normifrgain(rl),_normifrgain(lr)) 
    for (p,q),(rr,rl,lr,ll) in ig.iteritems() if not _is_unity(rl,lr) ];
  have_offdiag = bool(crl);
  if have_offdiag:
    _IFRGAIN_TYPE = ("%s%s-offset"%feeds[:2]).lower();
    plot_xy(crl,"IFR offset amplitude (%s vs. %s)"%feeds[2:]);

  # plot RR and LL complex
  if have_diag:
    crl = [ ("%s-%s:%s"%(p,q,feeds[0].upper()),"%s-%s:%s"%(p,q,feeds[1].upper()),_complexifrgain(rr),_complexifrgain(ll)) 
      for (p,q),(rr,rl,lr,ll) in ig.iteritems() if not _is_unity(rr,ll) ];
    _IFRGAIN_TYPE = ("%s%s-complex"%feeds[:2]).lower();
    plot_complex(crl,"IFR complex %s %s offsets"%feeds[:2]);

  if have_offdiag:
    crl = [ ("%s-%s:%s"%(p,q,feeds[0].upper()),"%s-%s:%s"%(p,q,feeds[1].upper()),_complexifrgain(rl),_complexifrgain(lr)) 
      for (p,q),(rr,rl,lr,ll) in ig.iteritems() if not _is_unity(rl,lr) ];
    _IFRGAIN_TYPE = ("%s%s-complex"%feeds[2:]).lower();
    plot_complex(crl,"IFR complex %s %s offsets"%feeds[2:]);

  # plot per antenna
  if have_diag:
    igpa = {}
    for (p,q),(rr,rl,lr,ll) in ig.iteritems():
      rr,ll = _normifrgain(rr),_normifrgain(ll);
      igpa.setdefault(p,[]).append(("%s:%s"%(q,feeds[0]),'blue',rr));
      igpa.setdefault(q,[]).append(("%s:%s"%(p,feeds[0]),'blue',rr));
      igpa.setdefault(p,[]).append(("%s:%s"%(q,feeds[1]),'red',ll));
      igpa.setdefault(q,[]).append(("%s:%s"%(p,feeds[1]),'red',ll));
    content = [ (p,igpa[p]) for p in sorted(igpa.keys(),cmp=_cmp_antenna) ];
    _IFRGAIN_TYPE = ("%s%s-perant"%feeds[:2]).lower();
    plot_ants(content,"IFR %s %s offset amplitudes per antenna"%feeds[:2]);

  if have_offdiag:
    igpa = {}
    for (p,q),(rr,rl,lr,ll) in ig.iteritems():
      rr,ll = _normifrgain(rl),_normifrgain(lr);
      igpa.setdefault(p,[]).append(("%s:%s"%(q,feeds[2]),'blue',rr));
      igpa.setdefault(q,[]).append(("%s:%s"%(p,feeds[2]),'blue',rr));
      igpa.setdefault(p,[]).append(("%s:%s"%(q,feeds[3]),'red',ll));
      igpa.setdefault(q,[]).append(("%s:%s"%(p,feeds[3]),'red',ll));
    content = [ (p,igpa[p]) for p in sorted(igpa.keys(),cmp=_cmp_antenna) ];
    _IFRGAIN_TYPE = ("%s%s-perant"%feeds[:2]).lower();
    plot_ants(content,"IFR %s %s offset amplitudes per antenna"%feeds[2:]);

document_globals(make_ifrgain_plots,"IFRGAIN_PLOT* STEFCAL_IFRGAIN_SAVE FEED_TYPE");


_GAIN_PREFIX = "G"
_GAIN_TYPE = "summary"

define("GAIN_PLOT_Template","$OUTFILE.${_GAIN_PREFIX}-${_GAIN_TYPE}${-<SUFFIX}.png","filename for ifrgain plots")
define("GAIN_PLOT_AMPL_YLIM",(0,2),"explicit Y axis limits for amplitude plot, as (y1,y2) tuple, or None for auto-scale")

define("GAIN_PLOT_AMPL_STYLE",'dot',"amplitude plot style. Choose between 'fill' or 'dot'.")

def make_gain_plots (filename="$STEFCAL_GAIN_SAVE",prefix="G",ylim=None,ant=None):
  """Makes a set of gain plots from the specified saved file ('filename'). 
  'ylim' can be used to override GAIN_PLOT_AMPL_YLIM setting.
  'ant' can be set to a whitespace-separated list of antennas. Wildcard patterns are allowed."""
  import pylab
  from Timba.Meq import meq
  
  filename,prefix = interpolate_locals("filename prefix");

  global _GAIN_PREFIX,_GAIN_TYPE;
  _GAIN_PREFIX = prefix;

  info("loading gain solutions from $filename");
  G0 = cPickle.load(file(filename))

  G = G0['gains'][prefix]['solutions'];
  antennas = antennas0 = sorted(G.keys(),_cmp_antenna);
  info("loaded solutions for %d antennas"%len(antennas));
   
  ylim = ylim or GAIN_PLOT_AMPL_YLIM;

  if ant:
    antpatts = str(ant).split();
    antennas = [ x for x in antennas if any([ fnmatch.fnmatch(str(x),patt) for patt in antpatts ]) ];
    info("applied subset '$ant' to antennas %s, selecting %d"%(" ".join(map(str,antennas0)),len(antennas)));

  if not antennas:
    info("no antennas to plot");
    return;

  # feed labels
  feeds = "RR","RL","LR","LL" if FEED_TYPE.upper() == "RL" else "XX","XY","YX","YY";  

  ncols = 8
  nrows = len(antennas)

  # find amplitude axis ranges for parallel and cross-hand plots
  ppminmax = 1e+99,-1e+99
  xhminmax = 1e+99,-1e+99
  for ant in antennas: 
    for j in range(4):
      gg = G[ant][j];
      valid = (gg!=0)&(gg!=1);  # mask trivial or unfilled solutions
      if valid.any():
        a = abs(gg);
        amin,amax = a[valid].min(),a[valid].max();
        if j in (0,3):
          ppminmax = (min(ppminmax[0],amin),max(ppminmax[1],amax));
        else:
          xhminmax = (min(xhminmax[0],amin),max(xhminmax[1],amax));

  for xaxis,yaxis,label in (0,1,"time"),(1,0,"freq"):
    pylab.figure(figsize=(5*ncols,3*nrows))
    for row,ant in enumerate(antennas):
      sols = G[ant]
      for icol,j in enumerate([0,3,1,2]):
        feed = feeds[j];
        gg = numpy.transpose(sols[j],(xaxis,yaxis));
        # get plot axis and averaging axis
        nx,ny = gg.shape;
        x = numpy.zeros((nx,ny));
        x[...] = numpy.arange(nx)[:,numpy.newaxis];
        pylab.subplot(nrows,ncols,row*ncols+icol*2+1);
        pylab.subplots_adjust(bottom=.01,left=.01,top=.99,right=.99);
        valid = (gg!=0)&(gg!=1);  # mask trivial or unfilled solutions
        amp = abs(gg)
        ampwh = numpy.ma.masked_array(amp,mask=~valid);
        amid  = ampwh.mean(1);
        xvalid = ~amid.mask;
        if GAIN_PLOT_AMPL_STYLE == 'fill':
          pylab.fill_between(x[:,0],ampwh.min(1),ampwh.max(1),where=xvalid,color='grey')
          pylab.plot(x[xvalid,0],amid[xvalid],'.',mec='blue',mfc='blue')
        else:
          pylab.plot(x[valid],amp[valid],'.',ms=0.5,mec='grey',mfc='grey')
          pylab.plot(x[xvalid,0],amid[xvalid],'.-',ms=0.5,mec='blue',mfc='blue')
        pylab.title("%s:%s:ampl"%(ant,feed));
        pylab.xticks([]);
        pylab.xlim(0,nx-1)
        pylab.ylim(*(ppminmax if j in (0,3) else xhminmax))
        
        pylab.subplot(nrows,ncols,row*ncols+icol*2+2)
        ph0 = numpy.angle(gg)*180/math.pi;
        pylab.plot(x[valid],ph0[valid],'.',ms=0.5,mec='grey',mfc='grey')
        ph = ph0[:,ny/2];
        xvalid = valid[:,ny/2];
        if xvalid.any():
          pylab.plot(x[xvalid,0],ph[xvalid],'.-',ms=0.5,mec='blue',mfc='blue')
        
        pylab.title("%s:%s:phase (deg)"%(ant,feed));
        pylab.xticks([]);
        pylab.xlim(0,nx-1)

    _GAIN_TYPE = label;
    pylab.savefig(II("$GAIN_PLOT"),dpi=150);  
    info("generated plot $GAIN_PLOT")


  
document_globals(make_gain_plots,"GAIN_PLOT* STEFCAL_GAIN_SAVE FEED_TYPE");

