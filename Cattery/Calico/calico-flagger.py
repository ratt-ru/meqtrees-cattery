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

from Timba.TDL import *
from Timba.Meq import meq
from Timba.array import *
import os
import time
import signal

import Meow

from Meow import Bookmarks,Context
import Meow.StdTrees
import Calico.Flagger

# MS options first
mssel = Context.mssel = Meow.MSUtils.MSSelector(has_input=True,read_flags=True,write_flags=True,has_output=False,tile_sizes=[10,100,200]);\
mssel.enable_write_flags(False);  # enabled in define_forest() if needed
# MS compile-time options
TDLCompileOptions(*mssel.compile_options());
TDLCompileOption("run_purr","Start Purr on this MS",True);

# add a subset selector for the "flag subset" option
flag_subset = mssel.make_subset_selector('mssel_flag');

TDLCompileMenu("Enable in-tree flaggers",
  TDLOption("tree_flag_norm","Use matrix norm rather than element amplitudes",False),
  TDLOption("tree_flag_ampl","Flag on amplitudes >",[None],more=float,
      doc="""<P>If selected, your tree will flag visibility points where the complex amplitude exceeds the given value.</P>
      """),
  TDLOption("tree_flag_mean_ampl","Flag on mean amplitudes (over all IFRs) >",[None],more=float,
      doc="""<P>If selected, your tree will flag visibility points where the mean complex amplitude over all IFRs exceeds the given value.</P>
      """),
  TDLOption("tree_flag_mean_freq_ampl","Flag on mean amplitudes (over all IFRs and freqs) >",[None],more=float,
      doc="""<P>If selected, your tree will flag visibility points where the mean complex amplitude over all IFRs and frequencies exceeds the given value.</P>
      """),
  toggle='tree_flag_enable',
);

# The view_ms job simply cycles through the measurement set, without writing flags
def view_ms_job (mqs,parent,wait=False,**kw):
  req = mssel.create_io_request(write_flags=False);
  mqs.execute('VisDataMux',req,wait=wait);

# The flag_ms job cycles through the measurement set and generates flags
def treeflag_ms_job (mqs,parent,wait=False,**kw):
  req = mssel.create_io_request();
  mqs.execute('VisDataMux',req,wait=wait);

# Handle SIGCHLD signals (to manage external autoflag processes)
old_sigchld_handler = None;
sigchld_handler_installed = False;
glish_pid = None;
autoflag_msgbox = None;
progress_dialog = None;

def init_progress_dialog (parent,label):
  global progress_dialog;
  if progress_dialog is None:
    progress_dialog = GUI.ProgressDialog(label);
  else:
    progress_dialog.setLabelText(label);
  progress_dialog.show();
    
def progress_callback (current,total):
  if progress_dialog:
    progress_dialog.setMaximum(total);
    progress_dialog.setValue(current);

def _sigchild_handler (signal,frame):
  # check if glish process has exited
  global glish_pid;
  if glish_pid:
    try:
      pid,stat = os.waitpid(glish_pid,os.WNOHANG);
    except:
      # glish is already dead, but someone else reaped it for us...
      pid,stat = 0,0;
    if pid: 
      glish_pid = None;
      # re-enable the relevant options
      for opt in ms_job_options:
        opt.enable();
      # update status dialog, if present
      if autoflag_msgbox:
        # stat!=0: error, update dialog text
        if stat:
          exitcode = (stat&0x8F00)>>8;
          signal   = stat&0x00FF;
          if signal:
            status = "killed by signal %d"%signal;
          else:
            status = "exit code %d"%exitcode;
          autoflag_msgbox.setText("""<P>Oops, it appears the autoflag process has not
                exited cleanly (%s). This is usually a sign of something gone wrong. Please
                refer to the text console for more information.</P>
                """%status);
          autoflag_msgbox.setBoxType(GUI.Warning);
          autoflag_msgbox.show();
        # stat==0: normal exit
        else:
          autoflag_msgbox.setText("""<P>Autoflagging complete.</P>
                      <P>More information may be available on the text console.</P>""");
          autoflag_msgbox.setBoxType(GUI.Information);
  # call previous handler
  if callable(old_sigchld_handler):
    old_sigchld_handler(signal,frame);
  
# The run_autoflagger job is a wrapper around glish+autoflag tool
def run_autoflagger (mqs,parent,wait=False,**kw):
  af = flagger.autoflagger();
  # get column ID in autoflag terms
  column = dict(DATA='DATA',MODEL_DATA='MODEL',CORRECTED_DATA='CORR')[autoflag_column];
  # setup call to setdata()
  if autoflag_setdata:
    field = (autoflag_setdata_fieldid != "all" and autoflag_setdata_fieldid) or None;
    # channels is "all" (in which case we set "None" to use default), or "first,last" or "first last"
    start,nchan = None,None;
    if autoflag_setdata_chan != "all":
#      try:
        chan = map(int,re.split("[ ,]",autoflag_setdata_chan,1));
        start,nchan = chan[0],chan[1]-chan[0]+1;
 #     except:
        pass;
    af.setdata(spwid=autoflag_setdata_spwid,field=field,start=start,nchan=nchan,msselect=autoflag_setdata_msselect);
  # setup calls to individual flagging methods
  if autoflag_timemed:
    af.settimemed(thr=autoflag_timemed_thr,hw=autoflag_timemed_hw,
                  norow=not autoflag_timemed_flagrow,
                  rowhw=autoflag_timemed_rowhw,rowthr=autoflag_timemed_rowthr,
                  expr=autoflag_timemed_expr,fignore=autoflag_timemed_fignore,column=column);
  if autoflag_freqmed:
    af.setfreqmed(thr=autoflag_freqmed_thr,hw=autoflag_freqmed_hw,
                  norow=not autoflag_freqmed_flagrow,
                  rowhw=autoflag_freqmed_rowhw,rowthr=autoflag_freqmed_rowthr,
                  expr=autoflag_freqmed_expr,fignore=autoflag_freqmed_fignore,column=column);
  if autoflag_newtimemed:
    af.setnewtimemed(thr=autoflag_newtimemed_thr,
                  expr=autoflag_newtimemed_expr,fignore=autoflag_newtimemed_fignore,column=column);
  if autoflag_sprej:
    spwid = (autoflag_sprej_spwid != "all" and autoflag_sprej_spwid) or None;
    af.setsprej(rowthr=autoflag_sprej_rowthr,rowhw=autoflag_sprej_rowhw,
                ndeg=autoflag_sprej_ndeg,spwid=spwid,chan=[autoflag_sprej_chan0,autoflag_sprej_chan1],
                expr=autoflag_sprej_expr,fignore=autoflag_sprej_fignore,column=column);
  if autoflag_uvbin:
    plotchan = autoflag_uvbin_plotchan;
    if plotchan is None:
      plotchan = False;
    elif plotchan == "middle":
      plotchan = True;
    af.setuvbin(thr=autoflag_uvbin_thr,minpop=autoflag_uvbin_minpop,
                nbins=[autoflag_uvbin_nbin_uv,autoflag_uvbin_nbin_val],
                plotchan=plotchan,econoplot=autoflag_uvbin_econoplot,
                expr=autoflag_uvbin_expr,fignore=autoflag_uvbin_fignore,column=column);
  ## PGPlot windows don't play well with the browser, so disabling this for now
  # if autoflag_plotscr:
  #   plotscr = [autoflag_plotscr_nx,autoflag_plotscr_ny];
  # else:
  #   plotscr = False;
  if autoflag_plotdev:
    plotdev = [autoflag_plotdev_nx,autoflag_plotdev_ny];
  else:
    plotdev = False;
  # install a SIGCHLD handler if needed
  global sigchld_handler_installed;
  if not sigchld_handler_installed:
    old_sigchld_handler = signal.signal(signal.SIGCHLD,_sigchild_handler);
    sigchld_handler_installed = True;
  # pop up status window
  global autoflag_msgbox;
  if not autoflag_msgbox:
    autoflag_msgbox = GUI.MessageBox("Running autoflag","",GUI.Information,GUI.Button.Ok,GUI.Button.Ok);
    autoflag_msgbox.setButtonText(GUI.Button.Ok,"Hide this window");
  autoflag_msgbox.setText("""<P>We're running the autoflag process now. Please refrain from 
        any other operations on this MS until this window tells you that the process is complete.</P>
        <P>More information may be available on the text console.</P>
        """);
  autoflag_msgbox.setBoxType(GUI.Information);
  autoflag_msgbox.show();
  # form up command file
  cmdfile = time.strftime("autoflag-%m%d-%M%S.g",time.localtime(time.time()));
  cmdfile = os.path.join(mssel.msname,cmdfile);
  # save commands
  if autoflag_save:
    af.save(autoflag_save);
  # run the autoflagger
  for opt in ms_job_options:
    opt.disable();
  global glish_pid;
  glish_pid = af.run(plotdev=plotdev,devfile=autoflag_plotdev_file,         # plotscr=plotscr,
                     reset=autoflag_reset,trial=autoflag_trial,cmdfile=cmdfile,wait=wait);


def add_bitflags (mqs,parent,**kw):
  if not flagger:
    raise RuntimeError,"Flagger not available, perhaps MS was not specified?";
  flagger.add_bitflags();
  mssel.reload();

# The flag_ms job uses Calico.Flagger to flag subsets of the MS
def flag_ms (mqs,parent,**kw):
  if not flagger:
    raise RuntimeError,"Flagger not available, perhaps MS was not specified?";
  arg = dict();
  # form up list of antennas
  if flag_antennas:
    arg['antennas'] = Meow.MSUtils.parse_antenna_subset(flag_antennas,mssel.ms_antenna_names);
  # form up list of baselines
  if flag_baselines:
    arg['baselines'] = baselines = [];
    for spec in re.split("[\s,]+",flag_baselines):
      spec1,spec2 = spec.split('-',1);
      index1 = Meow.MSUtils.parse_antenna_subset(spec1,mssel.ms_antenna_names);
      index2 = Meow.MSUtils.parse_antenna_subset(spec2,mssel.ms_antenna_names);
      if index1 > index2:
        raise ValueError,"Illegal baseline specifier '%s'"%spec;
      baselines.append((index1[0],index2[0]));
  # form up time range
  if flag_reltime:
    arg['reltime'] = (flag_time0,flag_time1);
  # form up field, ddid & channel selection
  arg['ddid'] = flag_subset.get_ddid();
  arg['fieldid'] = flag_subset.get_field();
  arg['taql'] = flag_subset.get_taql_string();
  chans = flag_subset.get_channels();
  if chans:
    chan0,chan1,chanstep = chans;
    arg['channels'] = slice(chan0,chan1+1,chanstep);
  # form up flag/unflag arguments
  if flag_action == 'flag':
    arg['flag'] = flag_flag_selector.get_flagmask();
  else:
    arg['unflag'] = flag_flag_selector.get_flagmask();
  # open progress meter if GUI available
  init_progress_dialog(parent,"Flagging");
  try:
    # execute flagging
    stats = flagger.flag(progress_callback=progress_callback,**arg);
    flagger.close();
  finally:
    progress_callback(100,100);
    progress_dialog.hide();
  GUI.info_box("Flag statistics",
    """<P>%.4f%% of rows are flagged.</P>
      <P>%.4f%% of individual correlations are flagged.</P>"""%(stats[0]*100,stats[1]*100));
    
# The transfer_legacy_flags job uses Calico.Flagger to transfer legacy flags into a bitflag
stat_msgbox = None;
def transfer_legacy_flags (mqs,parent,**kw):
  if not flagger:
    raise RuntimeError,"Flagger not available, perhaps MS was not specified?";
  replace = (transfer_flag_policy == "replace");
  # open progress meter if GUI available
  init_progress_dialog(parent,"Transferring legacy flags");
  # execute flagging
  try:
    stats = flagger.transfer(flag=transfer_flag_selector.get_flagmask(),replace=replace,
                     progress_callback=progress_callback);
    flagger.close();
  finally:
    progress_callback(100,100);
    progress_dialog.hide();
  GUI.info_box("Flag statistics",
    """<P>%.4f%% of rows are flagged.</P>
      <P>%.4f%% of individual correlations are flagged.</P>"""%(stats[0]*100,stats[1]*100));
  
# The get_flag_stats job uses Calico.Flagger to get flag statistics and display them
def get_flag_stats (mqs,parent,**kw):
  if not flagger:
    raise RuntimeError,"Flagger not available, perhaps MS was not specified?";
  # open progress meter if GUI available
  init_progress_dialog(parent,"Getting flag statistics");
  # execute flagging
  try:
    stats = flagger.get_stats(flag=stat_flag_selector.get_flagmask(),
                              legacy=stat_flag_selector.read_legacy_flags,
                              progress_callback=progress_callback);
    flagger.close();
  finally:
    progress_callback(100,100);
    progress_dialog.hide();
  GUI.info_box("Flag statistics",
      """<P>%.4f%% of rows are flagged.</P>
      <P>%.4f%% of individual correlations are flagged.</P>"""%(stats[0]*100,stats[1]*100));

# The fill_legacy_flags job uses Calico.Flagger to fill legacy flags from bitflags
def fill_legacy_flags (mqs,parent,**kw):
  if not flagger:
    raise RuntimeError,"Flagger not available, perhaps MS was not specified?";
  # open progress meter if GUI available
  init_progress_dialog(parent,"Filling legacy flag column");
  # execute flagging
  try:
    flagger.set_legacy_flags(flags=fill_flag_selector.get_flagmask(),
                            progress_callback=progress_callback);
    flagger.close();
  finally:
    progress_callback(100,100);
    progress_dialog.hide();

def clear_flagset (mqs,parent,**kw):
  if not flagger:
    raise RuntimeError,"Flagger not available, perhaps MS was not specified?";
  fsets = remove_flag_selector.selected_flagsets();
  if GUI.warning_box("Clearing flagsets","""<P>This will clear all 
	flags in the selected flagsets (%s). Click OK to proceed.</P>"""%" ".join(fsets),
	GUI.Button.Ok|GUI.Button.Cancel,GUI.Button.Cancel) != GUI.Button.Ok:
    return;
  # open progress meter if GUI available
  init_progress_dialog(parent,"Clearing flags");
  # execute flagging
  try:
    flagger.unflag(remove_flag_selector.get_flagmask(),
                   progress_callback=progress_callback);
    flagger.close();
  finally:
    progress_callback(100,100);
    progress_dialog.hide();

def clear_bitflags (mqs,parent,**kw):
  if not flagger:
    raise RuntimeError,"Flagger not available, perhaps MS was not specified?";
  fsets = remove_flag_selector.selected_flagsets();
  if GUI.warning_box("Clearing bitflags","""<P>This will clear all bitflags in
	all flagsets. Click OK to proceed.</P>""",
	GUI.Button.Ok|GUI.Button.Cancel,GUI.Button.Cancel) != GUI.Button.Ok:
    return;
  # open progress meter if GUI available
  init_progress_dialog(parent,"Clearing bitflags");
  # execute flagging
  try:
    flagger.unflag(-1,
                   progress_callback=progress_callback);
    flagger.close();
  finally:
    progress_callback(100,100);
    progress_dialog.hide();

def clear_legacy_flags (mqs,parent,**kw):
  if not flagger:
    raise RuntimeError,"Flagger not available, perhaps MS was not specified?";
  if GUI.warning_box("Clearing FLAG/FLAG_ROW","""<P>This will clear all 
	flags from the FLAG/FLAG_ROW columns. Click OK to proceed.</P>""",
	GUI.Button.Ok|GUI.Button.Cancel,GUI.Button.Cancel) != GUI.Button.Ok:
    return;
  # open progress meter if GUI available
  init_progress_dialog(parent,"Clearing flags");
  # execute flagging
  try:
    flagger.clear_legacy_flags(progress_callback=progress_callback);
    flagger.close();
  finally:
    progress_callback(100,100);


def remove_flagset (mqs,parent,**kw):
  if not flagger:
    raise RuntimeError,"Flagger not available, perhaps MS was not specified?";
  fsets = remove_flag_selector.selected_flagsets();
  if GUI.warning_box("Removing flagsets","""<P>This will completely 
           remove the selected flagsets (%s). Click OK to proceed.</P>"""%" ".join(fsets),
	GUI.Button.Ok|GUI.Button.Cancel,GUI.Button.Cancel) != GUI.Button.Ok:
    return;
  # open progress meter if GUI available
  init_progress_dialog(parent,"Clearing flags");
  flagger.remove_flagset(*remove_flag_selector.selected_flagsets());
  # execute flagging
  try:
    flagger.unflag(remove_flag_selector.get_flagmask(),
                   progress_callback=progress_callback);
    flagger.close();
  finally:
    progress_callback(100,100);
    progress_dialog.hide();

view_ms_opt         = TDLRuntimeJob(view_ms_job,"View MS visibilities and flags",job_id="view_ms");
treeflag_ms_opt     = TDLRuntimeJob(treeflag_ms_job,"Flag MS using the in-tree flaggers",job_id="treeflag_ms");
add_bitflag_opt     = TDLRuntimeJob(add_bitflags,"Initialize bitflag columns in this MS",
  doc="""This MS does not contain any bitflag columns. Once you have initialized
  these columns, advanced flagging options will become available.""",job_id="add_bitflags");
run_autoflagger_opt = TDLRuntimeJob(run_autoflagger,"Run the autoflagger",job_id="autoflag");
flag_ms_opt         = TDLRuntimeJob(flag_ms,"Run the flagger",job_id="flag_ms");
transfer_opt        = TDLRuntimeJob(transfer_legacy_flags,"Transfer FLAG/FLAG_ROW column into this flagset",job_id="xfer_legacy")
fill_opt            = TDLRuntimeJob(fill_legacy_flags,"Fill FLAG/FLAG_ROW from these flagsets",job_id="fill_legacy")
get_stat_opt        = TDLRuntimeJob(get_flag_stats,"Get statistics",job_id="get_stats")
clear_bf_opt        = TDLRuntimeJob(clear_bitflags,"Clear all bitflags",
  doc="""Clears all bitflag columns completely. Use this option if your BITFLAG columns are in error somehow""",job_id="clear_bitflags");
clear_fs_opt        = TDLRuntimeJob(clear_flagset,"Clear selected flagset(s)",
  doc="""Clears all flags in the selected flagset(s), but does not remove the flagsets.""",job_id="clear_flagsets");
clear_legacy_opt    = TDLRuntimeJob(clear_legacy_flags,"Clear flags from FLAG/FLAG_ROW columns",job_id="clear_legacy");
remove_fs_opt       = TDLRuntimeJob(remove_flagset,"Remove selected flagset(s)",
  doc="""Completely removes the selected flagsets.""",job_id="remove_flagsets");

ms_job_options = [ view_ms_opt,treeflag_ms_opt,run_autoflagger_opt,flag_ms_opt,transfer_opt,clear_legacy_opt,
                   get_stat_opt,clear_bf_opt,clear_fs_opt,remove_fs_opt,add_bitflag_opt ];

TDLRuntimeMenu("View or flag MS visibilities",
  *( mssel.runtime_options() +
     [ view_ms_opt,treeflag_ms_opt ]));
treeflag_ms_opt.show(False);          # will be shown by define_forest() if needed


flag_flag_selector = mssel.make_write_flag_selector(namespace='flag_fl');
flag_menu = TDLRuntimeMenu("Flag or unflag a subset of the data",
  *([ TDLOption('flag_antennas',"Antenna(s)",[None],more=str),
      TDLOption('flag_baselines',"Baseline(s)",[None],more=str),
      TDLMenu("Time range",
        TDLOption('flag_time0',"Starting time (secs relative to first timeslot)",[None],more=float),
        TDLOption('flag_time1',"Ending time (secs relative to first timeslot)",[None],more=float),
        toggle='flag_reltime')
    ]
    + flag_subset.option_list() 
    + flag_flag_selector.option_list() 
    + [ TDLOption('flag_action',"Flag or unflag",["flag","unflag"]),
        flag_ms_opt 
    ]
  )
);

remove_flag_selector = mssel.make_read_flag_selector(namespace='rm_fl',
    legacy=False,label="Flagset %s",doc="");
remove_menu = TDLRuntimeMenu("Remove or clear flagset(s)",
  *(remove_flag_selector.option_list() +
    [clear_fs_opt,remove_fs_opt,clear_bf_opt])
  );

fill_flag_selector = mssel.make_read_flag_selector(namespace='fill_fl',legacy=False);
change_menu = TDLRuntimeMenu("Change the FLAG/FLAG_ROW columns",
  *(fill_flag_selector.option_list() +
    [fill_opt,clear_legacy_opt])
  );
  
autoflag_menu = TDLRuntimeMenu("Run autoflagger (fills FLAG/FLAG_ROW column)",
  TDLOption('autoflag_column',"MS column",["DATA","MODEL_DATA","CORRECTED_DATA"]),
  TDLMenu("Use sliding time median method (\"timemed\")",toggle='autoflag_timemed',
     *[ TDLOption('autoflag_timemed_thr',"Flagging threshold (sigmas)",5.,more=float),
        TDLOption('autoflag_timemed_hw' ,"Time window half-width (timeslots)",10,more=int),
        TDLMenu("Enable row flagging",
          TDLOption('autoflag_timemed_rowthr',"Row flagging threshold (sigmas)",10.,more=float),
          TDLOption('autoflag_timemed_rowhw' ,"Row window half-width (timeslots)",6,more=int),
          toggle='autoflag_timemed_flagrow'),
        TDLOption('autoflag_timemed_expr' ,"Flagging expression",["ABS I","ABS Q"],more=str),
        TDLOption('autoflag_timemed_fignore' ,"Ignore existing flags",False)
     ]),
  TDLMenu("Use sliding freq median method (\"freqmed\")",toggle='autoflag_freqmed',
     *[ TDLOption('autoflag_freqmed_thr',"Flagging threshold (sigmas)",5.,more=float),
        TDLOption('autoflag_freqmed_hw' ,"Freq window half-width (channels)",10,more=int),
        TDLMenu("Enable row flagging",
          TDLOption('autoflag_freqmed_rowthr',"Row flagging threshold (sigmas)",10.,more=float),
          TDLOption('autoflag_freqmed_rowhw' ,"Row window half-width (timeslots)",6,more=int),
          toggle='autoflag_freqmed_flagrow'),
        TDLOption('autoflag_freqmed_expr' ,"Flagging expression",["ABS I","ABS Q"],more=str),
        TDLOption('autoflag_freqmed_fignore' ,"Ignore existing flags",False)
     ]),
  TDLMenu("Use global time median method (\"newtimemed\")",toggle='autoflag_newtimemed',
     *[ TDLOption('autoflag_newtimemed_thr',"Flagging threshold (sigmas)",3.,more=float),
        TDLOption('autoflag_newtimemed_expr' ,"Flagging expression",["ABS I","ABS Q"],more=str),
        TDLOption('autoflag_newtimemed_fignore' ,"Ignore existing flags",False)
     ]),
  TDLMenu("Use spectral rejection method (\"sprej\")",toggle='autoflag_sprej',
     *[ TDLOption('autoflag_sprej_rowthr',"Flagging threshold (sigmas)",3.,more=float),
        TDLOption('autoflag_sprej_rowhw' ,"Row window half-width (timeslots)",6,more=int),
        TDLOption('autoflag_sprej_ndeg' ,"Degree of polynomial fit",2,more=int),
        TDLOption('autoflag_sprej_spwid' ,"Spectral window ID",["all"],more=int),
        TDLOption('autoflag_sprej_chan0' ,"Starting channel for fit",[0],more=int),
        TDLOption('autoflag_sprej_chan1' ,"Ending channel for fit",[1],more=int),
        TDLOption('autoflag_sprej_expr' ,"Flagging expression",["ABS I","ABS Q"],more=str),
        TDLOption('autoflag_sprej_fignore' ,"Ignore existing flags",False)
     ]),
  TDLMenu("Use UV-binner method (\"uvbin\")",toggle='autoflag_uvbin',
     *[ TDLOption('autoflag_uvbin_thr',"Probability cutoff (i.e. ratio of points to be flagged",[0,0.01],default=0.01,more=float),
        TDLOption('autoflag_uvbin_minpop' ,"Bin count cutoff",[0,1,10],more=int),
        TDLOption('autoflag_uvbin_nbin_uv' ,"Number of bins along UV axis",[50,100],more=int),
        TDLOption('autoflag_uvbin_nbin_val' ,"Number of bins along value axis",[50,100],more=int),
        TDLOption('autoflag_uvbin_plotchan' ,"Produce flag report plot based on channel",[None,"middle",0],more=int,default="middle"),
        TDLOption('autoflag_uvbin_econoplot',"Make abbreviated plot",True),
        TDLOption('autoflag_uvbin_expr' ,"Flagging expression",["ABS I","ABS Q"],more=str),
        TDLOption('autoflag_uvbin_fignore',"Ignore existing flags",False)
     ]),
## PGPlot windows don't play well with the browser, so disabling this for now
#  TDLMenu("Produce on-screen flag report",toggle='autoflag_plotscr',
#     *[ TDLOption('autoflag_plotscr_nx',"Plots per window (horizontal)",[1,2,3,4],more=int,default=2),
#        TDLOption('autoflag_plotscr_ny',"Plots per window (vertical)",[1,2,3,4],more=int,default=2),
#     ]),
  TDLMenu("Restrict to subset of data",toggle='autoflag_setdata',
     *[ TDLOption('autoflag_setdata_fieldid',"Field ID",["all"],more=int),
        TDLOption('autoflag_setdata_spwid',"Spectral window ID",[0],more=int),
        TDLOption('autoflag_setdata_chan' ,"Channels (\"all\" or start,end)",["all"],more=str),
        TDLOption('autoflag_setdata_msselect' ,"Additional TaQL selection string",[None],more=str),
     ]),
  TDLMenu("Produce hardcopy flag report",toggle='autoflag_plotdev',
     *[ TDLOption('autoflag_plotdev_nx',"Plots per window (horizontal)",[1,2,3,4],more=int,default=2),
        TDLOption('autoflag_plotdev_ny',"Plots per window (vertical)",[1,2,3,4],more=int,default=2),
        TDLOption('autoflag_plotdev_file',"Output filename (or PGPlot device)",["flagreport.ps"],more=str)
     ]),
  TDLOption('autoflag_reset',"Reset existing flags",False),
  TDLOption('autoflag_trial',"Trial run only, do not write flags",False),
  TDLOption('autoflag_save',"Save flagging commands to file",[None,'default.af'],more=str),
  run_autoflagger_opt,
  );

transfer_flag_selector = mssel.make_write_flag_selector(namespace='xfer_fl');
transfer_menu = TDLRuntimeMenu("Transfer FLAG/FLAG_ROW column into a flagset",
  *(transfer_flag_selector.option_list() +
    [ TDLOption('transfer_flag_policy',"Replace or add to flagset",["replace","add to"]),
      transfer_opt ]
  )
);

stat_flag_selector = mssel.make_read_flag_selector(namespace='stat_fl',legacy=True);
stat_menu = TDLRuntimeMenu("Get flag statistics",
  *(stat_flag_selector.option_list() +
    [ get_stat_opt ]
  )
);
# insert add bitflags job last
# TDLRuntimeOptions(add_bitflag_opt);

# _open_ms will be called whenever a new MS is selected
flagger = None;
def _open_ms (msname):
  global flagger;
  # init a flagger
  flagger = Calico.Flagger.Flagger(msname,verbose=5,chunksize=50000);
  # show/hide bitflag-related options
  add_bitflag_opt.show(not flagger.has_bitflags);
  for opt in flag_menu,remove_menu,fill_opt,transfer_menu:
    opt.show(flagger.has_bitflags);
  # close flagger
  flagger.close();
mssel.when_changed(_open_ms);


import Purr.Pipe

def _define_forest(ns,parent=None,**kw):
  if run_purr:
    Timba.TDL.GUI.log_message("starting purr");
    Timba.TDL.GUI.purr(mssel.msname+".purrlog",[mssel.msname,'.']);

  mssel.setup_observation_context(ns);
  array = Meow.Context.array;

  # make spigot nodes
  outputs = spigots = array.spigots(corr=mssel.get_corr_index());
  
  # ...and an inspector for them
  Meow.StdTrees.vis_inspector(ns.inspector('input'),spigots,
                              bookmark="Inspect input visibilities");
  inspectors = [ ns.inspector('input') ];
  Bookmarks.make_node_folder("Input visibilities by baseline",
    [ spigots(p,q) for p,q in array.ifrs() ],sorted=True,ncol=2,nrow=2);

  # add in-tree flaggers, if appropriate
  if tree_flag_enable:
    flaggers = [];
    # make nodes to take complex amplitude, if any flaggers are going to use them
    if tree_flag_ampl is not None or tree_flag_mean_ampl is not None or tree_flag_mean_freq_ampl is not None:
      # use matrix norm
      if tree_flag_norm:
        for p,q in array.ifrs():
          node = outputs(p,q);
          nsq = node('sq') << Meq.MatrixMultiply(node,node('conj') << Meq.ConjTranspose(node));
          nsqtr = nsq('tr') << (nsq(11) << Meq.Selector(nsq,index=0)) + (nsq(22) << Meq.Selector(nsq,index=3));
          ns.ampl(p,q) << Meq.Sqrt(nsqtr('abs') << Meq.Abs(nsqtr));
      # use per-element amplitudes
      else:
        for p,q in array.ifrs():
          ns.ampl(p,q) << Meq.Abs(outputs(p,q));
      # ...and an inspector for them
      Meow.StdTrees.vis_inspector(ns.inspector('ampl'),ns.ampl,bookmark="Inspect mean-freq visibility amplitudes");
      inspectors.append(ns.inspector('ampl'));
    # make flagger for amplitudes
    if tree_flag_ampl is not None:
      for p,q in array.ifrs():
        ns.flagampl(p,q) << Meq.ZeroFlagger(ns.ampl(p,q)-tree_flag_ampl,oper='gt',flag_bit=Meow.MSUtils.FLAGMASK_OUTPUT);
      flaggers.append(ns.flagampl);
    # make flagger for mean-ifr amplitudes
    if tree_flag_mean_ampl is not None:
      ns.meanampl << Meq.Mean(*[ns.ampl(p,q) for p,q in array.ifrs()]);
      ns.flagmeanampl << Meq.ZeroFlagger(ns.meanampl-tree_flag_mean_ampl,oper='gt',flag_bit=Meow.MSUtils.FLAGMASK_OUTPUT);
      Meow.Bookmarks.Page("Mean-ifr visibility amplitudes").add(ns.meanampl,viewer="Result Plotter");
      flaggers.append(lambda p,q:ns.flagmeanampl);
    # make flagger for mean-fr-freq amplitudes
    if tree_flag_mean_freq_ampl is not None:
      ns.meanampl ** Meq.Mean(*[ns.ampl(p,q) for p,q in array.ifrs()]);
      ns.meanfreqampl << Meq.Mean(ns.meanampl,reduction_axes="freq");
      ns.flagmeanfreqampl << Meq.ZeroFlagger(ns.meanfreqampl-tree_flag_mean_freq_ampl,oper='gt',flag_bit=Meow.MSUtils.FLAGMASK_OUTPUT);
      Meow.Bookmarks.Page("Mean-ifr-freq visibility amplitudes").add(ns.meanfreqampl,viewer="Collections Plotter");
      flaggers.append(lambda p,q:ns.flagmeanfreqampl);
    
    # now merge all flaggers  
    if flaggers:
      for p,q in array.ifrs():
        ns.flagged(p,q) << Meq.MergeFlags(outputs(p,q),*[fl(p,q) for fl in flaggers]);
      outputs = ns.flagged;
      Meow.StdTrees.vis_inspector(ns.inspector('flagged'),outputs,
                                  bookmark="Inspect flagged visibilities");
      inspectors.append(ns.inspector('flagged'));
      Bookmarks.make_node_folder("Flagged visibilities by baseline",
        [ outputs(p,q) for p,q in array.ifrs() ],sorted=True,ncol=2,nrow=2);
      mssel.enable_write_flags(True);
      treeflag_ms_opt.show(True);

  if outputs is spigots:
    ns.VisDataMux << Meq.VisDataMux(post=(inspectors[0] if len(inspectors)<2 else Meq.ReqMux(*inspectors)));
  else:
    Meow.StdTrees.make_sinks(ns,outputs,spigots=spigots,post=inspectors);
  
  # add imaging options
  imsel = mssel.imaging_selector(npix=512,arcmin=120);
  TDLRuntimeMenu("Make an image from this MS",*imsel.option_list());



