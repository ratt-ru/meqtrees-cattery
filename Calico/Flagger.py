# -*- coding: utf-8 -*-
import numpy
import Timba.dmi
import re
import tempfile
import os

from Meow import MSUtils
import Purr.Pipe

# figure out which table implementation to use -- try pyrap/casacore first
try:
  import pyrap_tables
  TABLE = pyrap_tables.table
  print "Calico flagger: using the pyrap_tables module"
except:
  # else try the old pycasatable/aips++ thing
  try:
    import pycasatable
    TABLE = pycasatable.table
    print "Calico flagger: using the pycasatable module. WARNING: this is deprecated."
    print "Please install pyrap and casacore!"
  except:
    TABLE = None;
    print "Calico flagger: no tables module found"
    print "Please install pyrap and casacore!"
    
_gli = MSUtils.find_exec('glish');
if _gli:
  _GLISH = 'glish';
  print "Calico flagger: found %s, autoflag should be available"%_gli;
else:
  _GLISH = None;
  print "Calico flagger: glish not found, autoflag will not be available";

_addbitflagcol = MSUtils.find_exec('addbitflagcol');

# Various argument-formatting methods to use with the Flagger.AutoFlagger class
# below. These really should be static methods of the class, but that doesn't work
# with Python (specifically, I cannot include them into static member dicts)
def _format_nbins (bins,argname):
  if isinstance(bins,(list,tuple)) and len(bins) == 2:
    return str(list(bins));
  else:
    return "%d"%bins;
  raise TypeError,"invalid value for '%s' keyword (%s)"%(argname,bins);
  
def _format_bool (arg,argname):
  if arg: return "T";
  else:   return "F";

def _format_index (arg,argname):
  if isinstance(arg,int):
    if arg<0:
      return str(arg);
    else:
      return str(arg+1);
  raise TypeError,"invalid value for '%s' keyword (%s)"%(argname,arg);

def _format_list (arg,argname):
  if isinstance(arg,str):
    return "'%s'"%arg;
  elif isinstance(arg,(list,tuple)):
    return "[%s]"%(','.join([_format_list(x,argname) for x in arg]));
  elif isinstance(arg,bool):
    return _format_bool(arg,argname);
  elif isinstance(arg,(int,float)):
    return str(arg);
  raise TypeError,"invalid value for '%s' keyword (%s)"%(argname,arg);

def _format_ilist (arg,argname):
  if isinstance(arg,str):
    return "'%s'"%arg;
  elif isinstance(arg,(list,tuple)):
    return "[%s]"%(','.join([_format_ilist(x,argname) for x in arg]));
  elif isinstance(arg,bool):
    return _format_bool(arg,argname);
  elif isinstance(arg,int):
    return _format_index(arg,argname);
  elif isinstance(arg,float):
    return str(arg);
  raise TypeError,"invalid value for '%s' keyword (%s)"%(argname,arg);

def _format_plotchan (arg,argname):
  if isinstance(arg,bool):
    return _format_bool(arg,argname);
  else:
    return _format_index(arg,argname);
  
def _format_2N (arg,argname):
  if isinstance(arg,(list,tuple)):
    return "%s"%arg;
  elif Timba.dmi.is_array(arg):
    if arg.dtype == Timba.array.int32:
      a = arg + 1;
    else:
      a = arg;
    if a.ndim == 1:
      return str(list(a));
    elif a.ndim == 2 and a.shape[1] == 2:
      return "[%s]"%(','.join([ str(list(a[i,:])) for i in range(a.shape[0]) ]));
  raise TypeError,"invalid value for '%s' keyword (%s)"%(argname,arg);

def _format_clip (arg,argname):
  if isinstance(arg,(list,tuple)):
    if not isinstance(arg[0],(list,tuple)):
      arg = [ arg ];
    recfields = [];
    for i,triplet in enumerate(arg):
      if not isinstance(triplet,(list,tuple)) or len(triplet) != 3:
        raise TypeError,"invalid value for '%s' keyword (%s)"%(argname,arg);
      expr,minval,maxval = triplet;
      subfields = [ "expr='%s'"%expr ];
      if minval is not None:
        subfields.append("min=%g"%minval);
      if maxval is not None:
        subfields.append("max=%g"%maxval);
      recfields.append("a%d=[%s]"%(i,','.join(subfields)));
    return "[%s]"%','.join(recfields);
  raise TypeError,"invalid value for '%s' keyword (%s)"%(argname,arg);
    
class Flagger (Timba.dmi.verbosity):
  def __init__ (self,msname,verbose=0,chunksize=200000):
    Timba.dmi.verbosity.__init__(self,name="Flagger");
    self.set_verbose(verbose);
    if not TABLE:
      raise RuntimeError,"No tables module found. Please install pyrap and casacore";
    self.msname = msname;
    self.ms = None;
    self.chunksize = chunksize;
    self._reopen();
    
  def close (self):
    if self.ms:
      self.dprint(2,"closing MS",self.msname);
      self.ms.close();
    self.ms = self.flagsets = None;
    
  def _reopen (self):
    if self.ms is None:
      self.ms = ms = TABLE(self.msname,readonly=False);
      self.dprintf(1,"opened MS %s, %d rows\n",ms.name(),ms.nrows());
      self.has_bitflags = 'BITFLAG' in ms.colnames();
      self.flagsets = MSUtils.get_flagsets(ms);
      self.dprintf(1,"flagsets are %s\n",self.flagsets.names());
      self.purrpipe = Purr.Pipe.Pipe(self.msname);
    return self.ms;
    
  def add_bitflags (self,wait=True,purr=True):
    if not self.has_bitflags:
      global _addbitflagcol; 
      if not _addbitflagcol:
        raise RuntimeError,"cannot find addbitflagcol utility in $PATH";
      self.close();
      self.dprintf(1,"running addbitflagcol\n");
      if os.spawnvp(os.P_WAIT,_addbitflagcol,['addbitflagcol',self.msname]):
        raise RuntimeError,"addbitflagcol failed";
      # report to purr
      if purr:
        self.purrpipe.title("Flagging").comment("Adding bitflag columns.");
      # reopen and close to pick up new BITFLAG column
      self._reopen();
      self.close();
      
  def remove_flagset (self,*fsnames,**kw):
    self._reopen();
    mask = self.flagsets.remove_flagset(*fsnames);
    # report to purr
    if kw.get('purr',True):
      self.purrpipe.title("Flagging").comment("Removing flagset(s) %s."%(','.join(fsnames)));
    self.unflag(mask,purr=False);

  def flag (self,flag=1,**kw):
    kw = kw.copy();
    if kw.setdefault('purr',True):
      if isinstance(flag,int):
        fset = "bitflag %x"%flag;
      else:
        fset = "flagset %s"%flag;
      self.purrpipe.title("Flagging").comment("Flagging %s"%fset,endline=False);
    kw['flag'] = flag;
    kw['unflag'] = 0;
    return self._flag(**kw);
  
  def unflag (self,unflag=-1,**kw):
    kw = kw.copy();
    if kw.setdefault('purr',True):
      if isinstance(unflag,int):
        fset = "bitflag %x"%unflag;
      else:
        fset = "flagset %s"%unflag;
      self.purrpipe.title("Flagging").comment("Unflagging %s"%fset,endline=False);
    kw['flag'] = 0;
    kw['unflag'] = unflag;
    return self._flag(**kw);
  
  def transfer (self,flag=1,replace=False,*args,**kw):
    unflag = (replace and flag) or 0;
    kw.setdefault('purr',True);
    if kw['purr']:
      if isinstance(flag,int):
        fset = "bitflag %x"%flag;
      else:
        fset = "flagset %s"%flag;
      self.purrpipe.title("Flagging").comment("Transferring FLAG/FLAG_ROW to %s"%fset,endline=False);
    return self._flag(flag=flag,unflag=unflag,transfer=True,*args,**kw);
  
  def get_stats(self,flag=0,legacy=False,**kw):
    kw.setdefault('purr',True);
    if kw['purr']:
      if isinstance(flag,int) and flag:
        fset = "bitflag %x"%flag;
      else:
        fset = "flagset %s"%flag;
      if legacy:
        fset += ", plus FLAG/FLAG_ROW";
      self.purrpipe.title("Flagging").comment("Getting flag stats for %s"%fset,endline=False);
    stats = self._flag(flag=flag,get_stats=True,include_legacy_stats=legacy,**kw);
    msg = "%.2f%% of rows and %.2f%% of correlations are flagged."%(stats[0]*100,stats[1]*100);
    if kw['purr']:
      self.purrpipe.comment(msg);
    self.dprint(1,"stats: ",msg);
    return stats;
  
  def _get_bitflag_col (self,ms,row0,nrows,shape=None):
    """helper method. Gets the bitflag column at the specified location. On error (presumably,
    column is missing), returns zero array of specified shape. If shape is not specified,
    queries the DATA column to obtain it.""";
    try:
      return ms.getcol('BITFLAG',row0,nrows);
    except:
      if not shape:
        shape = ms.getcol('DATA',row0,nrows).shape;
      return numpy.zeros(shape,dtype=numpy.int32);
  
  def _get_submss (self,ms):
    """Helper method. Splits MS into subsets by DATA_DESC_ID. 
    Returns list of (ddid,nrows,subms) tuples, where subms is a subset of the MS with the given DDID,
    and nrows is the number of rows in preceding submss (which is needed for progress stats)
    """;
    nddids = TABLE(os.path.join(ms.name(),'DATA_DESCRIPTION'),readonly=True).nrows();
    # build up list of (ddid,nrows,subms) pairs.
    # subms is a subset of the MS with the given DDID
    # nrows is the number of rows in preceding submss (need for progress stats)
    sub_mss = [];
    nrows = 0;
    for ddid in range(nddids):
      subms = ms.query("DATA_DESC_ID==%d"%ddid);
      sub_mss.append((ddid,nrows,subms));
      nrows += subms.nrows();
    return sub_mss;

  def _flag (self,flag=1,unflag=0,create=False,transfer=False,
          get_stats=False,include_legacy_stats=False,
          ddid=None,fieldid=None,
          channels=None,multichan=None,corrs=None,
          antennas=None,
          baselines=None,
          time=None,reltime=None,
          taql=None,
          progress_callback=None,
          purr=False
          ):
    """Internal _flag method does the actual work.""";
    ms = self._reopen();
    if not self.has_bitflags:
      if transfer:
        raise TypeError,"MS does not contain a BITFLAG column, cannot use flagsets""";
      if get_stats and flag:
        raise TypeError,"MS does not contain a BITFLAG column, cannot get statistics""";
    if get_stats and not flag and not include_legacy_stats:
      flag = -1;
    # lookup flagset name, if so specified
    if isinstance(flag,str):
      flagname = flag;
      flag = self.flagsets.flagmask(flag,create=create);
      self.dprintf(2,"flagset %s corresponds to bitmask 0x%x\n",flagname,flag);
    # lookup same for unflag, except we don't create
    if isinstance(unflag,str):
      flagname = unflag;
      unflag = self.flagsets.flagmask(unflag);
      self.dprintf(2,"flagset %s corresponds to bitmask 0x%x\n",flagname,unflag);
    if self.flagsets.names() is not None:
      if flag:
        self.dprintf(2,"flagging with bitmask 0x%x\n",flag);
      if unflag:
        self.dprintf(2,"unflagging with bitmask 0x%x\n",unflag);
    else:
      self.dprintf(2,"no bitflags in MS, using legacy FLAG/FLAG_ROW columns\n",unflag);
    # form up list of TaQL expressions for subset selectors
    single_ddid = None;
    queries = [];
    if taql:
      queries.append(taql);
    if ddid is not None:
      queries.append("DATA_DESC_ID==%d"%ddid);
      single_ddid = ddid;
    if fieldid is not None:
      queries.append("FIELD_ID==%d"%ddid);
    if antennas is not None:
      antlist = str(list(antennas));
      queries.append("(ANTENNA1 in %s || ANTENNA2 in %s)"%(antlist,antlist));
    if time is not None:
      t0,t1 = time;
      if t0 is not None:
        queries.append("TIME>=%g"%t0);
      if t1 is not None:
        queries.append("TIME<=%g"%t1);
    if reltime is not None:
      t0,t1 = reltime;
      time0 = self.ms.getcol('TIME',0,1)[0];
      if t0 is not None:
        queries.append("TIME>=%f"%(time0+t0));
      if t1 is not None:
        queries.append("TIME<=%f"%(time0+t1));
      
    # form up TaQL string, and extract subset of table
    if queries:
      query = ' && '.join(queries);
      purr and self.purrpipe.comment("; effective MS selection is \"%s\""%query,endline=False);
      self.dprintf(2,"selection string is %s\n",query);
      ms = ms.query(query);
      self.dprintf(2,"query reduces MS to %d rows\n",ms.nrows());
    else:
      self.dprintf(2,"no selection applied\n");
    # check list of baselines
    if baselines:
      baselines = [ (int(p),int(q)) for p,q in baselines ];
      purr and self.purrpipe.comment("; baseline subset is %s"%
        " ".join(["%d-%d"%(p,q) for p,q in baselines]),
        endline=False);
    
    # this will be true if only whole rows are being selected. If channel/correlation
    # criteria are supplied, this will be set to False below
    flagrows = True;
    # form up channel and correlation slicers
    # multichan may specify multiple channel subsets. If not specified,
    # then channels specifies a single subset. In any case, in the end multichan
    # will contain a list of the current channel selection
    if multichan is None:
      if channels is None:
        multichan = [ numpy.s_[:] ];
      else:
        flagrows = False;
        multichan = [ channels ];
        purr and self.purrpipe.comment("; channels are %s"%channels,endline=False);
    else:
      purr and self.purrpipe.comment("; channels are %s"%(', '.join(map(str,multichan))),endline=False);
      flagrows = False;
    self.dprintf(2,"channel selection is %s\n",multichan);
    if corrs is None:
      corrs = numpy.s_[:];
    else:
      purr and self.purrpipe.comment("; correlations %s"%corrs,endline=False);
      flagrows = False;
    print "purr is",purr;
    purr and self.purrpipe.comment(".");
    self.dprintf(2,"correlation selection is %s\n",corrs);
    stat_rows_nfl = stat_rows = stat_pixels = stat_pixels_nfl = 0;
    # make list of sub-MSs
    if single_ddid is not None:
      sub_mss = [ (single_ddid,0,ms) ];
    else:
      sub_mss = self._get_submss(ms);
    nrow_tot = ms.nrows();
    # go through rows of the MS in chunks
    for ddid,irow_prev,ms in sub_mss:
      self.dprintf(2,"processing MS subset for ddid %d\n",ddid);
      if progress_callback:
        progress_callback(irow_prev,nrow_tot);
      for row0 in range(0,ms.nrows(),self.chunksize):
        if progress_callback:
          progress_callback(irow_prev+row0,nrow_tot);
        nrows = min(self.chunksize,ms.nrows()-row0);
        self.dprintf(2,"flagging rows %d:%d\n",row0,row0+nrows-1);
        # get mask of matching baselines
        if baselines:
          # init mask of all-false
          rowmask = numpy.zeros(nrows,dtype=numpy.bool);
          a1 = ms.getcol('ANTENNA1',row0,nrows);
          a2 = ms.getcol('ANTENNA2',row0,nrows);
          # update mask
          for p,q in baselines:
            rowmask |= (a1==p) & (a2==q);
          self.dprintf(2,"baseline selection leaves %d rows\n",len(rowmask.nonzero()[0]));
        # else select all rows
        else:
          rowmask = numpy.s_[:];
          self.dprintf(2,"no baseline selection applied, flagging %d rows\n",nrows);
        # form up subsets for channel/correlation selector
        subsets = [ (rowmask,ch,corrs) for ch in multichan ];
        # first, handle statistics mode
        if get_stats:
          # collect row stats
          if include_legacy_stats:
            lfr  = ms.getcol('FLAG_ROW',row0,nrows)[rowmask];
            lf   = ms.getcol('FLAG',row0,nrows);
          else:
            lfr = lf = 0;
          if flag:
            bfr = ms.getcol('BITFLAG_ROW',row0,nrows)[rowmask];
            lfr = lfr + ((bfr&flag)!=0);
            bf = self._get_bitflag_col(ms,row0,nrows);
          # size seems to be a method or an attribute depending on numpy version :(
          stat_rows     += (callable(lfr.size) and lfr.size()) or lfr.size;
          stat_rows_nfl += lfr.sum();
          for subset in subsets:
            if include_legacy_stats:
              lfm = lf[subset];
            else:
              lfm = 0;
            if flag:
              lfm = lfm + (bf[subset]&flag)!=0;
            # size seems to be a method or an attribute depending on numpy version :(
            stat_pixels     += (callable(lfm.size) and lfm.size()) or lfm.size;
            stat_pixels_nfl += lfm.sum();
        # second, handle transfer-flags mode
        elif transfer:
          bf = ms.getcol('BITFLAG_ROW',row0,nrows);
          bfm = bf[rowmask];
          if unflag:
            bfm &= ~unflag;
          lf = ms.getcol('FLAG_ROW',row0,nrows)[rowmask];
          bf[rowmask] = numpy.where(lf,bfm|flag,bfm);
            # size seems to be a method or an attribute depending on numpy version :(
          stat_rows     += (callable(lf.size) and lf.size()) or lf.size;
          stat_rows_nfl += lf.sum();
          ms.putcol('BITFLAG_ROW',bf,row0,nrows);
          lf = ms.getcol('FLAG',row0,nrows);
          bf = self._get_bitflag_col(ms,row0,nrows,lf.shape);
          for subset in subsets:
            bfm = bf[subset];
            if unflag:
              bfm &= ~unflag;
            lfm = lf[subset]
            bf[subset] = numpy.where(lfm,bfm|flag,bfm);
            # size seems to be a method or an attribute depending on numpy version :(
            stat_pixels     += (callable(lfm.size) and lfm.size()) or lfm.size;
            stat_pixels_nfl += lfm.sum();
          ms.putcol('BITFLAG',bf,row0,nrows);
        # else, are we flagging whole rows?
        elif flagrows:
          if self.has_bitflags:
            bfr = ms.getcol('BITFLAG_ROW',row0,nrows);
            if unflag:
              bfr[rowmask] &= ~unflag;
              # if unflaging, also have to clear flags from BITFLAG too
              bf = self._get_bitflag_col(ms,row0,nrows);
              bf[rowmask,:,:] &= ~unflag;
              ms.putcol('BITFLAG',bf,row0,nrows);
            if flag:
              bfr[rowmask] |= flag;
            ms.putcol('BITFLAG_ROW',bfr,row0,nrows);
          else:
            lf = ms.getcol('FLAG_ROW',row0,nrows);
            lf[rowmask] = (flag!=0);
            ms.putcol('FLAG_ROW',lf,row0,nrows);
        # else flagging individual correlations or channels
        else: 
          # apply flagmask
          if self.has_bitflags:
            bf = self._get_bitflag_col(ms,row0,nrows);
            if unflag:
              for subset in subsets:
                bf[subset] &= ~unflag;
            if flag:
              for subset in subsets:
                bf[subset] |= flag;
            ms.putcol('BITFLAG',bf,row0,nrows);
          else:
            bf = ms.getcol('FLAG',row0,nrows);
            for subset in subsets:
              bf[subset] = (flag!=0);
            ms.putcol('FLAG',bf,row0,nrows);
    if progress_callback:
      progress_callback(nrow_tot,nrow_tot);
    stat0 = (stat_rows and stat_rows_nfl/float(stat_rows)) or 0;
    stat1 = (stat_pixels and stat_pixels_nfl/float(stat_pixels)) or 0;
    return stat0,stat1;
      
  def set_legacy_flags (self,flags,progress_callback=None,purr=True):
    """Fills the legacy FLAG/FLAG_ROW column by applying the specified flagmask
    to bitflags.
    * if flags is an int, it is used as a bitflag mask
    * if flags is a str, it is treated as the name of a flagset
    * if flags is a list or tuple, it is treated as a list of flagsets
    """;
    ms = self._reopen();
    if not self.has_bitflags:
      raise TypeError,"MS does not contain a BITFLAG column, cannot use bitflags""";
    if isinstance(flags,str):
      flagmask = self.flagsets.flagmask(flags);
      self.dprintf(1,"filling legacy FLAG/FLAG_ROW using flagset %s\n",flags);
      purr and self.purrpipe.title("Flagging").comment("Filling FLAG/FLAG_ROW from flagset %s."%flags);
    elif isinstance(flags,(list,tuple)):
      flagmask = 0;
      for fl in flags:
        flagmask |= self.flagsets.flagmask(fl);
      self.dprintf(1,"filling legacy FLAG/FLAG_ROW using flagsets %s\n",flags);
      purr and self.purrpipe.title("Flagging").comment("Filling FLAG/FLAG_ROW from flagsets %s."%','.join(flags));
    elif isinstance(flags,int):
      flagmask = flags;
      self.dprintf(1,"filling legacy FLAG/FLAG_ROW using flagsets %s\n",flags);
      purr and self.purrpipe.title("Flagging").comment("Filling FLAG/FLAG_ROW from bitflags %x."%flags);
    else:
      raise TypeError,"flagmask argument must be int, str or sequence";
    self.dprintf(1,"filling legacy FLAG/FLAG_ROW using bitmask 0x%d\n",flagmask);
    # now go through MS and fill the column
    # get list of per-DDID subsets
    sub_mss = self._get_submss(ms);
    nrow_tot = ms.nrows();
    # go through rows of the MS in chunks
    for ddid,irow_prev,ms in sub_mss:
      self.dprintf(2,"processing MS subset for ddid %d\n",ddid);
      if progress_callback:
        progress_callback(irow_prev,nrow_tot);
      for row0 in range(0,ms.nrows(),self.chunksize):
        if progress_callback:
          progress_callback(irow_prev+row0,nrow_tot);
        nrows = min(self.chunksize,ms.nrows()-row0);
        self.dprintf(2,"filling rows %d:%d\n",row0,row0+nrows-1);
        bf = self._get_bitflag_col(ms,row0,nrows,ms.getcol('FLAG').shape);
        bfr = ms.getcol('BITFLAG_ROW',row0,nrows);
        ms.putcol('FLAG',(bf&flagmask).astype(Timba.array.dtype('bool')),row0,nrows);
        ms.putcol('FLAG_ROW',(bfr&flagmask).astype(Timba.array.dtype('bool')),row0,nrows);
    if progress_callback:
      progress_callback(nrow_tot,nrow_tot);
      
  def clear_legacy_flags (self,progress_callback=None,purr=True):
    """Clears the legacy FLAG/FLAG_ROW columns.
    """;
    ms = self._reopen();
    self.dprintf(1,"clearing legacy FLAG/FLAG_ROW column\n");
    purr and self.purrpipe.title("Flagging").comment("Clearing FLAG/FLAG_ROW columns");
    # now go through MS and fill the column
    # go through rows of the MS in chunks
    shape = list(ms.getcol('FLAG',0,1).shape);
    shape[0] = self.chunksize
    fzero  = Timba.array.zeros(shape,dtype='bool');
    frzero = Timba.array.zeros((self.chunksize,),dtype='bool');
    # now go through MS and fill the column
    # get list of per-DDID subsets
    sub_mss = self._get_submss(ms);
    nrow_tot = ms.nrows();
    # go through each sub-MS, and through rows of the sub-MS in chunks
    for ddid,irow_prev,ms in sub_mss:
      self.dprintf(2,"processing MS subset for ddid %d\n",ddid);
      if progress_callback:
        progress_callback(irow_prev,nrow_tot);
      for row0 in range(0,ms.nrows(),self.chunksize):
        if progress_callback:
          progress_callback(row0+irow_prev,ms.nrows());
        nrows = min(self.chunksize,ms.nrows()-row0);
        self.dprintf(2,"filling rows %d:%d\n",row0,row0+nrows-1);
        fl = ms.getcol('FLAG',row0,nrows);
        fl[:,:,:] = False;
        ms.putcol('FLAG',fl,row0,nrows);
        ms.putcol('FLAG_ROW',Timba.array.zeros((nrows,),dtype='bool'),row0,nrows);
    if progress_callback:
      progress_callback(nrow_tot,nrow_tot);
  
  def autoflagger (self,*args,**kw):
    return Flagger.AutoFlagger(self,*args,**kw);
  
  class AutoFlagger (object):
    def __init__ (self,flagger,load=False):
      self.flagger = flagger;
      self._cmds = [];
      if load:
        if isinstance(load,bool):
          load = 'default.af';
        self.load(load);
      
    def _cmd (self,cmd):
      self._cmds.append(cmd);
      
    def reset (self):
      self._cmds = [];
      
    def setdata (self,chanstart=None,chanend=None,chanstep=None,spwid=None,fieldid=None,msselect=None):
      args = [];
      if chanstart is not None:
        chanstep = chanstep or 1;
        chanend  = chanend or chanstart;
        totchan = chanend-chanstart+1;
        nchan = totchan/chanstep;
        if totchan%chanstep:
          nchan += 1;
        args += [ "mode='channel'","nchan=%d"%nchan,"start=%d"%(start+1),"step=%d"%step ];
      if spwid is not None:
        args.append("spwid=%d"%(spwid+1));
      if fieldid is not None:
        args.append("fieldid=%d"%(fieldid+1));
      if msselect is not None:
        args.append("msselect='%s'"%msselect);
      self._cmd("af.setdata(%s);"%(','.join(args))); 
    
    _settimemed_dict    = dict(thr="%g",hw="%d",rowthr="%g",rowhw="%d",norow=_format_bool,
                                column="'%s'",expr="'%s'",fignore=_format_bool);
    _setnewtimemed_dict = dict(thr="%g",
                                column="'%s'",expr="'%s'",fignore=_format_bool);
    _setfreqmed_dict    = dict(thr="%g",hw="%d",rowthr="%g",rowhw="%d",norow=_format_bool,
                                column="'%s'",expr="'%s'",fignore=_format_bool);
    _setuvbin_dict      = dict(thr="%g",minpop="%d",nbins=_format_nbins,
                                plotchan=_format_plotchan,econoplot=_format_bool,
                                column="'%s'",expr="'%s'",fignore=_format_bool);
    _setsprej_dict      = dict(ndeg="%d",rowthr="%g",rowhw="%d",norow=_format_bool,
                                spwid=_format_ilist,fq=_format_2N,chan=_format_ilist,
                                column="'%s'",expr="'%s'",fignore=_format_bool);
    _setselect_dict     = dict(spwid=_format_ilist,field=_format_ilist,
                               fq=_format_2N,chan=_format_2N,corr=_format_list,
                               ant=_format_ilist,baseline=_format_ilist,timerng=_format_list,
                               autocorr=_format_bool,timeslot=_format_list,dtime="%g",
                               quack=_format_bool,unflag=_format_bool,clip=_format_clip);
    _setdata_dict       = dict(spwid=_format_ilist,field=_format_ilist,
                               nchan=_format_ilist,start=_format_ilist,
                               mode="'%s'",msselect="'%s'");
    _run_dict           = dict(plotscr=_format_list,plotdev=_format_list,devfile="'%s'",
                                reset=_format_bool,trial=_format_bool);
    
    def _setmethod (self,methodname,kwargs):
      argdict = getattr(self,'_%s_dict'%methodname);
      args = [];
      for kw,value in kwargs.iteritems():
        if value is not None:
          format = argdict.get(kw,None);
          if format is None:
            raise TypeError,"Autoflagger: invalid keyword '%s' passed to method %s()"%(kw,methodname);
          elif callable(format):
            args.append("%s=%s"%(kw,format(value,kw)));
          else:
            args.append("%s=%s"%(kw,format%value));
      return "af.%s(%s);"%(methodname,','.join(args)); 
      
    def settimemed (self,**kw):
      self._cmd(self._setmethod('settimemed',kw));
    def setfreqmed (self,**kw):
      self._cmd(self._setmethod('setfreqmed',kw));
    def setnewtimemed (self,**kw):
      self._cmd(self._setmethod('setnewtimemed',kw));
    def setsprej (self,**kw):
      self._cmd(self._setmethod('setsprej',kw));
    def setuvbin (self,**kw):
      self._cmd(self._setmethod('setuvbin',kw));
    def setselect (self,**kw):
      self._cmd(self._setmethod('setselect',kw));
    def setdata (self,**kw):
      if kw.get("nchan",None) is not None:
        kw["mode"] = "channel";
      # else if any other arguments are specified, set mode to "spwids", as this effectively
      # causes the flagger to select on everything except channels
      elif [ x for x in kw.itervalues() if x is not None ]:
        kw["mode"] = "spwids";
      self._cmd(self._setmethod('setdata',kw));
      
    def run (self,wait=True,cmdfile=None,purr=True,**kw):
      runcmd = self._setmethod('run',kw);
      # init list of command strings
      cmds = [ "include 'autoflag.g'","af:=autoflag('%s');"%self.flagger.msname ];
      # add default setdata() if not specified
      if not filter(lambda x:x.startswith('af.setdata'),self._cmds):
        cmds.append("af.setdata();");
      # add specified command set
      cmds += self._cmds;
      # add the run command
      cmds += [ runcmd ];
      self.flagger.dprint(2,"running autoflag with the following command set:")
      for cmd in cmds:
        self.flagger.dprint(2,cmd);
      if _GLISH is None:
        raise RuntimeError,"glish not found, so cannot run autoflagger";
      # write commands to temporary file and run glish
      if cmdfile:
        fobj = file(cmdfile,"wt");
      else:
        fh,cmdfile = tempfile.mkstemp(prefix="autoflag",suffix=".g");
        fobj = os.fdopen(fh,"wt");
#      cmds.append("shell('rm -f %s')"%cmdfile);
      cmds.append("exit");
      fobj.writelines([line+"\n" for line in cmds]);
      fobj.close();
      # write to pipe
      purr and self.flagger.purrpipe.title("Flagging").comment("Running autoflagger.").pounce(cmdfile);
      # tell flagger to detach from MS
      self.flagger.close();
      # spawn glish
      self.flagger.dprintf(3,"temp glish command file is %s\n"%cmdfile);
      waitcode = (wait and os.P_WAIT) or os.P_NOWAIT;
      return os.spawnvp(waitcode,_GLISH,['glish','-l',cmdfile]);
    
    def save (self,filename='default.af'):
      file(filename,"w").writelines([line+"\n" for line in self._cmds]);
      self.flagger.dprintf(2,"saved autoflag command sequence to file %s\n",filename);
      
    def load (self,filename='default.af'):
      self._cmds = file(filename).readlines();
      self.flagger.dprintf(2,"loaded autoflag command sequence from file %s\n",filename);
      self.flagger.dprint(2,"sequence is:");
      for cmd in self._cmds:
        self.flagger.dprint(3,cmd);
      