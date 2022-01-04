#!/usr/bin/env python3
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

import time
import os.path
import sys
import traceback
import pickle
import copy
import numpy
import numpy.ma

from Timba import dmi
from Timba import mequtils
from Timba.parmtables import FastParmTable
from Timba.Meq import meq

# debug printing
import Kittens.utils
verbosity = Kittens.utils.verbosity(name="ParmTables");
dprint = verbosity.dprint;
dprintf = verbosity.dprintf;

def verbose (level):
  """Changes the verbosity level""";
  verbosity.set_verbose(level);

def _int_or_str(x):
  """helper function, converts string to integer if possible""";
  try: return int(x);
  except: return x;

def cmp_qualified_names (a,b):
  """compares two qualified names: splits them at ':', and comparing the subfields in a numeric sense""";
  from past.builtins import cmp
  return cmp(list(map(_int_or_str,a.split(':'))),list(map(_int_or_str,b.split(':'))));

def sort_qualified_names (inlist):
  """Sorts qualified names as follows: splits into fields at ':', and sorts according to fields.
  If field is a numbers, sorts in the numeric sense."""
  from functools import cmp_to_key
  return sorted(inlist,key=cmp_to_key(cmp_qualified_names));

class _AxisStats (object):
  """_AxisStats represents information about one axis in the parmtable. It is created internally
  by ParmTab.""";
  def __init__ (self,name):
    self.name = str(name).lower();
    self.cells = {};

  def empty (self):
    return not self.cells;

  def add_cell (self,x1,x2):
    assert x1 is int and x2 is int
    x0 = (x1+x2)//2;
    self.cells[x0] = max(x2-x1,self.cells.get(x0,0));

  def update (self):
    if self.cells:
      # axis range
      self.minmax = min([x0-dx//2 for x0,dx in self.cells.items()]),max([x0+dx//2 for x0,dx in self.cells.items()]);
      # grid is simply a sorted list of cell values
      self.grid = list(self.cells.keys());
      self.grid.sort();
      # cell_index gives the orfinal number of each cell in the grid
      self.cell_index = dict([(x0,i) for i,x0 in enumerate(self.grid)]);

  def lookup_cell (self,x1,x2):
    return self.cell_index[(x1+x2)/2];

class DomainSlicing (list):
  """A DomainSlicing represents a slicing of the ParmTable domain.
  Basically, it is a list of slice_index tuples, with each tuple corresponding to a slice through
  the domain.
  """;
  def __init__ (self,slice_axes,domain_subset):
    """The constructor makes a slicing from a list of axes to be incorporated in the slice.
    All other axes will be iterated over.""";
    list.__init__(self);
    self.append([]);
    for iaxis,axis_subset in enumerate(domain_subset):
      if axis_subset is None or iaxis in slice_axes or mequtils.get_axis_id(iaxis) in slice_axes:
        self[:] = [ sl + [None] for sl in self ];
      else:
        self[:] = [ sl + [i] for sl in self for i in axis_subset ];

class FunkSlice (object):
  """FunkSlice represents a slice of funklets from a parmtable."""
  def __init__ (self,parmtab,name,funklist,index,iaxes,axes=None):
    self.pt = parmtab;
    self.name = name;
    self.funklets = funklist;
    self.slice_index = index;
    self.slice_iaxes = iaxes;
    self.slice_axes = axes or list(map(mequtils.get_axis_id,iaxes));
    self.rank = len(iaxes);
  def __len__ (self):
    return len(self.funklets);
  def __getitem__ (self,key):
    return self.funklets[key];
  def __iter__ (self):
    return iter(self.funklets);
  def array (self,coeff=0,fill_value=0,masked=True,collapse=True):
    """Returns funklet coefficients arranged into a hypercube.
    'coeff' is applied as an index into each funklet's coeff array, so coeff=0 or coeff=(0,0) selects
    c00, coeff=(0,1) selects c01, etc. IndexError is raised when a non-existing coeff is selected.
    The return value is an array with the same dimensions as the slice. If funklets for some
    domains are missing, the missing values will be set to fill_value. If masked=True, the return value
    will also be a masked array, with missing values masked. If masked=False, return value
    if an ordinary array.
    If 'collapse' is False, axes of the array will have a one-to-one correspondence with domain axes
        (so e.g. a slice for time=0,freq=0,l=*,m=* will have shape [1,1,nl,nm]
    If 'collapse' is True, axes not in this slice will be eliminated.
        (so e.g. the same slice will have shape [nl,nm]
    """
    # convert things like (0,0) into None
    if not coeff or not any(coeff):
      coeff = None;
    # initially array is of uncollapsed: one axis for each dimension (up to max_axies), with those not in the
    # slice having a size of 1. Extra axrs will be trimmed later.
    shape = [1]*mequtils.max_axis;
    for iaxis in self.slice_iaxes:
      shape[iaxis] = len(self.pt.axis_stats(iaxis).grid);
    # init empty arrays. Mask is all True initially, cleared as filled
    arr = numpy.zeros(shape,float);
    arr[...] = fill_value;
    mask = numpy.ones(shape,bool);
    # index object: modified and used below
    idx0 = [0]*mequtils.max_axis;
    # now go through funklets and fill them in
    for funk in self.funklets:
      if numpy.isscalar(funk.coeff):
        if coeff is not None:
          raise IndexError("invalid coeff index %s (funklet is scalar)"%coeff);
        val = funk.coeff;
      else:
        try:
          val = funk.coeff[coeff];
        except:
          raise IndexError("invalid coeff index %s (funklet coeffs are %s)"%(coeff,funk.coeff.shape));
      # funk.slice_index is the global index of this funklet. We want to use just the axes
      # found in our slice for assigning to the array
      for iaxis in self.slice_iaxes:
        idx0[iaxis] = funk.slice_index[iaxis];
      idx = tuple(idx0);
      try:
        arr[idx] = val;
        mask[idx] = False;
      except:
        raise;
    # make masked array if needed
    if masked:
      arr = numpy.ma.masked_array(arr,mask,fill_value=fill_value);
    # reshape
    if collapse:
      arr.shape = [ shape[iaxis] for iaxis in self.slice_iaxes ];
    else:
      arr.shape = shape[:(self.slice_iaxes[-1]+1)];
    return arr;

class FunkSet (object):
  """FunkSet represents a set of funklets with the same name""";
  def __init__ (self,parmtab,name):
    self.pt = parmtab;
    self.name = name;

  def get_slice (self,*index,**axes):
    """get_slice() returns a FunkSlice of all funklets matching a specified slice.
    Slice can be specified in one of two ways: by a vector of axis indices, or by keywords.
    E.g.:
      get_slice(None,0) or get_slice(freq=0)
    returns all funklets for 'name' and the first frequency (first form assumes freq is axis 1,
    so time=None and freq=0)
      get_slice(0) or get_slice(time=0)
    returns all funklets for 'name' and the first timeslot (first form assumes time is axis 0),
      get_slice(1,2) or get_slice(time=1,freq=2)
    returns all funklets for 'name', timeslot 1, frequency 2.
    """;
    # make sure index is a vector of max_axes length
    if len(index) < mequtils.max_axis:
      index = list(index) + [None]*(mequtils.max_axis - len(index));
    else:
      index = list(index);
    # set additional indices from keywords
    for axis,num in axes.items():
      index[mequtils.get_axis_number(axis)] = num;
    # build up list of full indices corresponding to specified slice
    slice_iaxis = [];
    indices = [[]];
    for iaxis,axis_idx in enumerate(index):
      stats = self.pt.axis_stats(iaxis);
      # for empty axes, or axes specifed in our call, append number to all indices as is
      if axis_idx is not None or stats.empty():
        list(map(lambda idx:idx.append(axis_idx),indices));
      # non-empty axes NOT specified in our call will be part of the slice
      else:
        slice_iaxis.append(iaxis);
        indices = [ idx + [i] for idx in indices for i in range(len(stats.grid)) ];
    # now make funklet list
    funkslice = FunkSlice(self.pt,self.name,[],index,slice_iaxis);
    for idx in indices:
      idx = tuple(idx);
      idom = self.pt._domain_reverse_index.get(idx,None);
      if idom is not None:
        funk = self.pt.parmtable().get_funklet(self.name,idom);
        if funk:
          funk.domain_index = idom;
          funk.slice_index = idx;
          funkslice.funklets.append(funk);
    return funkslice;

  def __call__ (self,*index,**axes):
    """The () operator on a FunkSet is equivalent to get_slice()""";
    return self.get_slice(*index,**axes);

  def array (self,coeff=0,fill_value=0,masked=True,collapse=True):
    """Makes array corresponding to whole FunkSet. This function will also maintain a disk cache
    of the array, and read it in or regenerate it as needed (unlike FunkSlice.array(), which
    always builds its arrays from scratch.)
    \n\n""" + FunkSlice.array.__doc__;
    # see if we have a cached array
    cachefile = os.path.join(self.pt.filename,"array.%s.%s.cache"%(self.name,coeff));
    arr = None;
    if os.path.exists(cachefile) and os.path.getmtime(cachefile) >= self.pt.mtime:
      try:
        arr = pickle.load(open(cachefile, "rb"));
        dprintf(2,"read cache %s\n"%cachefile);
      except:
        dprintf(0,"error reading cached array %s, will regenerate\n"%cachefile);
        arr = None;
    # regenerate array if not read
    if arr is None:
      dprintf(2,"filling array for %s.%s\n",self.name,coeff);
      fullslice = self.get_slice();
      # generate full, uncollapsed masked array
      arr = fullslice.array(coeff,fill_value=0,masked=True,collapse=False);
      # write to cache
      try:
        pickle.dump(arr,open(cachefile,"wb"));
      except:
        if verbosity.get_verbose() > 0:
          traceback.print_exc();
        dprintf(0,"error writing cache array %s, but proceeding anyway\n"%cachefile);
    # now apply the masked and fill_value properties
    if not masked:
      arr = arr.filled(fill_value);
    else:
      arr.fill_value = fill_value;
    # and collapse axes if asked
    if collapse:
      arr.shape = [ n for i,n in enumerate(arr.shape) if not self.pt.axis_stats(i).empty() ];
    return arr;

  def apply (self,op_func,slicing=[],outtab=None,remove=False):
    """For each funklet in the subset, takes all funklets along the designated slicing axis
    (i.e. for each slice along the non-listed axes), creates a FunkSlice, and calls
    op_func(slice).
    The return value of op_func() should be either None if no operation was performed, or a list of
    funklets to be written to the output table. This list may also contain strings, which are
    interpreted as funklet names. The default name is the same as the current FunkSet name; a string at
    any position in the funklet list applies to subsequent funklets.
    if 'outtab' is None, a new output table is created. Otherwise set outtab to a filename, or a
    ParmTab, or a FastParmTable.
    If 'remove' is True, input funklets will be removed if an output funklet is returned.
    """;
    # resolve tables
    outtab = self.pt.resolve_output_table(outtab);
    pt = self.pt.parmtable();
    # resolve slicing
    slicing = self.pt.make_slicing(slicing);
    # loop over funklets and slices
    num_infunk = num_outfunk = num_slices = 0;
    for sl0 in slicing:
      funklets = self.get_slice(*sl0);
      # call reduction function if we find any
      if funklets:
        num_slices += 1;
        num_infunk += len(funklets);
        dprintf(4,"%s slice %s: %d funklets found\n",self.name,sl0,len(funklets));
        dprint(6,"first funklet:",funklets[0]);
        dprint(6,"last funklet:",funklets[-1]);
        try:
          outfunk = op_func(funklets);
          dprint(6,"output funklets:",outfunk);
        except:
          if verbosity.get_verbose() > 0:
            dprintf(1,"exception computing funklets for %s slice %s\n",self.name,sl0);
            traceback.print_exc();
            dprintf(1,"this slice will be ignored\n");
          outfunk = None;
        # if a new funklet was generated, write it to output
        if outfunk:
          # remove input funklets if successful
          if remove:
            dprintf(4,"%s slice %s: removing %d input funklets\n",self.name,sl0,len(funklist));
            for funk in funklist:
              try:
                self.mtime = time.time();
                self.parmtable(True).delete_funklet(self.name,funk.domain_index);
              except:
                if verbosity.get_verbose() > 0:
                  traceback.print_exc();
                dprintf(0,"error deleting funklet for %s slice %s\n",self.name,funk.slice_index);
          dprintf(4,"%s slice %s: writing %d output funklets\n",self.name,sl0,len(outfunk));
          name = self.name;
          for ff in outfunk:
            if isinstance(ff,str):
              name = ff;
            else:
              num_outfunk += 1;
              try:
                outtab.mtime = time.time();
                outtab.parmtable(True).put_funklet(name,ff);
              except:
                dprintf(0,"error saving funklet for %s slice %s\n",self.name,sl0);
                if verbosity.get_verbose() > 0:
                  traceback.print_exc();
                dprintf(0,"this slice will be ignored\n");
                funk = None;
                break;
    dprintf(3,"%s: %s() transformed %d input funklets over %d slices into %d output funklets\n",self.name,op_func.__name__,num_infunk,num_slices,num_outfunk);



class ParmTab (object):
  """A ParmTab is a wrapper around a FastParmTable providing high-level facilities
   for dealing with funklets."""
  def __init__ (self,filename,write=False,new=False):
    """ParmTab constructor.
     'filename' is a ParmTable name, this will be created if it doesn't exist.
     'write' is True to open the table for writing (can always pass False, and re-open for writing later)
     'new' is True to delete any existing parmtable by the same name and open a new one.
    """;
    self._pt = self._pt_write = None;
    self.load(filename,write=write,new=new);

  def _start_progress(self,label,maxval=100):
    self._progress_label,self._progress_maxval = label,maxval;
    sys.stderr.write("%s: 00%%\r"%label);

  def _report_progress(self,value,maxval=None):
    maxval = maxval if maxval is not None else self._progress_maxval;
    sys.stderr.write("%s: %-40s\r"%(self._progress_label,"%02d%%"%(100*value//maxval)));

  def _end_progress(self,keep=True):
    if keep:
      self._report_progress(self._progress_maxval);
      sys.stderr.write("\n");
    else:
      sys.stderr.write("\r%80s\n"%"");

  def load (self,filename,write=False,new=False):
    """Loads the specified parmtable
     'write' is True to open the table for writing (can always pass False, and re-open for writing later)
     'new' is True to delete any existing parmtable by the same name and open a new one.
    """;
    if new:
      if os.path.exists(filename):
        dprintf(1,"removing old table %s\n",filename);
        os.system("/bin/rm -fr %s"%filename);
      dprintf(1,"creating new table %s\n",filename);
      write = True;
    else:
      dprintf(1,"loading table %s (write=%d)\n",filename,write);
    self.filename = filename;
    self.parmtable(write);
    self._cachepath = os.path.join(filename,'ParmTab.cache');
    self._make_axis_index();

  def merge (self,filename):
    """Merges in the specified parmtable""";
    pt1 = FastParmTable(filename);
    t0 = time.time();
    dprintf(2,"reading funklet list from %s\n",filename);
    funklist = pt1.funklet_list();
    dprintf(2,"elapsed time: %f seconds\n",time.time()-t0); t0 = time.time();
    nfunk = len(funklist);
    dprintf(1,"merging in %d funklets from table %s\n",nfunk,filename);
    if funklist:
      pt = self.parmtable(True);
      self._start_progress("merging in parmtable %s"%filename,nfunk);
      try:
        for ifunk,(name,idom,domain) in enumerate(funklist):
          self._report_progress(ifunk);
          self.mtime = time.time();
          pt.put_funklet(name,pt1.get_funklet(name,idom));
        dprintf(2,"elapsed time: %f seconds\n",time.time()-t0); t0 = time.time();
        pt1 = None;
        self._make_axis_index();
      finally:
        self._end_progress();

  def close (self):
    """Detaches from FastParmTable object (needed when interacting with a meqserver). It is usually a good idea
    to detach from a aprmtable when not using it for a while (i.e. only reattach for the duration of a lengthy
    operation.) Table will be reopened on a subsequent call to parmtable() below.""";
    if self._pt:
      self._pt = None;

  def parmtable (self,write=False):
    """Returns FastParmTable object, reopening it if needed. Use write=True to reopen for writing.""";
    if not self.filename:
      return None;
    if self._pt and write and not self._pt_write:
      self._pt = None;
    if not self._pt:
      self._pt = FastParmTable(self.filename,write);
      self._pt_write = write;
    return self._pt;

  def funklet_names (self):
    """Returns list of funklet names""";
    return self._funklet_names;

  def funklet_name_components (self):
    """Returns list of sets of funklet name components""";
    return self._name_components;

  def envelope_domain (self):
    """Returns domain object which envelops all of the parmtable""";
    kw = {};
    for iaxis,stats in enumerate(self._axis_stats):
      if not stats.empty():
        kw[mequtils.get_axis_id(iaxis)] = stats.minmax;
    return meq.gen_domain(**kw);

  def envelope_cells (self,**num_cells):
    """Returns cells object which envelops all of the parmtable, and is regularly spaced. The number of points
    along each axis is equal to the number of subdomains along that axis, but the cells do not necessarily follow
    the structure of the subdomain if the subdomains are overlapping, spaced out, or irregular.
    """
    dom = self.envelope_domain();
    kw = {};
    for iaxis,stats in enumerate(self._axis_stats):
      if not stats.empty():
        kwname = 'num_'+str(mequtils.get_axis_id(iaxis)).lower();
        kw[kwname] = num_cells.get(kwname,len(stats.cells));
    return meq.gen_cells(dom,**kw);

  def subdomain_cells (self):
    """Returns cells object which envelops all of the parmtable, and contains a cells matching every subdomain.""";
    dom = self.envelope_domain();
    kw = {};
    cells = meq.gen_cells(dom);
    for iaxis,stats in enumerate(self._axis_stats):
      if not stats.empty():
        grid = list(stats.cells.keys());
        grid.sort();
        meq.add_cells_axis(cells,mequtils.get_axis_id(iaxis),grid=grid,cell_size=[
              stats.cells[x0] for x0 in grid ]);
    return cells;

  def axis_stats (self,iaxis):
    """Returns _AxisStats object for the specified parmtable.""";
    return self._axis_stats[iaxis];

  def _make_axis_index (self):
    """Builds up various indices based on content of the parmtable""";
    # check if cache is up-to-date
    cachepath = os.path.join(self.filename,'ParmTab.cache');
    funkpath = os.path.join(self.filename,'funklets');
    self.mtime = os.path.getmtime(funkpath) if os.path.exists(funkpath) else time.time();
    try:
      has_cache = os.path.getmtime(cachepath) >= self.mtime;
      if not has_cache:
        dprintf(2,"cache is out of date, will regenerate\n");
    except:
      dprintf(0,"%s: os.path.getmtime() throws exception, assuming cache is out of date\n",self.filename);
      has_cache = False;
    # try to load the cache if so
    t0 = time.time();
    if has_cache:
      try:
        dprintf(2,"loading index cache\n");
        self._funklet_names,self._domain_list,self._axis_stats,self._name_components, \
        self._domain_fullset,self._domain_cell_index,self._domain_reverse_index \
          = pickle.load(open(cachepath, "rb"));
        dprintf(2,"elapsed time: %f seconds\n",time.time()-t0); t0 = time.time();
        return;
      except:
        if verbosity.get_verbose() > 0:
          traceback.print_exc();
        dprintf(0,"%s: error reading cached stats, regenerating\n",self.filename);
        has_cache = False;
    # no cache, so regenerate everything
    if not has_cache:
      self._axis_stats = [ _AxisStats(mequtils.get_axis_id(i)) for i in range(mequtils.max_axis) ];
      pt = self.parmtable();
      dprintf(2,"loading domain list\n");
      self._domain_list = pt.domain_list();
      dprintf(2,"elapsed time: %f seconds\n",time.time()-t0); t0 = time.time();
      dprintf(2,"collecting axis stats\n");
      self._axes = {};
      for domain in self._domain_list:
        for axis,rng in domain.items():
          if str(axis) != 'axis_map':
            self._axis_stats[mequtils.get_axis_number(axis)].add_cell(*rng);
      dprintf(2,"elapsed time: %f seconds\n",time.time()-t0); t0 = time.time();
      dprintf(2,"finalizing axis stats\n");
      self._domain_fullset = [None]*mequtils.max_axis;
      for iaxis,stats in enumerate(self._axis_stats):
        stats.update();
        if not stats.empty():
          self._domain_fullset[iaxis] = list(range(len(stats.cells)));
          dprintf(2,"axis %s: %d unique cells from %g to %g\n",stats.name,len(stats.cells),*stats.minmax);
      dprintf(2,"elapsed time: %f seconds\n",time.time()-t0); t0 = time.time();
      dprintf(2,"making subdomain indices\n");
      # now make a subdomain index
      self._domain_cell_index = [0]*len(self._domain_list);
      self._domain_reverse_index = {};
      for idom,domain in enumerate(self._domain_list):
        index = [None]*mequtils.max_axis;
        for axis,rng in domain.items():
          if str(axis) != 'axis_map':
            iaxis = mequtils.get_axis_number(axis);
            index[iaxis] = self._axis_stats[iaxis].lookup_cell(*rng);
        # insert into domain_cells_index and domain_reverse_index
        index = tuple(index);
        self._domain_cell_index[idom] = index;
        self._domain_reverse_index[index] = idom;
      dprintf(2,"elapsed time: %f seconds\n",time.time()-t0); t0 = time.time();

      dprintf(2,"loading funklet name list\n");
      self._funklet_names = list(pt.name_list());
      dprintf(2,"elapsed time: %f seconds\n",time.time()-t0); t0 = time.time();
      dprintf(2,"computing funklet indices\n");
      self._name_components = {};
      for name in self._funklet_names:
        for i,token in enumerate(name.split(':')):
          self._name_components.setdefault(i,set()).add(token);
      self._name_components = [ self._name_components[i] for i in range(len(self._name_components)) ];
      for i,values in enumerate(self._name_components):
        dprintf(2,"component %d: %s\n",i,' '.join(values));
      dprintf(2,"elapsed time: %f seconds\n",time.time()-t0); t0 = time.time();

      dprintf(2,"writing cache\n");
      try:
        pickle.dump((
            self._funklet_names,self._domain_list,self._axis_stats,self._name_components, \
            self._domain_fullset,self._domain_cell_index,self._domain_reverse_index \
          ),open(cachepath,'wb')
        );
      except:
        if verbosity.get_verbose() > 0:
          traceback.print_exc();
        dprintf(0,"%s: error writing stats to cache, will probably regenerate next time\n",self.filename);
      dprintf(2,"elapsed time: %f seconds\n",time.time()-t0); t0 = time.time();


  def resolve_output_table (self,outtab,new=False):
    """Helper function. Resolves outtab to a ParmTab object as follows:
    * if it is a ParmTab, returns it directly
    * if it is a filename, creates ParmTab for that table
    * If it is None, creates a ParmTab based on the current table.
    """;
    if isinstance(outtab,ParmTab):
      return outtab;
    elif isinstance(outtab,str):
      return ParmTab(outtab,write=True,new=new);
    else:
      name,ext = os.path.splitext(self.filename);
      return ParmTab(name+"_out"+ext,write=True,new=new);

  def make_slicing (self,slicing):
    if not slicing:
      return DomainSlicing([],self._domain_fullset);
    elif isinstance(slicing,DomainSlicing):
      return slicing;
    elif isinstance(slicing,(str,int)):
      slicing = [slicing];
    return DomainSlicing(slicing,self._domain_fullset);

  def funkset (self,name):
    return FunkSet(self,name);

  def apply (self,op_func,slicing,outtab=None,remove=False,newtab=False):
    """For each funklet in our table, takes all funklets along the designated slicing axis
    (i.e. for each slice along the non-listed axes), creates a FunkSlice, and calls
    op_func(slice).
    The return value of func() should be a funklet list, which is written to the output table,
    or None to ignore.
    if 'outtab' is None, a new output table is created. Otherwise set outtab to a filename, or a
    ParmTab, or a FastParmTable.
    If 'remove' is True, input funklets will be removed if an output funklet is returned.
    """;
    t0 = time.time();
    nn = len(self.funklet_names());
    self._start_progress("applying operation '%s'"%op_func.__name__,nn);
    try:
      outtab = self.resolve_output_table(outtab,newtab);
      dprintf(1,"using output table %s\n",outtab.filename);
      pt = self.parmtable();
      # loop over funklets and slices
      dprintf(3,"input slicing is %s\n",slicing);
      slicing = self.make_slicing(slicing);
      dprintf(2,"%d slices will be iterated over\n",len(slicing));
      for iname,name in enumerate(sort_qualified_names(self.funklet_names())):
        self._report_progress(iname);
        self.funkset(name).apply(op_func,slicing,outtab=outtab,remove=remove);
    finally:
      dprintf(2,"elapsed time: %f seconds\n",time.time()-t0); t0 = time.time();
      self._end_progress(False);
      self.close();


def open (*args,**kw):
  """Opens a ParmTab. Arguments are passed to ParmTab constructor:

  """+ParmTab.__init__.__doc__;
  return ParmTab(*args,**kw);



if __name__ == '__main__':
  import sys
  from . import FunkOps
  if len(sys.argv) < 2:
    print("Pass in an fmep table to read its stats");
  else:
    verbose(3);
    pt = ParmTab(sys.argv[1]);
    if '-average' in sys.argv:
      pt.apply(FunkOps.average,"time",newtab=True);
    if '-interpol' in sys.argv:
      pt.apply(FunkOps.linear_interpol,"time",newtab=True);
    if '-rank0' in sys.argv:
      pt.apply(FunkOps.force_rank0,["time","freq"],newtab=True);
    if '-infdom' in sys.argv:
      pt.apply(FunkOps.make_infinite_domain,["time","freq"],newtab=True);

