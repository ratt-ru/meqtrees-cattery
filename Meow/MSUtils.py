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

import re
import traceback
import sys
import os
import os.path
import Meow
import sets
import math

# figure out which table implementation to use -- try pyrap/casacore first
try:
  import pyrap_tables
  TABLE = pyrap_tables.table
  print "Meow.MSUtils: using the pyrap_tables module"
except:
  try:
    import pyrap.tables
    TABLE = pyrap.tables.table
    print "Meow.MSUtils: using the pyrap.tables module"
  except:
    # else try the old pycasatable/aips++ thing
    try:
      import pycasatable
      TABLE = pycasatable.table
      print "Meow.MSUtils: using the pycasatable module. WARNING: this is deprecated."
      print "Please install pyrap and casacore!"
    except:
      TABLE = None;
      print "Meow.MSUtils: no tables module found, GUI functionality will be reduced"
      print "Please install pyrap and casacore!"

def find_exec (execname):
  path = os.environ.get('PATH') or os.defpath;
  for dirname in path.split(os.pathsep):
    fname = os.path.join(dirname,execname);
    if os.access(fname,os.R_OK|os.X_OK):
      return fname;
  return None;

# figure out if we have an imager
_lwimager = find_exec('lwimager');
if _lwimager:
  print "Meow.MSUtils: found %s, can use it for imaging"%_lwimager;

_glish = find_exec('glish');
if _glish:
  print "Meow.MSUtils: found %s, can use AIPS++ imager"%_glish;

if not _lwimager and not _glish:
  print "Meow.MSUtils: no imager found";
  
# figure out if we have a visualizer
#_VISUALIZER = None;
#if os.system('which kvis >/dev/null') == 0:
#  _VISUALIZER = 'kvis';
#  print "Meow.MSUtils: found kvis";


# queue size parameter for MS i/o record
ms_queue_size = 500;

FLAGBITS = range(31);
FLAG_ADD = "add to set";
FLAG_REPLACE = "replace set";
FLAG_REPLACE_ALL = "replace all sets";
FLAG_FULL = 0x7FFFFFFF;

def longest_prefix (*strings):
  """helper function, returns longest common prefix of a set of strings""";
  if not strings:
    return '';
  if len(strings) == 1:
    return strings[0];
  # sort by length
  lenlist = [ (len(s),s) for s in strings ];
  lenlist.sort();
  strings = [ pair[1] for pair in lenlist ];
  for i,ch in enumerate(strings[0]):
    for s1 in strings[1:]:
      if s1[i] != ch:
        return s1[:i];
  return strings[0];

class MSContentSelector (object):
  def __init__ (self,ddid=[0],field=None,channels=True,namespace='ms_sel'):
    """Creates options for selecting a subset of an MS.
    ddid:       list of suggested DDIDs, or false for no selector.
    field:      list of suggested fields. or false for no selector.
      NB: if TABLE is available, ddid/field is ignored, and selectors are always
      provided, based on the MS content.
    namespace:  the TDLOption namespace name, used to qualify TDL options created here.
        If making multiple selectors, you must give them different namespace names
    """;
    self.tdloption_namespace = namespace;
    self._opts = [];
    self.ddid_index = self.field_index = None;
    self.ms_spws = self.ms_field_names = self.ms_ddid_numchannels = None;
    # field/ddid selectors
    more_ids = int;
    if TABLE:
      ddid = [0];
      field = [0];
      more_ids = None;  # no more= option to ddid/field: they come from the MS
    if ddid:
      self.ddid_option = TDLOption('ddid_index',"Data description ID",
        ddid,more=more_ids,namespace=self,
        doc="""If the MS contains multiple spectral windows, etc., then use this option to select different DATA_DESCRIPTION_IDs. Default is 0.""");
      self._opts.append(self.ddid_option);
    if field:
      self.field_option = TDLOption('field_index',"Field ID",field,
        more=more_ids,namespace=self,
        doc="""If the MS contains multiple fields, then use this option to select different FIELD_IDs. Default is 0.""");
      self._opts.append(self.field_option);
    # channel selector
    if channels:
      self.channel_options = [
        TDLOption('ms_channel_start',"First channel",[0],more=int,namespace=self),
        TDLOption('ms_channel_end',"Last channel",[0],more=int,namespace=self),
        TDLOption('ms_channel_step',"Channel stepping",[1,2],more=int,namespace=self)
      ];
      chanmenu = TDLMenu("Channel selection",
                         toggle='select_channels',namespace=self,
                         *self.channel_options);
      self._opts.append(chanmenu);
    else:
      self.select_channels = False;
    # additional taql string
    self.taql_option = TDLOption('ms_taql_str',"Additional TaQL selection",
                                [None],more=str,namespace=self)
    self._opts.append(self.taql_option);
    self._nchan = None;
    # if TABLE exists, set up interactivity for ddid/channels
    if TABLE:
      if channels:
        self.channel_options[0].set_validator(self._validate_first_channel);
        self.channel_options[1].set_validator(self._validate_last_channel);
        self.channel_options[2].set_validator(self._validate_channel_step);
      self.ddid_option.when_changed(self._select_ddid);

  def option_list (self):
    """Returns list of all TDL options"""
    return self._opts;

  def get_channels (self):
    """Returns start,end,step tuple for selected channels, or None
    if we were created without a channel selector""";
    if self.select_channels:
      if self.ms_channel_start <= self.ms_channel_end:
        return self.ms_channel_start,self.ms_channel_end,self.ms_channel_step;
      else:
        return self.ms_channel_end,self.ms_channel_start,self.ms_channel_step;
    else:
      return None;
    
  def get_total_channels (self):
    """Returns total # of channels in current spectral window, or None
    if not known""";
    return self._nchan;

  def get_spectral_window (self):
    """Returns the spw corresponding to the selected DDID. If pycasatable is n/a,
    simply returns the ddid_index (kludgy kludgy)
    """;
    ddid = self.ddid_index or 0;
    if self.ms_spws:
      return self.ms_spws[ddid];
    else:
      return ddid;

  def get_field (self):
    """Returns the field number""";
    return self.field_index or 0;
  
  def get_ddid (self):
    """Returns the field number""";
    return self.ddid_index or 0;
  
  def get_taql_string (self):
    return self.ms_taql_str or '';

  def create_selection_record (self):
    """Forms up a selection record that can be added to an MS input record""";
    selection = record();
    if self.select_channels:
      selection.channel_start_index = min(self.ms_channel_start,self.ms_channel_end);
      selection.channel_end_index = max(self.ms_channel_start,self.ms_channel_end);
      selection.channel_increment = self.ms_channel_step;
    if self.ddid_index is not None:
      selection.ddid_index = self.ddid_index;
    if self.field_index is not None:
      selection.field_index = self.field_index;
    if self.ms_taql_str:
      selection.selection_string = self.ms_taql_str;
    return selection;

  def _select_new_ms (self,ms):
    """Called (from MSSelector) when a new MS is selected. ms is a pycasatable.table
    object. Fills ddid/field/channel selectors from the MS.
    """;
    # DDIDs
    self.ms_spws = list(TABLE(ms.getkeyword('DATA_DESCRIPTION'),
                                            lockoptions='autonoread') \
                                          .getcol('SPECTRAL_WINDOW_ID'));
    numchans = TABLE(ms.getkeyword('SPECTRAL_WINDOW'),
                                  lockoptions='autonoread').getcol('NUM_CHAN');
    self.ms_ddid_numchannels = [ numchans[spw] for spw in self.ms_spws ];
    # Fields
    field = TABLE(ms.getkeyword('FIELD'),lockoptions='autonoread');
    self.ms_field_names = list(field.getcol('NAME'));
    self.ms_field_phase_dir = field.getcol('PHASE_DIR');
    # update selectors
    self._update_ms_options();

  def _update_ms_options (self):
    # DDID selector
    if self.ms_ddid_numchannels is not None:
      self.ddid_option.set_option_list(range(len(self.ms_ddid_numchannels)));
      if len(self.ms_ddid_numchannels) < 2:
        self.ddid_option.hide();
    # Field selector
    if self.ms_field_names is not None:
      field_list = [ "%d:%s"%x for x in enumerate(self.ms_field_names) ];
      self.field_option.set_option_list(dict(enumerate(field_list)));
      if len(self.ms_field_names) < 2:
        self.field_option.hide();

  def _update_from_other (self,other):
    self.ms_spws = other.ms_spws;
    self.ms_ddid_numchannels = other.ms_ddid_numchannels;
    self.ms_field_names = other.ms_field_names;
    self._update_ms_options();

  def _select_ddid (self,value):
    """callback used when a new DDID is selected."""
    if self.ms_ddid_numchannels:
      nchan = self._nchan = self.ms_ddid_numchannels[value];
      self.channel_options[0].set_option_list([0,nchan-1]);
      self.channel_options[0].set_value(min(self.ms_channel_start,nchan-1),save=False);
      self.channel_options[0].set_doc("(max %d) select first channel"%(nchan-1));
      self.channel_options[1].set_option_list([0,nchan-1]);
      self.channel_options[1].set_value(min(self.ms_channel_end,nchan-1),save=False);
      self.channel_options[1].set_doc("(max %d) select last channel"%(nchan-1));
      self.channel_options[2].set_value(min(self.ms_channel_step,nchan),save=False);
      self.channel_options[2].set_doc("(max %d) select channel steping"%nchan);

  def _validate_first_channel (self,value):
    """validator for channel selectors."""
    if self.ms_ddid_numchannels is None:
      return True;
    return isinstance(value,int) and \
           value >= 0 and \
           value < self.ms_ddid_numchannels[self.ddid_index];

  def _validate_last_channel (self,value):
    """validator for channel selectors."""
    if self.ms_ddid_numchannels is None:
      return True;
    return isinstance(value,int) and \
          value >= 0 and \
          value < self.ms_ddid_numchannels[self.ddid_index];

  def _validate_channel_step (self,value):
    """validator for channel selectors."""
    if self.ms_ddid_numchannels is None:
      return True;
    return isinstance(value,int) and \
          value >= 1 and \
          value <= self.ms_ddid_numchannels[self.ddid_index];

class MSFlagSelector (object):
  """Abstract base class for flags electors""";
  def __init__ (self,namespace):
    """Creates options for selecting a subset of flags in an MS.
    legacy: if True, options to include legacy flags will be provided.
    namespace:  the TDLOption namespace name, used to qualify TDL options created here.
        If making multiple selectors, you must give them different namespace names
    """;
    self.tdloption_namespace = namespace;
    self._opts = [];
    self.flagsets = None;
    self.bitflag_labels = [];
    self.bitflag_bits   = [ 1<<bit for bit in FLAGBITS ];
    
  def option_list (self):
    """Returns list of all TDL options"""
    return self._opts;
    
  def update_flagsets (self,flagsets):
    self.flagsets = flagsets;
    flagsets.when_changed(self.update);
    
  def update(self):
    pass;

class MSReadFlagSelector (MSFlagSelector):
  def __init__ (self,legacy=True,namespace='ms_rfl',
    label="Include flagset %s",
    doc="Include this flagset in the overall data flags."):
    """Creates options for selecting a subset of flags in an MS.
    legacy: if True, options to include legacy flags will be provided.
    namespace:  the TDLOption namespace name, used to qualify TDL options created here.
        If making multiple selectors, you must give them different namespace names
    """;
    MSFlagSelector.__init__(self,namespace);
    self._label = label;
    # option to read legacy flags
    if legacy:
      self.read_legacy_flag_opt = TDLOption('read_legacy_flags',
          "Include legacy FLAG column",True,namespace=self,
          doc="""If enabled, then data flags in the standard AIPS++ FLAG column will
          be included in the overall set of data flags. 
          Normally you would only work with flagsets and ignore the FLAG column.
          """);
      self._opts.append(self.read_legacy_flag_opt);
    else:
      self.read_legacy_flags = False;
      self.read_legacy_flag_opt = None;
    self.read_bitflag_opts = [ 
        TDLOption('read_bitflag_%d'%bit,
                  "MS bitflag %d"%bit,True,namespace=self,
                  doc=doc) 
      for bit in FLAGBITS ];
    for opt in self.read_bitflag_opts:
      opt.hide(); 
    self._opts += self.read_bitflag_opts;
    
  def update (self):
    self.bitflag_labels = self.flagsets.names() or [];
    self.bitflag_bits   = [ self.flagsets.flagmask(name) for name in self.bitflag_labels ];
    nlab = min(len(self.read_bitflag_opts),len(self.bitflag_labels));
    for opt,label in zip(self.read_bitflag_opts,self.bitflag_labels):
      opt.set_name(self._label%label);
      opt.show();
    for opt in self.read_bitflag_opts[len(self.bitflag_labels):]:
      opt.hide();
      
  def get_flagmask (self):
    flagmask = 0;
    for bit in range(len(self.bitflag_bits)):
      if getattr(self,'read_bitflag_%d'%bit,True):
        flagmask |= self.bitflag_bits[bit];
    return flagmask;
      
  def selected_flagsets (self):
    return [ fs for i,fs in enumerate(self.bitflag_labels) 
              if getattr(self,'read_bitflag_%d'%i,True) ];
  

class MSWriteFlagSelector (MSFlagSelector):
  def __init__ (self,namespace='ms_wfl'):
    """Creates options for selecting an output flagset of an MS.
    namespace:  the TDLOption namespace name, used to qualify TDL options created here.
        If making multiple selectors, you must give them different namespace names
    """;
    MSFlagSelector.__init__(self,namespace);
    if TABLE is not None:
      # For the output, either we already have a BITFLAG column, in which
      # case we'll get a list of bitflag labels from there when the MS is selected,
      # or we don't, in which case we need only provide one default name
      self.write_bitflag_opt = TDLOption("write_bitflag","Output flagset",
                                          ["FLAG0"],more=str,namespace=self,
        doc="""If new data flags are generated within the tree, they may be stored in an existing
        or new flagset. Select the output flagset here""");
    else:
      # No table support, so we'll just ask for a bit number
      self.write_bitflag_opt = TDLOption("write_bitflag","Output bitflag",
                                          range(16),more=int,namespace=self,
        doc="""If new data flags are generated within the tree, they may be stored in an existing
        or new flagset. Select the output flagset here by specifying a bit number.""");
    self._has_new_bitflag_label = self.flagsets = None;
    self._opts.append(self.write_bitflag_opt);
    
  def update (self):
    self.bitflag_labels = self.flagsets.names() or [];
    self.bitflag_bits   = [ self.flagsets.flagmask(name) for name in self.bitflag_labels ];
    self.write_bitflag_opt.set_option_list(self.bitflag_labels);
      
  def get_flagmask (self):
    bitflag = self.write_bitflag;
    if isinstance(bitflag,str):
      return self.flagsets.flagmask(bitflag,create=True);
    else:
      return 1<<bitflag;



class MSSelector (object):
  _corr_1 = "1";
  _corr_2 = "2";
  _corr_2x2 = "2x2";
  _corr_2x2_diag = "2x2, diagonal terms only"
  _corr_2x2_offdiag = "2x2, off-diagonal terms only"
  _corr_index = {
    _corr_1:[0],
    _corr_2:[0,-1,-1,1],
    _corr_2x2:[0,1,2,3],
    _corr_2x2_diag:[0,-1,-1,3],
    _corr_2x2_offdiag:[-1,1,2,-1],
  };
  ms_corr_names = [ "XX","XY","YX","YY" ];
  
  """An MSSelector implements TDL options for selecting an MS and a subset therein""";
  def __init__ (self,
                pattern="*.ms *.MS",
                has_input=True,
                has_output=True,
                forbid_output=["DATA"],
                antsel=True,
                tile_sizes=[1,5,10,20,30,60],
                ddid=[0],
                field=[0],
                channels=True,
                flags=False,read_flags=False,write_flags=False,write_legacy_flags=False,
                hanning=False,invert_phases=False,
                namespace='ms_sel'
                ):
    """Creates an MSSelector object
    filter:     ms name filter. Default is "*.ms *.MS"
    has_input:  is an input column selector initially enabled.
    has_input:  is an output column selector initially enabled.
    forbid_output: a list of forbidden output columns. "DATA" by default"
    antsel:     if True, an antenna subset selector will be provided
    tile_sizes: list of suggested tile sizes. If false, no tile size selector is provided.
    ddid:       list of suggested DDIDs, or false for no selector.
    field:      list of suggested fields. or false for no selector.
      NB: if TABLE is available, ddid/field is ignored, and selectors are always
      provided, based on the MS content.
    channels:   if True, channel selection will be provided
    flags:      if flags or write_flags, a "write flags" option will be provided.
                if flags or read_flags, a "read flags" option will be provided.
                if flags or write_legacy_flags, a "fill legacy flags" option will be provided.
    hanning:    if True, an apply Hanning tapering option will be provided
    invert_phases: if True, an invert phases option will be provided
    namespace:  the TDLOption namespace name, used to qualify TDL options created here.
        If making multiple selectors, you must give them different namespace names.
    """;
    self.tdloption_namespace = namespace;
    self._content_selectors = [];
    self.ms_antenna_names = [];
    self.ms_antenna_sel = self.antsel_option = None;
    ms_option = self._ms_option = \
      TDLOption('msname',"MS",TDLDirSelect(pattern,default=True),namespace=self);
    self._compile_opts = [ ms_option ];
    self._opts = [];
    # antenna selector
    if antsel:
      if isinstance(antsel,str):
        antsel = [None,antsel];
      elif isinstance(antsel,(tuple,list)):
        antsel = [None]+antsel;
      else:
        antsel = [None];
      self.antsel_option = TDLOption("ms_antenna_sel","Antenna subset",
                                     [None],more=str,namespace=self,
        doc="""Selects a subset of antennas to use. You may specify individual indices
        (1-based) separated by commas or spaces, or ranges, e.g. "M:N" (M to N inclusive),
        or ":M" (1 to M), or "N:" (N to last).
        Example subset: ":3 5 8 10:12 16:"."""
      );
      self.antsel_option.set_validator(self._antenna_sel_validator);
      # hide until an MS is selected
      if TABLE:
        self.antsel_option.hide();
      self._compile_opts.append(self.antsel_option);
    # correlation options
    self.corrsel_option = TDLOption("ms_corr_sel","Correlations",
                                    [self._corr_2x2,self._corr_2x2_diag,
                                    self._corr_2,self._corr_1],
                                    namespace=self);
    self._compile_opts.append(self.corrsel_option);
    self.ms_data_columns = ["DATA","MODEL_DATA","CORRECTED_DATA"];
    if isinstance(forbid_output,str):
      self._forbid_output = sets.Set([forbid_output]);
    elif forbid_output:
      self._forbid_output = sets.Set(forbid_output);
    else:
      self._forbid_output = [];
    self.input_column = self.output_column = None;
    self.ms_has_input = has_input;
    self.ms_has_output = has_output;
    self.ms_antenna_positions = None;
    # if no access to tables, then allow more input columns to be entered
    if TABLE:
      more_col = None;
    else:
      more_col = str;
    self.input_col_option = TDLOption('input_column',"Input MS column",
                                      self.ms_data_columns,
                                      namespace=self,more=more_col);
    self.input_col_option.show(has_input);
    self._opts.append(self.input_col_option);
    self.output_col_option = TDLOption('output_column',"Output MS column",
                                      self.ms_data_columns + [None],
                                      namespace=self,more=str,default=2);
    self.output_col_option.show(has_output);
    self._opts.append(self.output_col_option);
    # tile sizes
    if tile_sizes:
      self._opts.append(TDLOption('tile_size',"Tile size (timeslots)",
                                  tile_sizes,more=int,namespace=self));
    else:
      self.tile_size = 1;
    if hanning:
      self._opts.append(TDLOption('ms_apply_hanning',"Apply Hanning taper to input",
                                  False,namespace=self));
    else:
      self.ms_apply_hanning = None;
    if invert_phases:
      self._opts.append(TDLOption('ms_invert_phases',"Invert phases in input",
                                  False,namespace=self));
    else:
      self.ms_invert_phases = None;
    # add a default content selector
    self._ddid,self._field,self._channels = ddid,field,channels;
    self.subset_selector = self.make_subset_selector(namespace);
    self._opts += self.subset_selector.option_list();
    # setup flag-related options
    self._flag_selectors = [];
    self.flagsets = None;
    if flags or read_flags:
      self.read_flag_selector = self.make_read_flag_selector('ms_rfl',legacy=True);
      # add into submenu
      self.read_flags_opt = TDLMenu("Read flags from MS",
        toggle="ms_read_flags",default=True,namespace=self,open=False,
        doc="""If the Measurement Set contains data flags, enable this option to propagate the flags
        into the tree. Flagged data is (normally) ignored in all calculations.
        There may exist multiple sets of data flags. One set is stored in the legacy FLAG column
        defined by AIPS++/casa. With MeqTrees, you may create a number of additional flagsets with
        arbitrary labels. Via the suboptions within this menu, you may choose to include or exclude 
        each flagset from the overall set of active data flags.
        """,
        *self.read_flag_selector.option_list()
        );
      self._opts.append(self.read_flags_opt);
    else:
      self.ms_read_flags = self.read_flag_selector = None;
    if flags or write_flags:
      self.write_flag_selector = self.make_write_flag_selector('ms_wfl');
      self.write_flags_opt = TDLMenu("Write output flagset",
        toggle='ms_write_flags',default=False,namespace=self,
        doc="""If new data flags are generated within the tree, they may be written out to the
        Measurement Set. Enable this option to write flags to the MS.""",
        *( self.write_flag_selector.option_list() +
           [ TDLOption("ms_write_flag_policy",
                    "Output flagset policy",[FLAG_ADD,FLAG_REPLACE],namespace=self,
                    doc="""Flags generated within the tree may be added to an existing flagset, or may
                    replace a flagset.""")]
         )
        );
      self._opts.append(self.write_flags_opt);
    else:
      self.ms_write_flags = False;
      self.ms_write_flag_selector = None;
      self.ms_write_replace_flags = False;
    # add callbacks
    self._when_changed_callbacks = [];
    # if TABLE exists, set up interactivity for MS options
    if TABLE:
      ms_option.set_validator(self._select_new_ms);
      
  def when_changed (self,callback):
    # if tables are available, callbacks will be called by _select_new_ms()
    if TABLE:
      self._when_changed_callbacks.append(callback);
      # and call it immediately if we already have an MS
      if getattr(self,'_msname',None):
        callback(self._msname);
    # otherwise just register the callback with the ms_option itself
    else:
      self._ms_option.when_changed(callback);

  def enable_input_column (self,enable=True):
    self.ms_has_input = enable;
    self.input_col_option.show(enable);

  def enable_output_column (self,enable=True):
    self.ms_has_output = enable;
    self.output_col_option.show(enable);

  def option_list (self):
    """Returns list of all TDL options. Note that the MS name is always
    the first option."""
    return self._compile_opts + self._opts;

  def compile_options (self):
    """Returns MS name option."""
    return self._compile_opts;

  def runtime_options (self):
    """Returns list of all TDL options. Note that the MS name is always
    the first option."""
    return self._opts;

  def get_antenna_names (self):
    """Returns the list of antenna names from the current MS. If pycasatable
    is not available, this will be empty""";
    if self.ms_antenna_sel:
      subset = parse_antenna_subset(self.ms_antenna_sel,self.ms_antenna_names);
      if self.ms_antenna_names:
        return [self.ms_antenna_names[index] for index in subset];
      else:
        return subset;
    else:
      return self.ms_antenna_names;
    
  def get_phase_dir (self):
    """Returns the phase direction of the currently selected pointing, or None if none
    is available."""
    # if no field info available, or more than 1 field in MS, then return None
    if not self.subset_selector.ms_field_names or len(self.subset_selector.ms_field_names) > 1:
      return None;
    # otherwise, return phase dir of first and only field
    return self.subset_selector.ms_field_phase_dir[0,0,0],self.subset_selector.ms_field_phase_dir[0,0,1];

  def get_antenna_set (self,default=None):
    """Returns the set of selected antenna indices from the current MS, or None
    if no selection info is available."""
    subset = self.ms_antenna_sel and parse_antenna_subset(self.ms_antenna_sel,self.ms_antenna_names);
    if subset:
      if self.ms_antenna_names:
        return [(index,self.ms_antenna_names[index]) for index in subset];
      else:
        return subset;
    elif self.ms_antenna_names:
      return list(enumerate(self.ms_antenna_names));
    else:
      return default;
    
  def get_corr_index (self):
    """Returns the set of selected antenna correlation indices
    """;
    return self._corr_index[self.ms_corr_sel];
  
  def get_correlations (self):
    return [ self.ms_corr_names[icorr] for icorr in self.get_corr_index() if icorr >= 0 ];

  def make_subset_selector (self,namespace,**kw):
    """Makes an MSContentSelector object connected to this MS selector."""
    sel = MSContentSelector(ddid=self._ddid,field=self._field,
                            channels=self._channels,namespace=namespace,**kw);
    # extra selectors needs to be initialized with existing info
    if TABLE and self._content_selectors:
      sel._update_from_other(self.subset_selector);
      sel._select_ddid(self.subset_selector.ddid_index or 0);
    self._content_selectors.append(sel);
    return sel;
    
  def make_read_flag_selector (self,namespace,**kw):
    """Makes an MSReadFlagSelector object connected to this MS selector."""
    sel = MSReadFlagSelector(namespace=namespace,**kw);
    # extra selectors needs to be initialized with existing info
    if self.flagsets:
      sel.update_flagsets(self.flagsets);
    self._flag_selectors.append(sel);
    return sel;
    
  def make_write_flag_selector (self,namespace,**kw):
    """Makes an MSWriteFlagSelector object connected to this MS selector."""
    sel = MSWriteFlagSelector(namespace=namespace,**kw);
    # extra selectors needs to be initialized with existing info
    if self.flagsets:
      sel.update_flagsets(self.flagsets);
    self._flag_selectors.append(sel);
    return sel;
    
  def get_input_flagmask (self):
    if not self.ms_read_flags:
      return 0;
    return self.read_flag_selector.get_flagmask();
      
  def get_output_bitflag (self):
    if not self.ms_write_flags:
      return 0;
    return self.write_flag_selector.get_flagmask();
  
  def reload (self):
    msname = getattr(self,'_msname',None);
    if msname:
      self._msname = None;
      self._select_new_ms(msname);
    
  def _select_new_ms (self,msname):
    """This callback is called whenever a new MS is selected. Returns False if
    table is malformed or n/a""";
    # do nothing if already read this MS
    if msname == getattr(self,'_msname',None):
      return True;
    try:
      ms = TABLE(msname,lockoptions='autonoread');
      # data columns
      self.ms_data_columns = [ name for name in ms.colnames() if name.endswith('DATA') ];
      self.input_col_option.set_option_list(self.ms_data_columns);
      outcols = [ col for col in self.ms_data_columns if col not in self._forbid_output ];
      self.output_col_option.set_option_list(outcols);
      # antennas
      anttable = TABLE(ms.getkeyword('ANTENNA'),lockoptions='autonoread');
      antnames = anttable.getcol('NAME');
      # if NAME column is missing, use indices
      if not antnames:
        print "Warning! This MS does not define ANTENNA names. Using antenna indices instead.";
        self.ms_antenna_names = map(str,range(anttable.nrows()));
      # else use name, but trim off longest common prefix (so that RT0,RT1,..RTD become 0,1,...,D
      else:
        prefix = len(longest_prefix(*antnames));
        self.ms_antenna_names = [ name[prefix:] for name in antnames ];
        # some broken MSs do not have unique antenna names -- replace them with indices if so
        if len(sets.Set(self.ms_antenna_names)) < len(self.ms_antenna_names):
          print "Warning! This MS does not define unique ANTENNA names. Using antenna indices instead.";
          self.ms_antenna_names = [ str(i) for i in range(len(self.ms_antenna_names)) ];
      if self.antsel_option:
        self.antsel_option.set_option_list([None,"0:%d"%(len(self.ms_antenna_names)-1)]);
        self.antsel_option.show();
      self.ms_antenna_positions = anttable.getcol('POSITION');
      # correlations
      corrs = TABLE(ms.getkeyword('POLARIZATION'),
                    lockoptions='autonoread').getcol('CORR_TYPE');
      ncorr = len(corrs[0]);
      if ncorr < 2:
        corrlist = [self._corr_1];
      elif ncorr < 4:
        corrlist = [self._corr_2,self._corr_1];
      else:
        corrlist = [self._corr_2x2,self._corr_2x2_diag,self._corr_2,self._corr_1];
      self.corrsel_option.set_option_list(corrlist);
      # get flagsets and notify flag selectors
      self.flagsets = get_flagsets(ms);
      self.flagsets.load(ms);
      for sel in self._flag_selectors:
        sel.update_flagsets(self.flagsets);
      # notify content selectors
      for sel in self._content_selectors:
        sel._select_new_ms(ms);
      self._msname = msname;
      # notify callbacks
      for cb in self._when_changed_callbacks:
        cb(msname);
      ms.close();
      return True;
    except:
      print "error reading MS",msname;
      traceback.print_exc();
      return False;
      
  def _antenna_sel_validator (self,value):
    try:
      antenna_subset = parse_antenna_subset(value,self.ms_antenna_names);
      return True;
    except:
      print "error parsing antenna selection '%s'"%value;
      traceback.print_exc();
      return False;

  def imaging_selector (self,*args,**kw):
    """Makes an ImagingSelector connected to this MS selector. All arguments
    are passed to the ImagingSelector constructor""";
    return ImagingSelector(self,*args,**kw);

  def create_inputrec (self,tiling=None):
    """Creates an input record with the selected options""";
    if self.msname is None:
      raise ValueError,"Measurement Set not specified";
    rec = record();
    rec.ms_name          = self.msname
    if self.input_column:
      rec.data_column_name = self.input_column;
    tiling = tiling or self.tile_size;
    if isinstance(tiling,(list,tuple)):
      if len(tiling) != 2:
        raise TypeError,"tiling: 2-list or 2-tuple expected";
      (tile_segments,tile_size) = tiling;
      if tile_segments is not None:
        rec.tile_segments    = tile_segments;
      if tile_size is not None:
        rec.tile_size        = tile_size;
    else:
      rec.tile_size = tiling;
    rec.selection = self.subset_selector.create_selection_record();
    if self.ms_apply_hanning is not None:
      rec.apply_hanning = self.ms_apply_hanning;
    if self.ms_invert_phases is not None:
      rec.invert_phases = self.ms_invert_phases;
    rec.flag_mask = 0;
    rec.legacy_bitflag = 0;
    if self.ms_read_flags:
      rec.tile_bitflag = 2;
      if self.read_flag_selector.read_legacy_flags:
        rec.legacy_bitflag = 1;
      rec.flag_mask = self.read_flag_selector.get_flagmask();
    # form top-level record
    iorec = record(ms=rec);
    iorec.python_init = 'Meow.ReadVisHeader';
    iorec.mt_queue_size = ms_queue_size;
    return iorec;

  def create_outputrec (self):
    """Creates an output record with the selected options""";
    rec = record();
    rec.write_bitflag = self.ms_write_flags;
    if self.ms_write_flags:
      # output masks
      rec.tile_bitflag = self.get_output_bitflag();
      rec.tile_flag_mask = FLAG_FULL & ~3; # bitflags 1|2 = input flags
      if self.ms_write_flag_policy == FLAG_ADD:
        rec.ms_flag_mask = FLAG_FULL;
      elif self.ms_write_flag_policy == FLAG_REPLACE:
        rec.ms_flag_mask = FLAG_FULL & ~rec.tile_bitflag;
    if self.ms_has_output and self.output_column:
      rec.data_column = self.output_column;
    return record(ms=rec,mt_queue_size=ms_queue_size);

  def create_io_request (self,tiling=None):
    """Creates an i/o record with the selected options, suitable for
    passing to a VisDataMux""";
    req = meq.request();
    req.input = self.create_inputrec(tiling);
    if self.ms_write_flags or self.output_column is not None:
      req.output = self.create_outputrec();
    return req;

  def run_solve_job (self,mqs,solvables,tiling=None,
                     solver_node="solver",vdm_node="VisDataMux",options=None,
                     wait=False):
    """helper method to run a solution with a bunch of solvables""";
    # set solvables list in solver
    solver_defaults = Meow.Utils.create_solver_defaults(solvables,options=options)
    mqs.setnodestate(mqs,solver_node,solver_defaults,sync=True,wait=wait)
    req = self.create_io_request(tiling);
    mqs.execute(vdm_node,req,wait=wait);

CHANMODE_NOAVG = "all (no averaging)";
CHANMODE_ALL   = "1 (average all)";
CHANMODE_MFS   = "1 (multi-freq synthesis)";
CHANMODE_MANUAL = "select manually";

class ImagingSelector (object):
  """ImagingSelector provides a set of TDL options for imaging""";
  def __init__ (self,mssel,npix=256,arcmin=5,cellsize=None,subset=True,namespace='img_sel'):
    """Creates an ImagingSelector.
    mssel: an MSSelector object
    npix:       default image size in pixels (or list of suggested sizes)
    arcmin:     default image size in arc minutes (or list of suggested sizes)
    cellsize:   default cell size (as an aips++ quantity string), or list. Note
                that either arcmin or cellsize should be specified, not both.
    subset:     if True, a separate MSContentSelector will be provided, allowing
                to image a different subset of the MS.
    namespace:  the TDLOption namespace name, used to qualify TDL options created here.
        If making multiple selectors, you must give them different namespace names.
    """;
    global _lwimager;
    global _glish;
    if not _lwimager and not _glish:
      def do_nothing (mqs,parent,**kw):
        None;
      self._opts = [ TDLJob(do_nothing,"No lwimager or glish found, please install casarest or aips++") ];
      return;
    self._opts = [];
    self.tdloption_namespace = namespace;
    self.mssel = mssel;
    # add imager
    if _lwimager and _glish:
      self._opts.append(TDLOption('imager_type',"Imager to use",["lwimager","AIPS++ imager"],namespace=self));
    elif _lwimager:
      self.imager_type = "lwimager";
    elif _glish:
      self.imager_type = "aips++";
    # add imaging column option
    self.img_col_option = TDLOption('imaging_column',"Image type or column",
                          list(mssel.ms_data_columns)+["psf"],more=str,namespace=self);
    mssel.output_col_option.when_changed(curry(self.img_col_option.set_value,save=False));
    self._opts.append(self.img_col_option);
    
    chan_menu = TDLMenu("Output channel selection",
        TDLOption('imaging_nchan',"Number of output channels",[1],more=int,namespace=self),
        TDLOption('imaging_chanstart',"Starting at channel",[0],more=int,namespace=self),
        TDLOption('imaging_chanstep',"Stepped by",[1],more=int,namespace=self),
    );
    docstr = """This selects the number of frequency channels in the output
                image. If "all" is used, then every selected channel in the
                MS is imaged independently. If "multi-freq synthesis" is selected, then
                all channels are averaged together in MFS mode. Otherwise, supply a
                number of output channels, the input channels will be averaged down
                to this number."""
    self.imaging_totchan = 1;
    if TABLE:
      chan_opt = TDLOption('imaging_chanmode',"Frequency channels in image",
                    [CHANMODE_NOAVG,CHANMODE_ALL,CHANMODE_MFS,CHANMODE_MANUAL],
                    more=int,namespace=self,doc=docstr);
    else:
      self._opts += [ TDLOption('imaging_totchan',"Frequency channels in input data",
                               [1],more=int,namespace=self,\
                      doc="""Since pycasatable (aips++) or pyrap_tables (pyrap) is not 
                      available, we don't know how many total channels there are in the MS.  
                      Please supply the number here.""") ];
      chan_opt = TDLOption('imaging_chanmode',"Frequency channels in image",
                    [CHANMODE_ALL,CHANMODE_MFS,CHANMODE_MANUAL],
                    more=int,namespace=self,doc=docstr);
    # manual selection menu only shown in manual mode
    def show_chansel_menu (value):
      chan_menu.show(value == CHANMODE_MANUAL);
    chan_opt.when_changed(show_chansel_menu);
    # 
    self._opts += [
      chan_opt,
      chan_menu,
      TDLOption('imaging_weight',"Imaging weights",
                ["default","natural","uniform","briggs","radial"],namespace=self),
      TDLMenu("Apply Gaussian taper to visibilities",toggle='imaging_taper_gauss',namespace=self,
            *( TDLOption('imaging_taper_bmaj',"Major axis (arcsec)",[12],more=float,namespace=self),
               TDLOption('imaging_taper_bmin',"Minor axis (arcsec)",[12],more=float,namespace=self),
               TDLOption('imaging_taper_bpa',"Position angle (deg)",[0],more=float,namespace=self),
            )),
      TDLOption('imaging_stokes',"Stokes parameters to image",
                ["I","IQUV"],namespace=self)
    ];
    if npix:
      if not isinstance(npix,(list,tuple)):
        npix = [ npix ];
      self._opts.append(TDLOption('imaging_npix',
                        "Image size, in pixels",npix,more=int,namespace=self));
    self.imaging_arcmin = self.imaging_cellsize = None;
    if arcmin or not cellsize:
      if cellsize:
        raise ValueError,"ImagingSelector: specify either 'arcmin' or 'cellsize', not both";
      if not arcmin:
        arcmin = [ 5 ];
      elif not isinstance(arcmin,(list,tuple)):
        arcmin = [ arcmin ];
      self._opts.append(TDLOption('imaging_arcmin',
                    "Image size, in arcmin",list(arcmin)+["all-sky"],more=float,namespace=self));
    elif cellsize:
      if not isinstance(cellsize,(list,tuple)):
        cellsize = [ cellsize ];
      self._opts.append(TDLOption('imaging_cellsize',
                        "Pixel size",cellsize,more=str,namespace=self));
    self._opts.append(TDLOption('imaging_padding',
                                "Image padding factor for FFTs",[1.0],more=float,namespace=self,
      doc="""When gridding and transforming, the array may be padded
      by this factor in the image plane. This reduces aliasing, especially
      in wide-field cleaning."""));
    # add w-projection option
    self._opts.append(TDLMenu("Enable w-projection",
        TDLOption('imaging_wprojplanes',"Number of convolution functions for w-projection",
                  [128],more=int,namespace=self),
      toggle='imaging_enable_wproj',namespace=self,
      doc="""This option enables the w-projection algorithm. Imaging will
      be <b>much</b> slower, but wide-field artefacts will be greately reduced or
      eliminated. Note that this option does not work right with all versions
      of the aips++ imager, but the only way to find out is to give it a try.""",
      ));
    # add center
    self._opts.append(TDLOption('imaging_phasecenter',"Phase center",["default"],namespace=self,more=str,
        doc="""You can center the image on a particular point in the sky. The default is
        the phase center of the observation. To override this, enter a direction string 
        of the form, e.g., 'J2000,05h35m10s,-30deg15m30s'"""));
    # add MS subset selector, if needed
    if subset:
      self.subset_selector = mssel.make_subset_selector(namespace);
      custom_sel_menu = TDLMenu("Use custom MS selection for imaging",
                          toggle="imaging_custom_ms_select",namespace=self,
                          *self.subset_selector.option_list());
      self._opts.append(custom_sel_menu);
    else:
      self.subset_selector = mssel.subset_selector;
    # add TDL job to make an image
    def job_make_image (mqs,parent,**kw):
      self.make_dirty_image();
    self._opts.append(TDLJob(job_make_image,"Make a dirty image"));

  def option_list (self):
    """Returns list of all TDL options"""
    return self._opts;

  def make_dirty_image (self,npix=None,cellsize=None,arcmin=None):
    """Runs glish script to make an image.
    The following parameters, if supplied, will override the option settings:
      npix:       image size in pixels
      arcmin:     image size in arc minutes
      cellsize:   default cell size (as an aips++ quantity string). Note
                  that either arcmin or cellsize may be specified, not both.
    """;
    # choose selector based on custom MS select option
    if self.imaging_custom_ms_select:
      selector = self.subset_selector;
    else:
      selector = self.mssel.subset_selector;
    # check MS name and column
    if not self.mssel.msname:
      raise ValueError,"make_dirty_image: MS not set up";
    col = self.imaging_column;
    if not col:
      raise ValueError,"make_dirty_image: output column not set up";
    # image size
    npix = npix or self.imaging_npix;
    # see if arcmin or cellsize is supplied; if neither, then use TDLOption values
    if arcmin or cellsize:
      if arcmin and cellsize:
          raise ValueError,"ImagingSelector: specify either 'arcmin' or 'cellsize', not both";
    else:
      arcmin = self.imaging_arcmin;
      cellsize = self.imaging_cellsize;
    # resolve to required cellsize, finally
    if arcmin is not None:
      if arcmin == 'all-sky':
       # for all-sky images, fudge: 120" cell gives full sky in ~3200 pixels
        cellsize = str(120*(3200./npix))+"arcsec";
      else:
        cellsize = str(float(arcmin*60)/npix)+"arcsec";
    if self.imaging_chanmode == CHANMODE_MFS:
      imgmode = "mfs";
    else:
      imgmode = "channel";
    # form up initial argument list to run imaging script
    if self.imager_type == "lwimager":
      _IMAGER = "python";
      script_name = os.path.join(Meow._meow_path,'make_dirty_image.py');
      args = [ 'python',script_name,'data='+col ];
      offset = 0;
    else:
      _IMAGER = "glish";
      script_name = os.path.join(Meow._meow_path,'make_dirty_image.g');
      script_name = os.path.realpath(script_name);  # glish don't like symlinks...
      args = [ 'glish','-l',
        script_name,
        col ];
      offset = 1;
    # these arguments are common to both imagers
    args += \
      [ 'ms='+self.mssel.msname,
        'mode='+imgmode,
        'weight='+self.imaging_weight,
        'stokes='+self.imaging_stokes,
        'npix=%d'%npix,
        'cellsize='+cellsize,
        'spwid=%d'%(selector.get_spectral_window()+offset),
        'field=%d'%(selector.get_field()+offset),
        'padding=%f'%self.imaging_padding,
      ];
    # add taper arguments
    if self.imaging_taper_gauss:
      if _IMAGER == "glish":
        args += [ 
          'filter_bmaj=%farcsec'%self.imaging_taper_bmaj,
          'filter_bmin=%farcsec'%self.imaging_taper_bmin,
          'filter_bpa=%fdeg'%self.imaging_taper_bpa
        ];
      else:
        args.append("filter=%farcsec,%farcsec,%fdeg"%
            (self.imaging_taper_bmaj,self.imaging_taper_bmin,self.imaging_taper_bpa));
    # add w-proj arguments
    if self.imaging_enable_wproj:
      args.append("wprojplanes=%d"%self.imaging_wprojplanes);
    # add phase center
    if self.imaging_phasecenter and self.imaging_phasecenter != "default":
      args.append("phasecenter=%s"%self.imaging_phasecenter.replace(" ",""));
    # add channel arguments for setdata
    chans = selector.get_channels();
    totchan = selector.get_total_channels() or self.imaging_totchan;
    if chans:
      nchan = chans[1]-chans[0]+1;
      chanstart = chans[0]+offset;
      if len(chans) > 2:
        chanstep = chans[2];
        nchan /= chanstep;
      else:
        chanstep = 1;
      args += [ 'chanmode=channel',
                'nchan='+str(nchan),
                'chanstart='+str(chanstart),
                'chanstep='+str(chanstep) ];
    else:
      args.append("chanmode=none");
      chanstart,nchan,chanstep = offset,totchan,1;
    # add channel arguments for setimage
    if self.imaging_chanmode == CHANMODE_MANUAL:
      img_nchan     = self.imaging_nchan;
      img_chanstart = self.imaging_chanstart+offset;
      img_chanstep  = self.imaging_chanstep;
    elif self.imaging_chanmode == CHANMODE_NOAVG:
      img_nchan     = nchan;
      img_chanstart = chanstart;
      img_chanstep  = chanstep;
    elif self.imaging_chanmode == CHANMODE_ALL:
      img_nchan     = 1;
      img_chanstart = chanstart;
      img_chanstep  = nchan;
    elif self.imaging_chanmode == CHANMODE_MFS:
      img_nchan     = 1;
      img_chanstart = chanstart;
      img_chanstep  = nchan;
    elif isinstance(self.imaging_chanmode,int):
      img_nchan     = self.imaging_chanmode;
      img_chanstart = chanstart;
      img_chanstep  = nchan/self.imaging_chanmode;
    args += [ 'img_nchan='+str(img_nchan),
              'img_chanstart='+str(img_chanstart),
              'img_chanstep='+str(img_chanstep) ];
    # add TaQL string
    taql = selector.get_taql_string();
    if taql:
      args.append("select=%s"%taql);
    # figure out an output FITS filename
    fitsname = self.mssel.msname;
    if fitsname.endswith('/'):
      fitsname = fitsname[:-1];
    fitsname = os.path.basename(fitsname);
    fitsname = "%s.%s.%s.%dch.fits"%(fitsname,col,imgmode,img_nchan);
    args += [ 'fits='+fitsname ];
    # if the fits file exists, clobber it
    if os.path.exists(fitsname):
      try:
        os.unlink(fitsname);
      except:
        pass; 
    print "MSUtils: imager args are",args;
    # run script
    os.spawnvp(os.P_NOWAIT,_IMAGER,args);

# keep a global map of flagsets associated with each MS,
# so that different "customers" can always deal with the same flagset object
_flagset_map = {};
def get_flagsets (ms):
  global _flagset_map;
  msname = ms.name();
  fs = _flagset_map.get(msname,None);
  if fs is None:
    fs = _flagset_map[msname] = Flagsets(ms);
  return fs;

class Flagsets (object):
  """Flagsets implements a class to manage an MS's flagsets""";
  def __init__ (self,ms):
    """Flagsets is constructed from an MS.""";
    self.msname = ms.name();
    
  def load (self,ms):
    self._wc_callbacks = [];
    if not 'BITFLAG' in ms.colnames():
      self.order = None;
      self.bits = {};
    else:
      kws = ms.colkeywordnames('BITFLAG');
      self.bits  = {};
      # scan FLAGSET_xxx keywords and populate name->bitmask mappings
      for kw in kws:
        match = re.match('^FLAGSET_(.*)$',kw);
        if match:
          name = match.group(1);
          bit = ms.getcolkeyword('BITFLAG',kw);
          if isinstance(bit,int):
            self.bits[name] = bit;
          else:
            print "Warning: unexpected type (%s) for %s keyword of BITFLAG column, ignoring"%(type(order),kw);
      # have we found any FLAGSET_ specs?
      if self.bits:
        order = 'FLAGSETS' in kws and ms.getcolkeyword('BITFLAG','FLAGSETS');
        if isinstance(order,str):
          order = order.split(',');
        else:
          print "Warning: unexpected type (%s) for FLAGSETS keyword of BITFLAG column, ignoring"%type(order);
          order = [];
        # form up "natural" order by comparing bitmasks
        bitwise_order = list(self.bits.iterkeys());
        bitwise_order.sort(lambda a,b:cmp(self.bits[a],self.bits[b]));
        # if an order is specified, make sure it is actually valid,
        # and add any elements from bitwise_order that are not present
        self.order = [ fs for fs in order if fs in self.bits ] + \
                     [ fs for fs in bitwise_order if fs not in order ];
        # if order was fixed compared to what was in MS, write back to MS
        if ms.iswritable() and self.order != order:
          ms._putkeyword('BITFLAG','FLAGSETS',-1,False,','.join(self.order));
          ms.flush();
      # else if no flagsets found, try the old-style NAMES keyword
      elif 'NAMES' in kws:
        names = ms.getcolkeyword('BITFLAG','NAMES');
        if isinstance(names,(list,tuple)):
          self.order = map(str,names);
          bit = 1;
          for name in self.order:
            self.bits[name] = bit;
            bit <<= 1;
          if ms.iswritable():
            ms._putkeyword('BITFLAG','FLAGSETS',-1,False,','.join(self.order));
            for name,bit in self.bits.iteritems():
              ms._putkeyword('BITFLAG','FLAGSET_%s'%name,-1,False,bit);
            ms.flush();
      else:
        self.order = [];
        
  def when_changed (self,callback):
    self._wc_callbacks.append(callback);
    callback();
  
  def names (self):
    """Returns a list of flagset names, in the order in which they were
    created. Returns None if BITFLAG column is missing (so flagsets are
    not available.)""";
    return self.order;
  
  def flagmask (self,name,create=False):
    """Returns flagmask corresponding to named flagset. If flagset does not exist:
      * if create is True, creates a new one
      * if create is False, raises exception
    """;
    # lookup flagbit, return if found
    if self.bits is None:
      raise TypeError,"MS does not contain a BITFLAG column, cannot use flagsets""";
    bit = self.bits.get(name,None);
    if bit is not None:
      return bit;
    # raise exception if not allowed to create a new one
    if not create:
      raise ValueError,"Flagset '%s' not found"%name;
    # find empty bit
    for bitnum in range(32):
      bit = 1<<bitnum;
      if bit not in self.bits.values():
        self.order.append(name);
        self.bits[name] = bit;
        ms = TABLE(self.msname,readonly=False);
        ms._putkeyword('BITFLAG','FLAGSETS',-1,False,','.join(self.order));
        ms._putkeyword('BITFLAG','FLAGSET_%s'%name,-1,False,bit);
        ms.flush();
        for callback in self._wc_callbacks:
          callback();
        return bit;
    # no free bit found, bummer
    raise ValueError,"Too many flagsets in MS, cannot create another one";
    
  def remove_flagset (self,*fsnames):
    """Removes the named flagset(s). Returns flagmask corresponding to the removed
    flagsets.""";
    # lookup all flagsets, raise error if any not found
    if self.bits is None:
      raise TypeError,"MS does not contain a BITFLAG column, cannot use flagsets""";
    removing = [];
    for fs in fsnames:
      bit = self.bits.get(fs,None);
      if bit is None:
        raise ValueError,"Flagset '%s' not found"%fs;
      removing.append((fs,bit));
    if not removing:
      return;
    # remove items, form up mask of bitflags to be cleared
    ms = TABLE(self.msname,readonly=False);
    mask = 0;
    for name,bit in removing:
      mask |= bit;
      del self.bits[name];
      del self.order[self.order.index(name)];
      ms.removecolkeyword('BITFLAG','FLAGSET_%s'%name);
    # write new list of bitflags
    ms._putkeyword('BITFLAG','FLAGSETS',-1,False,','.join(self.order));
    ms.flush();
    for callback in self._wc_callbacks:
      callback();
    return mask;
    
def parse_antenna_spec (spec,names=None):
  """Parses string containing antenna number or name, returns number""";
  match = re.match("^(\d+)|([^:,\s]+)$",spec);
  if match:
    number,name = match.group(1,2);
    if number is not None:
      number = int(number);
      if number >=0 and (not names or number < len(names)):
        return number;
    elif names and name in names:
      return names.index(name);
  raise ValueError,"illegal antenna specifier '%s'"%spec;

def parse_antenna_subset (value,names=None):
  if not value:
    return None;
  if names:
    nant = len(names);
  else:
    nant = 1;     # unreasonably large...
  subset = [];
  for spec in re.split("[\s,]+",value):
    # single element 
    if spec.find(':') < 0:
      index = parse_antenna_spec(spec);
      if not names:
        nant = max(nant,index+1);
      subset.append(index);
      continue;
    # range from:to
    else:
      spec1,spec2 = spec.split(':',1);
      index1 = 0;
      index2 = nant-1;
      if spec1:
        index1 = parse_antenna_spec(spec1);
      if spec2:
        index2 = parse_antenna_spec(spec2);
      if not names:
        nant = max(nant,index1+1,index2+1);
      if index1>index2:
        raise ValueError,"illegal antenna specifier '%s'"%spec;
      # add to subset
      subset += range(index1,index2+1);
  return subset;
