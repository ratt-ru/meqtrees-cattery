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
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from Timba.TDL import *
from Timba.Meq import meq
import Timba.Apps
import Timba.array
import math
import random
import Meow

try:
  import pycasatable
except:
  pycasatable = None;
  
msname = '';
input_column = output_column = imaging_column = None;
tile_size = None;
ddid_index = None;
field_index = None;
ms_channels = None;
ms_channel_start = None;
ms_channel_end = None;
ms_channel_step = None;
ms_data_columns = None;
ms_antenna_names = [];
ms_field_names = None;
ms_ddid_numchannels = None;
ms_write_flags = False;
ms_input_flag_bit = 1;
ms_has_output = False;

def include_ms_options (
    has_input=True,
    has_output=True,
    tile_sizes=[1,5,10,20,30,60],
    ddid=[0,1],
    channels=None,
    flags=False
  ):
  """Instantiates MS input/output options""";
  TDLRuntimeOptions(*ms_options(has_input,has_output,tile_sizes,ddid,channels,flags));

def ms_options (
    has_input=True,
    has_output=True,
    tile_sizes=[1,5,10,20,30,60],
    ddid=[0],
    field_id=None,
    channels=None,
    flags=False
  ):
  """Returns list of MS input/output options""";
  global ms_option,input_col_option,ddid_option,field_option,channel_options;
  ms_option = TDLOption('msname',"MS",TDLDirSelect("*.ms *.MS"));
  opts = [ ms_option ];
  if has_input:
    input_col_option = TDLOption('input_column',"Input MS column",
                        ["DATA","MODEL_DATA","CORRECTED_DATA"],default=0);
    opts.append(input_col_option);
  else:
    input_col_option = None;
  if has_output:
    global ms_has_output;
    ms_has_output = True;
    opts.append(TDLOption('output_column',"Output MS column",
                        ["DATA","MODEL_DATA","CORRECTED_DATA",None],default=2));
  if tile_sizes:
    opts.append(TDLOption('tile_size',"Tile size (timeslots)",tile_sizes,more=int));
  more_ids = int;
  if pycasatable:
    ddid = [0];
    field_id = [0];
    more_ids = None;  # no more= option to ddid/field: they come from the MS
  if ddid:
    ddid_option = TDLOption('ddid_index',"Data description ID",ddid,more=more_ids,
      doc="""If the MS contains multiple spectral windows, etc., then use this option to select different DATA_DESCRIPTION_IDs. Default is 0.""");
    opts.append(ddid_option);
  if field_id:
    field_option = TDLOption('field_index',"Field ID",field_id,more=more_ids,
      doc="""If the MS contains multiple fields, then use this option to select different FIELD_IDs. Default is 0.""");
    opts.append(field_option);
  channel_options = [];
  if channels:
    if pycasatable:
      channel_options = [
        TDLOption('ms_channel_start',"First channel",[0],more=int),
        TDLOption('ms_channel_end',"Last channel",[0],more=int),
        TDLOption('ms_channel_step',"Channel stepping",[1,2],more=int)
      ];
      opts.append(TDLMenu("Channel selection",*channel_options));
    else:
      if not isinstance(channels,(list,tuple)):
        channels = [[0,0,1]];
      opts.append(TDLOption('ms_channels',"Channel selection",channels));
  if flags:
    opts.append(TDLOption('ms_write_flags',"Write flags to output",False));
  # if pycasatable exists, set up interactivity for MS options
  if pycasatable:
    global _select_new_ms;
    ms_option.when_changed(_select_new_ms);
    ddid_option.when_changed(_select_ddid);
    if channel_options:
      channel_options[0].set_validator(_validate_first_channel);
      channel_options[1].set_validator(_validate_last_channel);
      channel_options[2].set_validator(_validate_channel_step);
  return opts;
  
def _select_new_ms (msname):
  if not msname:
    return;
  try:
    ms = pycasatable.table(msname);
    # data columns
    global ms_data_columns;
    ms_data_columns = [ name for name in ms.colnames() if name.endswith('DATA') ];
    if input_col_option:
      input_col_option.set_option_list(ms_data_columns);
    # antennas
    global ms_antenna_names;
    ms_antenna_names = pycasatable.table(ms.getkeyword('ANTENNA')).getcol('NAME');
    # DDIDs
    ddid_table = pycasatable.table(ms.getkeyword('DATA_DESCRIPTION'));
    spws = ddid_table.getcol('SPECTRAL_WINDOW_ID');
    numchans = pycasatable.table(ms.getkeyword('SPECTRAL_WINDOW')).getcol('NUM_CHAN');
    global ms_ddid_numchannels;
    ms_ddid_numchannels = [ numchans[spw] for spw in spws ];
    ddid_option.set_option_list(list(range(len(ms_ddid_numchannels))));
    if len(ms_ddid_numchannels) < 2:
      ddid_option.hide();
    # Fields
    global ms_field_names;
    ms_field_names = pycasatable.table(ms.getkeyword('FIELD')).getcol('NAME');
    field_option.set_option_list(dict(enumerate(ms_field_names)));
    if len(ms_field_names) < 2:
      field_option.hide();
  except:
    print("error reading MS",msname);
    traceback.print_exc();
    return False;
  
def _select_ddid (value):
  if ms_ddid_numchannels:
    nchan = ms_ddid_numchannels[value];
    if channel_options:
      channel_options[0].set_option_list([0,nchan-1]);
      channel_options[0].set_value(0);
      channel_options[0].set_doc("(max %d) select first channel"%(nchan-1));
      channel_options[1].set_option_list([0,nchan-1]);
      channel_options[1].set_value(nchan-1);
      channel_options[1].set_doc("(max %d) select last channel"%(nchan-1));
      channel_options[2].set_value(1);
      channel_options[2].set_doc("(max %d) select channel steping"%nchan);

def _validate_first_channel (value):
  return isinstance(value,int) and \
         value >=0 and \
         value < channel_options[1].value;

def _validate_last_channel (value):
  return isinstance(value,int) and \
         value >= channel_options[0].value and \
         value < ms_ddid_numchannels[ddid_option.value];

def _validate_channel_step (value):
  return isinstance(value,int) and \
         value >= 1 and \
         value <= ms_ddid_numchannels[ddid_option.value];
  
imaging_npix = 256;
imaging_cellsize = '1arcsec';
imaging_arcmin = None;
imaging_channels = [32,1,1];
imaging_channels_specified = False;

def include_imaging_options (npix=None,arcmin=5,cellsize=None,channels=None):
  """Instantiates imager options""";
  TDLRuntimeOptions(*imaging_options(npix,arcmin,cellsize,channels));

def imaging_options (npix=None,arcmin=5,cellsize=None,channels=None):
  """Instantiates imager options""";
  opts = [
    TDLOption('imaging_mode',"Imaging mode",["mfs","channel"]),
    TDLOption('imaging_weight',"Imaging weights",["natural","uniform","briggs"]),
    TDLOption('imaging_stokes',"Stokes parameters to image",["I","IQUV"]) 
  ];
  if npix:
    if not isinstance(npix,(list,tuple)):
      npix = [ npix ];
    opts.append(TDLOption('imaging_npix',"Image size, in pixels",npix,more=int));
  if arcmin:
    if cellsize:
      raise ValueError("include_imaging_options: specify cellsize or arcmin, not both");
    if not isinstance(arcmin,(list,tuple)):
      arcmin = [ arcmin ];
    opts.append(TDLOption('imaging_arcmin',"Image size, in arcmin",arcmin,more=float));
  elif cellsize:
    if not isinstance(cellsize,(list,tuple)):
      cellsize = [ cellsize ];
    opts.append(TDLOption('imaging_cellsize',"Pixel size",cellsize,more=str));
  if channels:
    opts.append(TDLOption('imaging_channels',"Imaging channels selection",channels));
    imaging_channels_specified = True;
  def job_make_image (mqs,parent,**kw):
    make_dirty_image();
  if ms_has_output:
    opts.append(TDLJob(job_make_image,"Make image from MS output column"));
  else:
    opts.append(TDLOption('imaging_column',"MS column to image",["DATA","MODEL_DATA","CORRECTED_DATA"]));
    opts.append(TDLJob(job_make_image,"Make image"));
  return opts;  
  

source_table = "sources.mep";
mep_table = "calib.mep";

def get_source_table ():
  return (msname or '.') + "/" + source_table;

def get_mep_table ():
  return (msname or '.') + "/" + mep_table;

_solver_opts = dict(
  debug_level  = 0,
  colin_factor = 1e-6,
  lm_factor    = .001,
  balanced_equations = False,
  epsilon      = 1e-4,
  num_iter     = 30,
  convergence_quota = 0.8
);

_solver_opts = {};

def solver_options (optionset='_solver_opts',namespace=None):
  """Returns list of solver options.
  Default places options into dict at Meow.Utils._solver_opts. To make another set of options,
  supply a different optionset.""";
  optionset += '.';
  return [
    TDLOption(optionset+'debug_level',"Solver debug level",[0,1,10],namespace=namespace),
    TDLOption(optionset+'colin_factor',"Collinearity factor",[1e-8,1e-6,1e-3,1e-1],default=1,more=float,namespace=namespace),
    TDLOption(optionset+'lm_factor',"Initial LM factor",[1,.1,.01,.001],default=3,more=float,namespace=namespace),
    TDLOption(optionset+'balanced_equations',"Assume balanced equations",False,namespace=namespace),
    TDLOption(optionset+'epsilon',"Convergence threshold",[.01,.001,.0001,1e-5,1e-6],default=2,more=float,namespace=namespace),
    TDLOption(optionset+'num_iter',"Max iterations",[10,30,50,100,1000],default=1,more=int,namespace=namespace),
    TDLOption(optionset+'convergence_quota',"Subtiling convergence quota",[.8,.9,1.],namespace=namespace) \
  ];

def parameter_options ():
  return [
    TDLOption('use_previous',"Reuse solution from previous time interval",True,
    doc="""If True, solutions for successive time domains will start with
  the solution for a previous domain. Normally this speeds up convergence; you
  may turn it off to re-test convergence at each domain."""),
    TDLOption('use_mep',"Reuse solutions from MEP table",True,
    doc="""If True, solutions from the MEP table (presumably, from a previous
  run) will be used as starting points. Turn this off to solve from scratch.""")
 ];

def create_solver_defaults (solvables,options=None):
  global _solver_opts;
  opts = dict(_solver_opts);
  if options:
    opts.update(options);
  print(opts);
  # copy all options into solver defaults with the same name
  solver_defaults = record(**opts);
  # additionally, set epsilon_deriv
  solver_defaults.epsilon_deriv      = opts["epsilon"];
  solver_defaults.save_funklets    = True
  solver_defaults.last_update      = True
  solver_defaults.solvable         = record(command_by_list=(record(name=solvables,
                                            state=record(solvable=True)),
                                            record(state=record(solvable=False))))
  return solver_defaults

def set_node_state (mqs,node,fields_record):
  """helper function to set the state of a node specified by name or
  nodeindex""";
  rec = record(state=fields_record);
  if isinstance(node,str):
    rec.name = node;
  elif isinstance(node,int):
    rec.nodeindex = node;
  else:
    raise TypeError('illegal node argument');
  # pass command to kernel
  mqs.meq('Node.Set.State',rec);
  pass

# various global MS I/O options
ms_output = True;
ms_queue_size = 500;
ms_selection = None;


def create_inputrec (tiling=None):
  global tile_size;
  global ddid_index;
  if msname is None:
    raise ValueError("MS not specified in options");
  boioname = "boio."+msname+".predict."+str(tiling or tile_size);
  # if boio dump for this tiling exists, use it to save time
  if not ms_selection and os.access(boioname,os.R_OK):
    rec = record(boio=record(boio_file_name=boioname,boio_file_mode="r"));
  # else use MS, but tell the event channel to record itself to boio file
  else:
    rec = record();
    rec.ms_name          = msname
    if input_column:
      rec.data_column_name = input_column;
    # use global tile_size setting if tiling not specified
    if not tiling:
      tiling = tile_size;
    if isinstance(tiling,(list,tuple)):
      if len(tiling) != 2:
        raise TypeError("tiling: 2-list or 2-tuple expected");
      (tile_segments,tile_size) = tiling;
      if tile_segments is not None:
        rec.tile_segments    = tile_segments;
      if tile_size is not None:
        rec.tile_size        = tile_size;
    else:  
      rec.tile_size = tiling;
    rec.selection = ms_selection or record();
    if ms_channels is not None:
      rec.selection.channel_start_index = ms_channels[0];
      rec.selection.channel_end_index = ms_channels[1];
      if len(ms_channels) > 2:
        rec.selection.channel_increment = ms_channels[2];
    elif channel_options is not None:
      rec.selection.channel_start_index = ms_channel_start;
      rec.selection.channel_end_index = ms_channel_end;
      rec.selection.channel_increment = ms_channel_step;
    if ddid_index is not None:
      rec.selection.ddid_index = ddid_index;
    if field_index is not None:
      rec.selection.field_index = field_index;
    rec = record(ms=rec);
  rec.python_init = 'Meow.ReadVisHeader';
  rec.mt_queue_size = ms_queue_size;
  return rec;

def create_outputrec ():
  rec = record();
  rec.mt_queue_size = ms_queue_size;
  if ms_output:
    rec.write_flags    = ms_write_flags;
    if output_column:
      rec.data_column = output_column;
    return record(ms=rec);
  else:
    rec.boio_file_name = "boio."+msname+".solve."+str(tile_size);
    rec.boio_file_mode = 'W';
    return record(boio=rec);
    
def create_io_request (tiling=None):
  req = meq.request();
  req.input  = create_inputrec(tiling);
  if ms_write_flags or output_column is not None:
    req.output = create_outputrec();
  return req;
  
def phase_parm (tdeg,fdeg):
  """helper function to create a t/f parm for phase, including constraints.
  Placeholder until Maaijke implements periodic constraints.
  """;
  polc = meq.polc(Timba.array.zeros((tdeg+1,fdeg+1))*0.0,
            scale=array([3600.,8e+8,0,0,0,0,0,0]));
  shape = [tdeg+1,fdeg+1];
  # work out constraints on coefficients
  # maximum excursion in freq is pi/2
  # max excursion in time is pi/2
  dt = .2;
  df = .5;
  cmin = [];
  cmax = [];
  for it in range(tdeg+1):
    for jf in range(fdeg+1):
      mm = math.pi/(dt**it * df**jf );
      cmin.append(-mm);
      cmax.append(mm);
  cmin[0] = -1e+9;
  cmax[0] = 1e+9;
  return Meq.Parm(polc,shape=shape,real_polc=polc,node_groups='Parm',
                  constrain_min=cmin,constrain_max=cmax,
                  table_name=get_mep_table());


def perturb_parameters (mqs,solvables,pert="random",
                        absolute=False,random_range=[0.2,0.3],constrain=None):
  for name in solvables:
    polc = mqs.getnodestate(name).real_polc;
    if absolute:  # absolute pert value given
      polc.coeff[0,0] += pert;
    elif pert == "random":  # else random pert
      polc.coeff[0,0] *= 1 + random.uniform(*random_range)*random.choice([-1,1]);
    else: # else perturb in relative terms
      polc.coeff[0,0] *= (1 + pert);
    parmstate = record(init_funklet=polc,
      use_previous=use_previous,reset_funklet=not use_mep);
    if constrain is not None:
      parmstate.constrain = constrain;
    set_node_state(mqs,name,parmstate);
  return solvables;
    
def reset_parameters (mqs,solvables,value=None,use_table=False,reset=False):
  for name in solvables:
    polc = mqs.getnodestate(name).real_polc;
    if value is not None:
      polc.coeff[()] = value;
    reset_funklet = reset or not (use_table or use_mep);
    set_node_state(mqs,name,record(init_funklet=polc,
      use_previous=use_previous,reset_funklet=reset_funklet));
  return solvables;

def run_solve_job (mqs,solvables,tiling=None,solver_node="solver",vdm_node="VisDataMux",
                   options=None):
  """common helper method to run a solution with a bunch of solvables""";
  # set solvables list in solver
  solver_defaults = create_solver_defaults(solvables,options=options)
  set_node_state(mqs,solver_node,solver_defaults)

  req = create_io_request(tiling);
  
  mqs.execute(vdm_node,req,wait=False);
  pass
  
def make_dirty_image (npix=None,cellsize=None,arcmin=None,channels=None,**kw):
  """Runs glish script to make an image
  npix is image size, in pixels
  cellsize is pixel size, as an aips++ Measures string (e.g. "0.5arcsec")
  arcmin is image size, in arcminutes. Either cellsize or arcmin must be 
        specified, but not both.
  """;
  col = output_column or imaging_column;
  if not col:
    raise ValueError("make_dirty_image: output column not set up");
  if not msname:
    raise ValueError("make_dirty_image: MS not set up");
  npix = (npix or imaging_npix);
  arcmin = (arcmin or imaging_arcmin);
  if arcmin is not None:
    cellsize = str(float(arcmin*60)/npix)+"arcsec";
  import os
  import os.path
  # if explicit channels are specified, use them 
  if channels:
    nchan,chanstart,chanstep = channels;
  # else if MS channels were specified use them
  elif ms_channels and not imaging_channels_specified:
    nchan = ms_channels[1]-ms_channels[0]+1;
    chanstart = ms_channels[0]+1;
    if len(ms_channels) > 2:
      chanstep = ms_channels[2];
    else:
      chanstep = 1;
  # else use the imaging_channels option
  else:
    nchan,chanstart,chanstep = imaging_channels;
  script_name = os.path.join(Meow._meow_path,'make_dirty_image.g');
  script_name = os.path.realpath(script_name);  # glish don't like symlinks...
  args = [ 'glish','-l',
    script_name,
    col,
    'ms='+msname,'mode='+imaging_mode,
    'weight='+imaging_weight,'stokes='+imaging_stokes,
    'npix='+str(npix),
    'cellsize='+(cellsize or imaging_cellsize),
    'nchan='+str(nchan),
    'chanstart='+str(chanstart),    
    'chanstep='+str(chanstep)
  ];
  print(args);
  Timba.Apps.spawnvp_nowait('glish',args);
  

class ListOptionParser (object):
  def __init__ (self,minval=None,maxval=None,name=""):
    self.minval,self.maxval = minval,maxval;
    self.name = name or "index";
    
  def set_min (self,minval):
    self.minval = minval;
  def set_max (self,maxval):
    self.maxval = maxval;
    
  def _validate_index (self,index,minval,maxval):
    if self.minval is None:
      minval = min(minval,index);
    if self.maxval is None:
      maxval = max(maxval,index);
    if index < minval or index > maxval:
      raise ValueError("illegal %s '%d'"%(self.name,index));
    return index,minval,maxval;
    
  def parse_list (self,value):
    """Parses string list of values, returns list of numbers""";
    if not value:
      return None;
    if self.minval is not None:
      minval = self.minval;
    else:
      minval = 0;
    if self.maxval is not None:
      maxval = self.maxval;
    else:
      maxval = 0;
    subset = [];
    for spec in re.split("[\s,]+",value):
      if spec:
        # single number
        match = re.match("^\d+$",spec);
        if match:
          index,minval,maxval = self._validate_index(int(spec),minval,maxval);
          subset.append(index);
          continue;
        # [number]:[number]
        match = re.match("^(\d+)?:(\d+)?$",spec);
        if not match:
          raise ValueError("illegal %s '%s'"%(self.name,spec));
        if match.group(1):
          index1 = int(match.group(1));
        elif self.minval is not None:
          index1 = self.minval;
        else:
          raise ValueError("illegal %s '%s'"%(self.name,spec));
        if match.group(2):
          index2 = int(match.group(2));
        elif self.maxval is not None:
          index2 = self.maxval;
        else:
          raise ValueError("illegal %s '%s'"%(self.name,spec));
        index1,minval,maxval = self._validate_index(index1,minval,maxval);
        index2,minval,maxval = self._validate_index(index2,minval,maxval);
        # add to subset
        subset += list(range(index1,index2+1));
    return subset;

  def validator (self,value):
    self.parse_list(value);
    return True;
