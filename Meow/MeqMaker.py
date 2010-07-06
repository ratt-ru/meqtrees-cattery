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
import inspect
import types

import Meow
from Meow import StdTrees,ParmGroup,Parallelization,MSUtils

_annotation_label_doc = """The following directives may be embedded in the format string:
  '%N':       source name
  '%#':       source ordinal number
  '%Rr':      source distance from phase centre (if position is fixed), in radians
  '%Rd':        same, in degrees
  '%Rm':        same, in arcmin
  '%Rs':        same, in arcsec
  '%I','%Q','%U','%V':
              source fluxes
  '%attr':    value of source attribute 'attr', empty string if no such attribute
  '%(attr)format':  value of source attribute 'attr' ('attr' can also be one of N,#,Rr,Rd, etc.),
              converted using the specified %format a-la printf().
  '%%': '%' character""";

def _modname (obj):
  if hasattr(obj,'name'):
    name = obj.name;
  elif hasattr(obj,'__name__'):
    name = obj.__name__;
  else:
    name = obj.__class__.__name__;
  return name;

def _modopts (mod,opttype='compile'):
  """for the given module, returns a list of compile/runtime options suitable for passing
  to TDLMenu. If the module implements a compile/runtime_options() method, uses that,
  else steals options from the module itself.""";
  modopts = getattr(mod,opttype+'_options',None);
  # if module defines an xx_options() callable, use that
  if callable(modopts):
    return list(modopts());
  # else if this is a true module, it may have options to be stolen, so insert as is
  elif inspect.ismodule(mod):
    return TDLStealOptions(mod,is_runtime=(opttype!='compile'));
  # else item is an object emulating a module, so insert nothing
  else:
    return [];

SKYJONES_LM       = "lm plane";
SKYJONES_RADEC    = "RA-Dec full sky";
SKYJONES_AZEL     = "Az-El half sky";
SKYJONES_AZEL_FULL = "Az-El full sky";

class MeqMaker (object):
  def __init__ (self,namespace='me',solvable=False,
      use_decomposition=None,
      use_jones_inspectors=None,
      use_skyjones_visualizers=True):
    self.tdloption_namespace = namespace;
    self._uv_jones_list = [];
    self._sky_jones_list = [];
    self._uv_vpm_list = [];
    self._sky_vpm_list = [];
    self._compile_options = [];
    self._runtime_options = None;
    self._runtime_vis_options  = [];
    self._sky_models = None;
    self._source_list = None;
    self.export_kvis = False;
    self._solvable = solvable;
    self._inspectors = [];
    self._skyjones_visualizer_mux = [];
    self._skyjones_visualizer_source = False;
    other_opt = [];
    if use_decomposition is None:
      self.use_decomposition = False;
      self.use_decomposition_opt = \
        TDLOption('use_decomposition',"Use source coherency decomposition, if available",False,namespace=self,
          doc="""If your source models are heavy on point sources, then an alternative form of the M.E. --
          where the source coherency is decomposed into per-station contributions -- may produce
          faster and/or more compact trees. Check this option to enable. Your mileage may vary."""
        );
      other_opt.append(self.use_decomposition_opt);
    else:
      self.use_decomposition = use_decomposition;

    if use_jones_inspectors is None:
      self.use_jones_inspectors = True;
      self.use_jones_inspectors_opt = \
        TDLOption('use_jones_inspectors',"Enable inspectors for Jones modules",True,namespace=self,
          doc="""If enabled, then your trees will automatically include inspector nodes for all
          Jones terms. This will slow things down somewhat -- perhaps a lot, in an MPI configuration --
          so you might want to disable this in production trees."""
        );
      other_opt.append(self.use_jones_inspectors_opt);
    else:
      self.use_jones_inspectors = use_jones_inspectors;

    if use_skyjones_visualizers:
      self.use_skyjones_visualizers = True;
      self.use_skyjones_visualizers_opt = \
        TDLMenu("Include visualizers for sky-Jones terms",
            TDLOption("_skyvis_frame","Coordinate grid",
              [SKYJONES_LM,SKYJONES_RADEC,SKYJONES_AZEL,SKYJONES_AZEL_FULL],namespace=self),
          doc="""If enabled, then your trees will automatically include visualizer nodes for all
          sky-Jones terms. This does not affect normal performance if you do not use the visualizations.
          <B>WARNING:</B> this feature is still experimental, and may not support all sky-Jones modules.
          Turn this off if you experience unexpected compilation errors.
          """,
          toggle='use_skyjones_visualizers',namespace=self
        );
      other_opt.append(self.use_skyjones_visualizers_opt);
    else:
      self.use_skyjones_visualizers = False;

    other_opt.append(
      TDLMenu("Include time & bandwidth smearing",
        TDLOption('smearing_count',"Apply to N brightest sources only",["all",10,100],more=int,namespace=self),
        namespace=self,toggle='use_smearing',default=True,
      )
    );

    other_opt += Meow.IfrArray.compile_options();
    self._compile_options += [ TDLMenu("Measurement Equation options",*other_opt),TDLOptionSeparator("Image-plane components") ];
    # will be True once a uv-plane module has been added
    self._have_uvplane_modules = False;

  def compile_options (self):
    return self._compile_options;

  def add_sky_models (self,modules):
    if self._sky_models:
      raise RuntimeError,"add_sky_models() may only be called once";
    self._compile_options.append(
        self._module_selector("Sky model","sky",modules,use_toggle=False,exclusive=False));
    self._sky_models = modules;
    global _annotation_label_doc;
    self._compile_options.append(
      TDLMenu("Export sky model as kvis annotations",toggle="export_kvis",
                  default=False,namespace=self,
              *(TDLOption("export_karma_filename","Filename",
                          TDLFileSelect("*.ann",exist=False,default=None),namespace=self),
                TDLOption("export_karma_label","Label format",
                          [None,"%N","%N I=%(I).3g","%N Ia=%(Iapp).3g"],more=str,namespace=self,
                          doc="This determines the format of source labels."+_annotation_label_doc),
                TDLOption("export_karma_nlabel","Put labels on brightest N sources only",
                          ["all",10,50],more=int,namespace=self),
                TDLOption("export_karma_incremental","Export each N sources as a separate file",
                          ["all",10,50],more=int,namespace=self)
              ))
    );

  def add_uv_jones (self,label,name=None,modules=None,pointing=None,use_flagger=None):
    return self._add_jones_modules(label,name,modules,is_uvplane=True,pointing=pointing);

  def add_vis_proc_module (self,label,name=None,modules=None):
    return self._add_vpm_modules(label,name,True,modules);

  def add_sky_jones (self,label,name=None,modules=None,pointing=None,use_flagger=None):
    return self._add_jones_modules(label,name,modules,is_uvplane=False,pointing=pointing);

  class VPMTerm (object):
    """Essentially a record representing one VPM in an ME.
    A VPM has the following fields:
      label:      a label (e.g. "G", "E", etc.)
      name:       a descriprive name
      modules:    a list of possible modules implementing the VPM
    """;
    def __init__ (self,label,name,modules):
      self.label      = label;
      self.name       = name;
      self.modules    = modules;
      self.base_node  = None;

  class JonesTerm (object):
    """Essentially a record representing one Jones term in an ME.
    A Jones term has the following fields:
      label:      a label (e.g. "G", "E", etc.)
      name:       a descriptive name
      modules:    a list of possible modules implementing the Jones set
    """;
    def __init__ (self,label,name,modules):
      self.label      = label;
      self.name       = name;
      self.modules    = modules;
      self.base_node  = None;
      self.solvable   = False;

  class SkyJonesTerm (JonesTerm):
    """SkyJonesTerm represents a sky-Jones term.
    It adds a pointing_modules field, holding a list of poiting error modules.
    """;
    def __init__ (self,label,name,modules,pointing_modules=None):
      MeqMaker.JonesTerm.__init__(self,label,name,modules);
      self.pointing_modules = pointing_modules;
      self.subset_selector = None;

  def _add_module_to_compile_options (self,option,is_uvplane):
    if is_uvplane:
      if not self._have_uvplane_modules:
        self._have_uvplane_modules = True;
        self._compile_options.append(TDLOptionSeparator("UV-plane components"));
    else:
      if self._have_uvplane_modules:
        raise RuntimeError,"cannot add sky-plane terms after uv-plane terms";
    self._compile_options.append(option);

  def _add_vpm_modules (self,label,name,is_uvplane,modules):
    if type(label) is types.ModuleType:
      modules = label;
      label = getattr(modules,'__default_label__',None);
      name = getattr(modules,'__default_name__','');
      if label is None:
        raise RuntimeError,"""Module '%s' does not provide a __default_label__ attribute,
              so must be added with an explicit label"""%modules.__name__;
    elif not modules:
      raise RuntimeError,"No modules specified for %s"%name;
    if not isinstance(modules,(list,tuple)):
      modules = [ modules ];
    # Add to internal list
    term = self.VPMTerm(label,name,modules);
    if is_uvplane:
      self._uv_vpm_list.append(term);
    else:
      self._sky_jones_list.append(term);
    # make option menus for selecting a visiblity processor
    mainmenu = self._module_selector("Use %s"%name,label,modules);
    self._add_module_to_compile_options(mainmenu,is_uvplane);

  def _add_jones_modules (self,label,name,modules,is_uvplane=True,pointing=None):
    # see if only a module is specified instead of the full label,name,modules
    # syntax
    if type(label) is types.ModuleType:
      modules = label;
      label = getattr(modules,'__default_label__',None);
      name = getattr(modules,'__default_name__','');
      if label is None:
        raise RuntimeError,"""Module '%s' does not provide a __default_label__ attribute,
              so must be added with an explicit label"""%modules.__name__;
      elif not modules:
        raise RuntimeError,"No modules specified for %s Jones"%label;
    if not isinstance(modules,(list,tuple)):
      modules = [ modules ];
    # extra options for pointing
    extra_options = [];
    if pointing:
      if not isinstance(pointing,(list,tuple)):
        pointing = [ pointing ];
      extra_options += [ self._module_selector("Apply pointing errors to %s"%label,
                                              label+"pe",pointing,nonexclusive=True) ];
    # extra options for per-station Jones term
    sta_opt = TDLOption(self._make_attr('all_stations',label),"Use same %s-Jones for all stations"%label,False,namespace=self,nonexclusive=True,
        doc="""If checked, then the same %s-Jones term will be used for all stations in the array. This may
        make your trees smaller and/or faster. Whether this is a valid approximation or not depends on your
        physics; typically this is valid if your array is small, and the effect arises far from the receivers."""%label);
    if is_uvplane:
      extra_options.append(sta_opt);
      setattr(self,self._make_attr('all_sources',label),False);
      jt = self.JonesTerm(label,name,modules);
    else:
      jt = self.SkyJonesTerm(label,name,modules,pointing);
      tdloption_namespace = "%s.subset_%s"%(self.tdloption_namespace,label);
      jt.subset_selector = SourceSubsetSelector ("Apply %s-Jones to subset of sources"%label,tdloption_namespace,annotate=False);
      opts =  [ sta_opt, TDLOption(self._make_attr('all_sources',label),"Use same %s-Jones for all sources"%label,False,namespace=self,nonexclusive=True,
          doc="""If checked, then the same %s-Jones term will be used for all source in the model. This is a valid
          approximation for narrow fields, and will make your trees smaller and/or faster."""%label) ];
      opts += jt.subset_selector.options;
      extra_options.append(TDLMenu("Advanced options",*opts));
    # make option menus for selecting a jones module
    mainmenu = self._module_selector("Use %s Jones (%s)"%(label,name),
                                     label,modules,extra_opts=extra_options);
    self._add_module_to_compile_options(mainmenu,is_uvplane);
    # Add to internal list
    # each jones module set is represented by a list of:
    #   'label','full name',module_list,pointing_module_list,base_node
    # base_node is inintially None; when a module is ultimately invoked to
    # build the trees, the base node gets stashed here for later reuse
    if is_uvplane:
      self._uv_jones_list.append(jt);
    else:
      self._sky_jones_list.append(jt);

  def runtime_options (self,nest=True):
    if self._runtime_options is None:
      self._runtime_options = [];
      # build list of currently selected modules
      mods = [];
      for module in self._get_selected_modules('sky',self._sky_models):
        mods.append((module,'%s model options'%_modname(module),None,None));
      for jt in self._sky_jones_list:
        mod,name = self._get_selected_module(jt.label,jt.modules), \
                   "%s Jones (%s) options"%(jt.label,jt.name);
        if jt.pointing_modules:
          pemod,pename = self._get_selected_module(jt.label+"pe",jt.pointing_modules), \
                         "Pointing error options";
        else:
          pemod,pename = None,None;
        mods.append((mod,name,pemod,pename));
      for jt in self._uv_jones_list:
        mod,name = self._get_selected_module(jt.label,jt.modules), \
                   "%s Jones (%s) options"%(jt.label,jt.name);
        mods.append((mod,name,None,None));
      for vpm in self._uv_vpm_list:
        mod,name = self._get_selected_module(vpm.label,vpm.modules), \
                   "%s options"%vpm.name;
        mods.append((mod,name,None,None));
      # now go through list and pull in options from each active module
      for mod,name,submod,subname in mods:
        if mod:
          modopts = _modopts(mod,'runtime');
          # add submenu for submodule
          if submod:
            submodopts = _modopts(submod,'runtime');
            if submodopts:
              modopts.append(TDLMenu(subname,*submodopts));
          if nest and modopts:
            self._runtime_options.append(TDLMenu(name,*modopts));
          else:
            self._runtime_options += modopts;
      # add list of visualization options
      if self._runtime_vis_options:
        if self._skyvis_frame == SKYJONES_LM:
          self._runtime_vis_options = [
              TDLOption("_skyvis_xgrid","X size (arcmin)",[60],more=float,namespace=self),
              TDLOption("_skyvis_ygrid","Y size (arcmin)",[60],more=float,namespace=self)
            ] + self._runtime_vis_options;
        self._runtime_vis_options = [
            TDLOption("_skyvis_npix","Number of x/y pixels",[100],more=int,namespace=self),
            TDLOption("_skyvis_timespan","Time span (hours)",[12],more=float,namespace=self),
            TDLOption("_skyvis_ntime","Number of time planes",[13],more=int,namespace=self),
            TDLOption("_skyvis_freq0","Lowest frequency (MHz)",[20],more=float,namespace=self),
            TDLOption("_skyvis_freq1","Highest frequency (MHz)",[2000],more=float,namespace=self),
            TDLOption("_skyvis_nfreq","Number of freq planes",[10],more=int,namespace=self),
          ] + self._runtime_vis_options;
        self._runtime_options.append(TDLMenu("Visualize sky-Jones terms",
          *self._runtime_vis_options));
    return self._runtime_options;

  def get_inspectors (self):
    """Returns list of inspector nodes created by this MeqMaker""";
    return self._inspectors or [];

  def _make_attr (*comps):
    """Forms up the name of an 'enabled' attribute for a given module group (e.g. a Jones term)"""
    return "_".join(comps);
  _make_attr = staticmethod(_make_attr);

  def _module_togglename (self,label,mod):
    return self._make_attr("enable",label,_modname(mod).replace('.',"_"));

  def _group_togglename (self,label):
    return self._make_attr("enable",label);

  def is_group_enabled (self,label):
    """Returns true if the given group is enabled"""
    return getattr(self,self._group_togglename(label),False);

  def is_module_enabled (self,label,mod):
    """Returns true if the given module is enabled for the given group label"""
    return self.is_group_enabled(label) and getattr(self,self._module_togglename(label,mod),False);

  def _module_selector (self,menutext,label,modules,extra_opts=[],use_toggle=True,exclusive=True,**kw):
    # Forms up a "module selector" submenu for the given module set.
    # for each module, we either take the options returned by module.compile_options(),
    # or pass in the module itself to let the menu "suck in" its options

    # Overall toggle attribute passed to outer option. If we don't use a toggle, set it to True always
    toggle = self._group_togglename(label);
    if not use_toggle:
      setattr(self,toggle,True);
      toggle = None;
    # exclusive option. If True, modules are mutually exclusive
    if exclusive:
      exclusive = self._make_attr("use",label,"module");
    else:
      exclusive = None;
    if len(modules) == 1:
      doc = getattr(modules[0],'__doc__',None);
      modname = _modname(modules[0]);
      mainmenu = TDLMenu(menutext,toggle=toggle,namespace=self,name=modname,doc=doc,
                         *(_modopts(modules[0],'compile')+list(extra_opts)),**kw);
      # set the module's toggle to always on (will be controlled by outer toggle attribute)
      setattr(self,self._module_togglename(label,modules[0]),True);
    else:
      submenus = [ TDLMenu("Use '%s' module"%_modname(mod),name=_modname(mod),
                            toggle=self._module_togglename(label,mod),
                            namespace=self,
                            doc=getattr(mod,'__doc__',None),
                            *_modopts(mod,'compile'))
                    for mod in modules ];
      mainmenu = TDLMenu(menutext,toggle=toggle,exclusive=exclusive,namespace=self,
                          *(submenus+list(extra_opts)),**kw);
    return mainmenu;

  def _get_selected_module (self,label,modules):
    """Get selected module from list.
    Checks the self.<label>_module attribute against module names. This attribute
    is set by the TDLMenu() invocation in _module_selector above, when called
    in exclusive mode.
    """;
    # check global toggle for this module group
    for mod in modules:
      if self.is_module_enabled(label,mod):
        return mod;
    return None;

  def _get_selected_modules (self,label,modules):
    """Get subset of selected modules from list.
    Checks the self.<label>_enable_<module> attribute against module names. This attribute
    is set by the TDLMenu() invocations in _module_selector above, when called
    in non-exclusive mode.
    """;
    # get list of check global toggle for this module group
    if self.is_group_enabled(label):
      return [ mod for mod in modules if self.is_module_enabled(label,mod) ];
    else:
      return [];

  def estimate_image_size (self):
    max_estimate = None;
    for module in self._get_selected_modules('sky',self._sky_models):
      estimate = getattr(module,'estimate_image_size',None);
      if callable(estimate):
        max_estimate = max(max_estimate,estimate());
    return max_estimate;

  def get_source_list (self,ns):
    if self._source_list is None:
      modules = self._get_selected_modules('sky',self._sky_models);
      if not modules:
        return [];
      self._source_list = [];
      for mod in modules:
        self._source_list += mod.source_list(ns);
      # add indices to source attrubutes (used by the %# directive when exporting  annotations)
      for isrc,src in enumerate(self._source_list):
        src.set_attr("#",isrc);
    return self._source_list;

  def add_inspector_node (self,inspector_node,name=None):
    """adds an inspector node to internal list, and creates a bookmark page for it.""";
    self._inspectors.append(inspector_node);
    name = name or inspector_node.initrec().get('title',None) or inspector_node.name.replace('_',' ');
    Meow.Bookmarks.Page(name).add(inspector_node,viewer="Collections Plotter");

  def make_bookmark_set (self,nodes,quals,title_inspector,title_folder=None,inspector_node=None,labels=None,freqmean=True):
    """Makes a standard set of bookmarks for a series of nodes: i.e. an inspector, and a folder of
    individual bookmarks""";
    labels = labels or [ ':'.join([getattr(q,'name',None) or str(q) for q in qq]) for qq in quals ];
    if freqmean:
      nodelist = [ nodes('freqmean',*qq) << Meq.Mean(nodes(*qq),reduction_axes="freq") for qq in quals ];
    else:
      nodelist = [ nodes(*qq) for qq in quals ];
    # make inspector + bookmark
    insp = (inspector_node or nodes('inspector')) << Meq.Composer(dims=[0],plot_label=labels,mt_polling=True,*nodelist);
    self.add_inspector_node(insp,title_inspector);
    # make folder of per-baseline plots
    if title_folder is not None:
      Meow.Bookmarks.make_node_folder(title_folder,[nodes(*qq) for qq in quals],sorted=True,ncol=2,nrow=2);

  def make_per_ifr_bookmarks (self,nodes,title,ifrs=None,freqmean=True):
    """Makes a standard set of bookmarks for a per-ifr node: i.e. an inspector, and a folder of
    per-baseline bookmarks""";
    ifrs = ifrs or Meow.Context.array.ifrs();
    return self.make_bookmark_set(nodes,ifrs,
        "%s: inspector plot"%title,"%s: by interferometer"%title,
        labels=[ "%s-%s"%(p,q) for p,q in ifrs ],freqmean=freqmean);

  def make_per_station_bookmarks (self,nodes,title,stations=None,freqmean=True):
    """Makes a standard set of bookmarks for a station node: i.e. an inspector, and a folder of
    per-station bookmarks""";
    stations = stations or Meow.Context.array.stations();
    return self.make_bookmark_set(nodes,[ (p,) for p in stations],
        "%s: inspector plot"%title,"%s: by station"%title,
        freqmean=freqmean);

  def make_per_source_per_station_bookmarks (self,nodes,title,sources,stations=None,freqmean=True):
    """Makes a standard set of bookmarks for a station node: i.e. an inspector, and a folder of
    per-station bookmarks""";
    stations = stations or Meow.Context.array.stations();
    return self.make_bookmark_set(nodes,[(src,p) for src in sources for p in stations],
        "%s: inspector plot"%title,"%s: by station"%title,
        freqmean=freqmean);

  def _make_skyjones_visualizer_source (self,ns):
    # create visualization nodes
    if not self._skyjones_visualizer_source:
      if self._skyvis_frame == SKYJONES_LM:
        l = ns.lmgrid_l << Meq.Grid(axis="l");
        m = ns.lmgrid_m << Meq.Grid(axis="m");
        visdir = Meow.LMDirection(ns,"visgrid",l,m);
      elif self._skyvis_frame == SKYJONES_RADEC:
        ra = ns.lmgrid_ra << Meq.Grid(axis="m");
        dec = ns.lmgrid_dec << Meq.Grid(axis="l");
        visdir = Meow.Direction(ns,"visgrid",ra,dec);
      elif self._skyvis_frame == SKYJONES_AZEL or self._skyvis_frame == SKYVIS_AZEL_FULL:
        az = ns.lmgrid_az << Meq.Grid(axis="m");
        el = ns.lmgrid_el << Meq.Grid(axis="l");
        visdir = Meow.AzElDirection(ns,"visgrid",az,el);
      # make a source with this direction object
      self._skyjones_visualizer_source = Meow.PointSource(ns,"visgrid",visdir);
    return self._skyjones_visualizer_source;

  def _add_sky_visualizer (self,visnode,visnode_sq,label,name):
    Meow.Bookmarks.Page("%s (%s) visualizer"%(label,name)).add(visnode,viewer="Result Plotter");
    Meow.Bookmarks.Page("|%s^2| (%s) visualizer"%(label,name)).add(visnode_sq,viewer="Result Plotter");
    def visualize_node (mqs,parent,**kw):
      if self._skyvis_ntime > 1:
        dt = self._skyvis_timespan/((self._skyvis_ntime-1)*2.);
      else:
        dt = self._skyvis_ntime = 1;
      if self._skyvis_nfreq > 1:
        df = (self._skyvis_freq1-self._skyvis_freq0)/((self._skyvis_nfreq-1)*2.);
      else:
        df = self._skyvis_nfreq = 1;
      time = [-dt*3600,(self._skyvis_timespan+dt)*3600];
      freq = [(self._skyvis_freq0-df)*1e+6,(self._skyvis_freq1+df)*1e+6];
      if self._skyvis_frame == SKYJONES_LM:
        ARCMIN = math.pi/(180*60);
        dx = self._skyvis_xgrid*ARCMIN/2;
        dy = self._skyvis_ygrid*ARCMIN/2;
        domain = meq.gen_domain(time=time,freq=freq,l=[-dx,dx],m=[-dy,dy]);
      elif self._skyvis_frame in [SKYJONES_RADEC,SKYJONES_AZEL_FULL]:
        domain = meq.gen_domain(time=time,freq=freq,l=[-math.pi/2,math.pi/2],m=[0,2*math.pi]);
      elif self._skyvis_frame == SKYJONES_AZEL:
        domain = meq.gen_domain(time=time,freq=freq,l=[0,math.pi/2],m=[0,2*math.pi]);
      cells = meq.gen_cells(domain,num_time=self._skyvis_ntime,num_freq=self._skyvis_nfreq,
                                   num_l=self._skyvis_npix,num_m=self._skyvis_npix);
      request = meq.request(cells, rqtype='ev');
      mqs.clearcache(visnode_sq.name,recursive=True);
      mqs.clearcache(visnode.name,recursive=True);
      mqs.execute(visnode.name,request);
      mqs.execute(visnode_sq.name,request);
    self._runtime_vis_options.append(TDLJob(visualize_node,"Visualize %s (%s)"%(label,name)));

  def get_uvjones_nodes (self,label):
    """Returns the Jones nodes associated with the given Jones label.
    Returns unqualified node that should be qualified with a station index.
    If this Jones term is not yet initialized, or is disabled via compile-time
    options, returns None. If Jones term is not found, raises KeyError.""";
    for jt in self._uv_jones_list:
      if jt.label == label:
        return jt.base_node;
    raise KeyError,"uv-Jones term %s not defined"%label;

  def get_skyjones_nodes (self,label):
    """Returns the Jones nodes associated with the given Jones label.
    Returns unqualified node that should be qualified with a source and station index.
    If this Jones term is not yet initialized, or is disabled via compile-time
    options, returns None. If Jones term is not found, raises KeyError.""";
    for jt in self._sky_jones_list:
      if jt.label == label:
        return jt.base_node;
    raise KeyError,"sky-Jones term %s not defined"%label;

  def _get_jones_nodes (self,ns,jt,stations,sources=None,solvable_sources=set()):
    """Returns the Jones nodes associated with the given JonesTerm ('jt'). If
    the term has been disabled (through compile-time options), returns None.
    'stations' is a list of stations.
    'sources' is a list of sources (for sky-Jones only).
    Return value is tuple of (basenode,solvable). Basenode should be qualified
    with source (if sky-Jones) and station to obtain the Jones term. Solvable is True
    if the Jones set is to be treated as potentially solvable.
    If basenode==None, the given JonesTerm is not implemented and may be omitted from the ME.
    """;
    if jt.base_node is None:
      if isinstance(jt,self.SkyJonesTerm):
        sources = jt.subset_selector.filter(sources);
        if not sources:
          return None,False;
      all_sources,all_stations = sources,(stations or Context.array.stations());
      # for sky-jones terms, apply subset selector to finalize list of sources
      # are we doing a per-source and per-station Jones term?
      same_all_sources = sources and getattr(self,self._make_attr('all_sources',jt.label),False);
      same_all_stations = getattr(self,self._make_attr('all_stations',jt.label),False);
      if same_all_sources:
        sources = [ all_sources[0] ];
      if same_all_stations:
        stations = [ all_stations[0] ];
      # now get the module
      module = self._get_selected_module(jt.label,jt.modules);
      if not module:
        return None,False;
      prev_num_solvejobs = ParmGroup.num_solvejobs();  # to keep track of whether module creates its own
      jones_inspectors = [];
      jt.solvable = False;
      sources0 = sources;
      skyvis = None;
      # For sky-Jones terms, see if this module has pointing offsets enabled
      if isinstance(jt,self.SkyJonesTerm):
        # is visualization enabled? insert extra source then for the visualization
        if self.use_skyjones_visualizers:
          skyvis = self._make_skyjones_visualizer_source(ns);
          sources = [ skyvis ] + list(sources);
        dlm = None;
        if jt.pointing_modules:
          pointing_module = self._get_selected_module(jt.label+"pe",jt.pointing_modules);
          # get pointing offsets
          if pointing_module:
            dlm = ns[jt.label]('dlm');
            inspectors = [];
            dlm = pointing_module.compute_pointings(dlm,stations=stations,tags=jt.label,
                                                    label=jt.label+'pe',
                                                    inspectors=inspectors,meqmaker=self);
            pointing_solvable = getattr(pointing_module,'solvable',True);
            if inspectors:
              for node in inspectors:
                self.add_inspector_node(node);
            elif dlm and self.use_jones_inspectors:
              self.make_per_station_bookmarks(dlm,"%s pointing errors"%jt.label,stations,freqmean=False);
      else:
        dlm = None;
      # now make the appropriate matrices
      # note that Jones terms are computed using the original source list.
      # this is to keep the extra corruption qualifiers from creeping into their names
      Jj = ns[jt.label];    # create name hint for nodes
      inspectors = [];
      jt.base_node = Jj = module.compute_jones(Jj,sources=sources,stations=stations,
                                          pointing_offsets=dlm,solvable_sources=solvable_sources,
                                          tags=jt.label,label=jt.label,
                                          meqmaker=self,
                                          inspectors=inspectors);
      # ignoring the pointing-solvable attribute for now
      jt.solvable = getattr(module,'solvable',True);
      # Jj will be None if module is not active for some reason
      if Jj:
        # add sky-Jones visualizer if asked to
        if skyvis:
          jones = Jj(skyvis);
          if jones(stations[0]).initialized():
            # since the jones term may be the same for all stations, check the first two for identity.
            # If they are identical, use only the first station instead of a full station list.
            if len(stations)<2 or jones(stations[0]) is jones(stations[1]):
              vis_stations = stations[0:1];
            else:
              vis_stations = stations;
            vis = ns.visualizer(jt.label) << Meq.Composer(dims=[0],*[jones(p) for p in vis_stations]);
            for p in vis_stations:
              jp = jones(p);
              jpsq = jp('sq') << Meq.MatrixMultiply(jp,jp('t') << Meq.ConjTranspose(jones(p)));
              jj = jp('2tr') << (jpsq(11) << Meq.Selector(jpsq,index=0)) + (jpsq(22) << Meq.Selector(jpsq,index=3));
              jp('a2tr') << Meq.Abs(jj);
            vissq = ns.visualizer_sq(jt.label) << Meq.Composer(dims=[0],*[jones(p)('a2tr') for p in vis_stations]);
            self._add_sky_visualizer(vis,vissq,jt.label,jt.name);
          # remove skyvis from list of sources
          sources = sources[1:];
        if inspectors:
          for node in inspectors:
            self.add_inspector_node(node);
        # if module does not make its own inspectors, add automatic ones
        elif self.use_jones_inspectors:
          if sources:
            quals = [ (src,p) for src in sources0 for p in stations if Jj(src,p).initialized() ];
            self.make_bookmark_set(Jj,quals,"%s-Jones: inspector"%jt.label);
          else:
            self.make_per_station_bookmarks(Jj,"%s-Jones"%jt.label,stations);
        # expand into per-station, per-source terms as needed
        if same_all_sources and same_all_stations:
          jj0 = Jj(all_sources[0],all_stations[0]);
          for isrc,src in enumerate(all_sources):
            for ip,p in enumerate(all_stations):
              if isrc or ip:
                Jj(src,p) << Meq.Identity(jj0);
        elif same_all_sources:
          for p in stations:
            jj0 = Jj(all_sources[0],p);
            for src in all_sources[1:]:
              Jj(src,p) << Meq.Identity(jj0);
        elif same_all_stations:
          if sources:
            for src in sources:
              jj0 = Jj(src,all_stations[0]);
              for p in all_stations[1:]:
                Jj(src,p) << Meq.Identity(jj0);
          else:
            jj0 = Jj(all_stations[0]);
            for p in all_stations[1:]:
              Jj(p) << Meq.Identity(jj0);
    return jt.base_node,jt.solvable;

    ## create TDL job
    #def tdljob_compute_visualization (parent,mqs,**kw):
      #cells = meq.gen_cells(domain,num_time=10,num_l=100,num_m=100);
      #request = meq.request(cells,rqtype='ev')
      #mqs.execute('',request);

  def make_predict_tree (self,ns,sources=None,uvdata=None,ifrs=None):
    """makes predict trees using the sky model and ME.
    'ns' is a node scope
    'sources' is a list of sources; the current sky model is used if None.
    'uvdata' is a basenode for additional model uv-data (which should be qualified with a station pair). If supplied, it will be
      added to the sky model, before applying any uv terms.
    Returns a base node which should be qualified with a station pair.
    """;
    stations = Meow.Context.array.stations();
    ifrs = ifrs or Meow.Context.array.ifrs();
    # use sky model if no source list is supplied
    sources = sources or self.get_source_list(ns);
    if self.use_smearing:
      count = self.smearing_count;
      if count == "all":
        count = len(sources);
      for src in sources[0:count]:
        src.enable_smearing();
    # are we using decomposition? Then form up an alternate tree for decomposable sources.
    dec_sky = None;
    if self.use_decomposition:
      # first, split the list into decomposable and non-decomposable sources
      dec_sources = [ src for src in sources if src.is_station_decomposable() ];
      sources = [ src for src in sources if not src.is_station_decomposable() ];
      if dec_sources:
        # for every decomposable source, build a jones chain, and multiply
        # it by that source's sqrt-visibility
        skychain = {};
        for jt in self._sky_jones_list:
          Jj,solvable = self._get_jones_nodes(ns,jt,stations,sources=dec_sources);
          # if this Jones is enabled (Jj not None), add it to chain of each source
          if Jj:
            for src in dec_sources:
              jones = Jj(src.name);
              # only add to chain if this Jones term is initialized
              if jones(stations[0]).initialized():
                chain = skychain.setdefault(src.name,[]);
                chain.insert(0,Jj(src.name));
        # now, form up a chain of uv-Jones terms
        uvchain = [];
        for jt in self._uv_jones_list:
          Jj,solvable = self._get_jones_nodes(ns,jt,stations);
          if Jj:
            uvchain.insert(0,Jj);
        if uvchain:
          # if only one uv-Jones, will use it directly
          if len(uvchain) > 1:
            for p in stations:
              ns.uvjones(p) << Meq.MatrixMultiply(*[j(p) for j in uvchain]);
            uvchain = [ ns.uvjones ];
        # now, form up the corrupt sqrt-visibility (and its conjugate) of each source,
        # then the final visibility
        corrvis = ns.corrupt_vis;
        sqrtcorrvis = ns.sqrt_corrupt_vis;
        sqrtcorrvis_conj = sqrtcorrvis('conj');
        for src in dec_sources:
          # terms is a list of matrices to be multiplied
          sqrtvis = src.sqrt_visibilities();
          C = sqrtcorrvis(src);
          Ct = sqrtcorrvis_conj(src);
          # make a skyjones(src,p) node containing a product of all the sky-Jones
          skchain = skychain.get(src.name,[]);
          if len(skchain) > 1:
            for p in stations:
              ns.skyjones(src,p) << Meq.MatrixMultiply(*[j(p) for j in skchain]);
            skchain = [ ns.skyjones(src) ];
          jones_chain = uvchain + skchain;
          # if there's a real jones chain, multiply all the matrices
          if jones_chain:
            for p in stations:
              C(p) << Meq.MatrixMultiply(*([j(p) for j in jones_chain]+[sqrtvis(p)]));
              if p is not stations[0]:
                Ct(p) << Meq.ConjTranspose(C(p));
          # else use an identity relation
          else:
            for p in stations:
              C(p) << Meq.Identity(sqrtvis(p));
              if p is not stations[0]:
                Ct(p) << Meq.ConjTranspose(C(p));
          # ok, now get the visiblity of each source by multiplying its two per-station contributions
          for p,q in ifrs:
            ns.corrupt_vis(src,p,q) << Meq.MatrixMultiply(C(p),Ct(q));
        # finally, sum up all the source contributions
        # if no non-decomposable sources, then call the output 'visibility:sky', since
        # it already contains everything
        if sources:
          dec_sky = ns.visibility('sky1');
        else:
          dec_sky = ns.visibility('sky');
        Parallelization.add_visibilities(dec_sky,[ns.corrupt_vis(src) for src in dec_sources],ifrs);
        if not sources and not uvdata:
          return self._apply_vpm_list(ns,dec_sky);

    # Now, proceed to build normal trees for non-decomposable sources.
    # An important optimization when solving for sky terms is to put all the sources containing
    # solvable corruptions in one patch, and the rest in another
    # (this makes optimal use of cache when solving for the corruption).
    # This list will contain a list of (name,src) tuples, and will be updated as each sky-Jones term is applied
    sourcelist = [ (src.name,src) for src in sources ];
    # This set will contain names of sources to which solvable sky-Jones were been applied
    solvable_skyjones = set();

    # apply all sky Jones terms
    if sourcelist:
      for jt in self._sky_jones_list:
        Jj,solvable = self._get_jones_nodes(ns,jt,stations,sources=sources,solvable_sources=solvable_skyjones);
        # if this Jones is enabled (Jj not None), corrupt each source
        if Jj:
          newlist = [];
          # print jt.label,solvable,len(solvable_skyjones);
          for name,src in sourcelist:
            # if Jones term is initialized, corrupt and append to corr_sources
            # if not initialized, then append to uncorr_sources
            jones = Jj(name);
            if jones(stations[0]).initialized():
              src = src.corrupt(jones);
              if solvable:
                solvable_skyjones.add(name);
            newlist.append((name,src));
          # update list
          sourcelist = newlist;
    # now make two separate lists of sources with solvable corruptions, and sources without
    corrupted_sources = [ src for name,src in sourcelist if name in solvable_skyjones ];
    uncorrupted_sources = [ src for name,src in sourcelist if name not in solvable_skyjones ];
    # if we're solvable, and both lists are populated, make two patches
    if self._solvable and corrupted_sources and uncorrupted_sources:
      sky_sources = [
        Meow.Patch(ns,'sky-c',Meow.Context.observation.phase_centre,components=corrupted_sources),
        Meow.Patch(ns,'sky-nc',Meow.Context.observation.phase_centre,components=uncorrupted_sources)
      ];
    else:
      sky_sources = corrupted_sources + uncorrupted_sources;

    if uvdata:
      sky_sources.append(Meow.KnownVisComponent(ns,'uvdata',uvdata));

    # If >1 source, form up patch. Call it "sky1" if this is not the final output
    if len(sky_sources) == 1:
      allsky = sky_sources[0];
    else:
      allsky = Meow.Patch(ns,'sky1' if dec_sky else 'sky',Meow.Context.observation.phase_centre,components=sky_sources);

    # add uv-plane effects
    for jt in self._uv_jones_list:
      Jj,solvable = self._get_jones_nodes(ns,jt,stations);
      if Jj:
        allsky = allsky.corrupt(Jj);

    # form up list of visibility contributions
    terms = [ allsky.visibilities() ];
    if dec_sky:
      terms.append(dec_sky);

    # add a Meq.Add node on top of that, if needed
    if len(terms) > 1:
      vis = ns.visibility('sky');
      for p,q in ifrs:
        vis(p,q) << Meq.Add(*[term(p,q) for term in terms]);
    else:
      vis = terms[0];

    # now chain up any visibility processors
    return self._apply_vpm_list(ns,vis,ifrs=ifrs);

  def _apply_vpm_list (self,ns,vis,ifrs=None):
    # chains up any visibility processors, and applies them to the visibilities.
    # Returns new visibilities.
    for vpm in self._uv_vpm_list:
      module = self._get_selected_module(vpm.label,vpm.modules);
      if module:
        inspectors = [];
        nodes = vis(vpm.label);
        if module.process_visibilities(nodes,vis,ns=getattr(ns,vpm.label).Subscope(),
              ifrs=ifrs,meqmaker=self,
              tags=vpm.label,label=vpm.label,inspectors=inspectors) is not None:
          # add inspectors to internal list
          for insp in inspectors:
            self.add_inspector_node(insp);
          vis = nodes;
    return vis;

  def make_tree (self,*args,**kw):
    print "Your script uses the deprecated MeqMaker.make_tree() method. Please change it to use make_predict_tree()";
    return self.make_predict_tree(*args,**kw); # alias for compatibility with older code

  def corrupt_uv_data (self,ns,uvdata,ifrs=None,label="uvdata"):
    """Corrupts the visibilities given by the 'uvdata' nodes by the current uv-Jones chain.
    'label' is a label that will be applied to any intermediate visibilty nodes, default is
    'uvdata'.
    Returns an unqualified node that must be qualified with a station pair to get visibilities.
    """;
    stations = Meow.Context.array.stations();
    ifrs = ifrs or Meow.Context.array.ifrs();
    src = Meow.KnownVisComponent(ns,label,uvdata);
    # apply all uv-Jones
    for jt in self._uv_jones_list:
      Jj,solvable = self._get_jones_nodes(ns,jt,stations);
      if Jj:
        src = src.corrupt(Jj);
    # now chain up any visibility processors
    return self._apply_vpm_list(ns,src.visibilities(),ifrs=ifrs);

  def correct_uv_data (self,ns,inputs,outputs=None,sky_correct=None,inspect_ifrs=None,flag_jones_minmax=None):
    """makes subtrees for correcting the uv data given by 'inputs'.
    If 'outputs' is given, then it will be qualified by a jones label and by stations pairs
    to derive the output nodes. If it is None, then ns.correct(jones_label) is used as a base name.
    By default only uv-Jones corrections are applied, but if 'sky_correct' is set to
    a source object (or source name), then sky-Jones corrections for this particular source
    are also put in.
    If 'flag_jones_minmax' is given, it must be a (min,max) tuple. Flags will be assigned to Jones corrections that
    are above or below the indicated values.
      NB: The source/name given by 'sky_correct' here should have been present in the source list
      used to invoke make_predict_tree().
    Returns an unqualified node that must be qualified with a station pair to get visibilities.
    """;
    stations = Meow.Context.array.stations();
    ifrs = Meow.Context.array.ifrs();
    inspect_ifrs = inspect_ifrs or ifrs;

    # apply vpm corrections, if any
    for vpm in self._uv_vpm_list:
      module = self._get_selected_module(vpm.label,vpm.modules);
      if module and hasattr(module,'correct_visibilities'):
        inspectors = [];
        nodes = inputs(vpm.label);
        if module.correct_visibilities(nodes,inputs,ns=getattr(ns,vpm.label).Subscope(),
             tags=vpm.label,label=vpm.label,inspectors=inspectors,meqmaker=self) is not None:
          # add inspectors to internal list
          for insp in inspectors:
            self.add_inspector_nodes(insp);
          inputs = nodes;

    # now build up a correction chain for every station
    correction_chains = dict([(p,[]) for p in stations]);

    # first, collect all sky Jones terms
    if sky_correct is not None:
      for jt in self._sky_jones_list:
        # if using coherency decomposition, we will already have defined a "uvjones" node
        # containing a product of all the sky-Jones terms, so use that
        skyjones = ns.skyjones(sky_correct);
        if skyjones(stations[0]).initialized():
          for p in stations:
            correction_chains[p].insert(0,skyjones(p));
        else:
          Jj,solvable = self._get_jones_nodes(ns,jt,stations,sources=[sky_correct]);
          if Jj:
            Jj = Jj(sky_correct);
            for p in stations:
              correction_chains[p].insert(0,Jj(p));

    # now collect all uv-Jones and add them to the chains
    # if using coherency decomposition, we will already have defined a "uvjones" node
    # containing a product of all the uv-Jones terms, so use that
    if ns.uvjones(stations[0]).initialized():
      for p in stations:
        correction_chains[p].insert(0,ns.uvjones(p));
    else:
      for jt in self._uv_jones_list:
        Jj,solvable = self._get_jones_nodes(ns,jt,stations);
        if Jj:
          for p in stations:
            correction_chains[p].insert(0,Jj(p));

    # get base node for output visibilities. The variable will be replaced by a new name
    if outputs is None:
      outputs = ns.correct;

    # invert and conjugate the products of the correction chains
    # NB: this really needs some more thought, since different combinations of time/freq dependencies
    # will change the optimal order in which matrices should be inverted and multiplied.
    # We'll use outputs:Jinv:p and outputs:Jtinv:p for the (Jn...J1)^{-1} and (Jn...J1)^{-1}^t
    # products
    Jinv = outputs('Jinv');
    Jtinv = outputs('Jtinv');
    if len(correction_chains[stations[0]]) > 1:
      Jprod = outputs('Jprod');
      for p in stations:
        jj = Jprod(p) << Meq.MatrixMultiply(*correction_chains[p]);
        if flag_jones_minmax:
          jj = StdTrees.make_jones_norm_flagger(jj,flag_jones_minmax[0],flag_jones_minmax[1],flagmask=MSUtils.FLAGMASK_OUTPUT);
        Jinv(p) << Meq.MatrixInvert22(jj);
        if p != stations[0]:
          Jtinv(p) << Meq.ConjTranspose(Jinv(p));
    elif correction_chains[stations[0]]:
      for p in stations:
        jj = correction_chains[p][0];
        if flag_jones_minmax:
          jj = StdTrees.make_jones_norm_flagger(jj,flag_jones_minmax[0],flag_jones_minmax[1],flagmask=MSUtils.FLAGMASK_OUTPUT);
        Jinv(p) << Meq.MatrixInvert22(jj);
        if p != stations[0]:
          Jtinv(p) << Meq.ConjTranspose(Jinv(p));
    else:
      for p,q in ifrs:
        outputs(p,q) << Meq.Identity(inputs(p,q));
      # make an inspector for the results
      self.make_per_ifr_bookmarks(outputs,"Output visibilities",ifrs=inspect_ifrs);
      return outputs;
    # now apply the correction matrices
    for p,q in ifrs:
      outputs(p,q) << Meq.MatrixMultiply(Jinv(p),inputs(p,q),Jtinv(q));

    # make an inspector for the results
    return outputs;

  def close (self):
    """cleans up the MeqMaker""";
    # export karma annotations if asked to
    if self.export_kvis and self._source_list and self.export_karma_filename:
      # how many sources get labels
      nlab = self.export_karma_nlabel;
      if nlab == "all":
        nlab = len(self._source_list);
      # incremental mode -- split sources into multiple annotation files
      if self.export_karma_incremental and self.export_karma_incremental != "all":
        filename,ext = os.path.splitext(self.export_karma_filename);
        for i0 in range(0,len(self._source_list),self.export_karma_incremental):
          i1 = min(i0+self.export_karma_incremental,len(self._source_list));
          export_karma_annotations(self._source_list[i0:i1],
            filename="%s-%d-%d%s"%(filename,i0,i1-1,ext),
            label_format=self.export_karma_label,maxlabels=nlab);
          nlab = max(nlab-self.export_karma_incremental,0);
      # normal mode -- all sources go into one file
      else:
        export_karma_annotations(self._source_list,
            self.export_karma_filename,label_format=self.export_karma_label,maxlabels=nlab);

class SourceSubsetSelector (object):
  def __init__ (self,title,tdloption_namespace=None,doc=None,annotate=True):
    if tdloption_namespace:
      self.tdloption_namespace = tdloption_namespace;
    self.subset_enabled = False;
    global _annotation_label_doc;
    opts = [
        TDLOption('source_subset',"Sources",["all"],more=str,namespace=self,
                  doc="""<P>You can enter source names separated by space. "all" selects all sources. "=<i>tagname</i>"
                  selects all sources with a given tag. "-<i>name</i>" deselects the named sources. "-=<i>tagname</i>" deselects sources by tag,</P>""")
    ];
    self.annotate = False;
    if annotate:
      opts.append(TDLMenu("Annotate selected sources",doc="""If you're exporting annotations for your sky
                    model, you may choose to mark the selected sources using different colours and/or labels""",toggle='annotate',namespace=self,
          *(
            TDLOption('annotation_symbol_color',"Colour of source markers",
                      ["default","red","green","blue","yellow","purple","cyan","white"],more=str,namespace=self),
            TDLOption('annotation_label_color',"Colour of labels",
                      ["default","red","green","blue","yellow","purple","cyan","white"],more=str,namespace=self),
            TDLOption('annotation_label',"Label format",
                      ["default",None,"%N","%N I=%(I).3g","%N Ia=%(Iapp).3g"],more=str,namespace=self,
                      doc="Overrides default annotation labels for selected sources.\n"+
                          _annotation_label_doc)
          )
      ));
    subset_opt = TDLMenu(title,toggle='subset_enabled',namespace=self,doc=doc,*opts);
    self._src_set = None;
    self.options = [ subset_opt ];

  def selectedTag (self):
    """If sources were selected by tag, returns the tagname. Else returns None""";
    if self.source_subset[0] == "=":
      return self.source_subset.split()[0][1:];
    return None;

  def filter (self,srclist0):
    if not self.subset_enabled or self.source_subset == "all":
      return srclist0;
    all = set([src.name for src in srclist0]);
    srcs = set();
    for ispec,spec in enumerate(self.source_subset.split()):
      spec = spec.strip();
      if spec:
        # if first spec is a negation, then implictly select all sources first
        if not ispec and spec[0] == "-":
            srcs = all;
        if spec == "all":
          srcs = all;
        elif spec.startswith("="):
          srcs.update([ src for src in srclist if src.get_attr(spec[1:]) ]);
        elif spec.startswith("-="):
          srcs.difference_update([ src for src in srclist if src.get_attr(spec[2:]) ]);
        elif spec.startswith("-"):
          srcs.discard(spec[1:]);
    # make list
    srclist = [ src for src in srclist0 if src.name in srcs ];
    if self.annotate:
      # apply annotations
      if self.annotation_symbol_color != "default":
        [ src.set_attr('ANNOTATION_SYMBOL_COLOR',self.annotation_symbol_color) for src in srclist ];
      if self.annotation_label_color != "default":
        [ src.set_attr('ANNOTATION_LABEL_COLOR',self.annotation_label_color) for src in srclist ];
      if self.annotation_label != "default":
        [ src.set_attr('ANNOTATION_LABEL',self.annotation_label) for src in srclist ];
    return srclist;

def export_karma_annotations (sources,filename,
      label_format="%N",
      sym_color='yellow',
      lbl_color='blue',maxlabels=None):
  """Exports a list of sources as a Karma annotations file.
  label_format is used to generate source labels.\n"""+_annotation_label_doc;
  if maxlabels is None:
    maxlabels = len(sources);
  f = file(filename,'wt');
  f.write('COORD W\nPA STANDARD\nCOLOR GREEN\nFONT hershey12\n');
  # calculate default size for crosses
  xcross=0.01;
  for isrc,src in enumerate(sources):
    # can only plot static-position sources
    radec = src.direction.radec_static();
    if not radec:
      continue;
    # convert position to degrees
    ra_d=radec[0]/math.pi*180.0;
    dec_d=radec[1]/math.pi*180.0;
    # write symbol for source
    color = src.get_attr("ANNOTATION_SYMBOL_COLOR",sym_color);
    f.write('COLOR %s\n'%color);
    # GaussianSources get ellipses
    if isinstance(src,Meow.GaussianSource):
      if src.is_symmetric():
        sx = sy = src.get_value("sigma") or 0;
      else:
        sx = src.get_value("sigma1") or 0;
        sy = src.get_value("sigma2") or 0;
      pa = src.get_value("phi") or 0;
      f.write('ELLIPSE %.12f %.12f %f %f %f\n'%(ra_d,dec_d,sx/math.pi*180.0,sy/math.pi*180.0,-pa/math.pi*180.0));
    # else assume point source and do a cross
    else:
      f.write('CROSS %.12f %.12f %f %f\n'%(ra_d,dec_d,xcross,xcross));
    # now generate label
    label = src.get_attr('ANNOTATION_LABEL',label_format);
    if label and isrc < maxlabels:
      attrs = dict(src.attrs);
      # populate attribute dict
      attrs['%'] = '%';
      attrs['N'] = src.name;
      for st in "IQUV":
        attrs[st] = src.get_value(st,default=0);
      lmn = src.direction.lmn_static();
      if lmn:
        R = math.sqrt(lmn[0]**2+lmn[1]**2);
        attrs['Rr'] = R;
        R = R/math.pi*180;
        attrs['Rd'] = R;
        attrs['Rm'] = R*60;
        attrs['Rs'] = R*3600;
      # sort attributes by decreasing length (to make sure that %xyz is replaced before %xy).
      attrkeys = attrs.keys();
      attrkeys.sort(lambda a,b:cmp(len(b),len(a)));
      # now for each attribute, replace entry in format string
      for key in attrkeys:
        # replace simple %attr directives
        label = label.replace('%'+key,str(attrs[key]));
      # now replace remaining entries
      try:
        label = label%attrs;
      except:
        pass;
      # newlines
      label = label.replace("\n","\\\n");
      if label:
        color = src.get_attr("ANNOTATION_LABEL_COLOR",lbl_color);
        f.write('COLOR %s\n'%color);
        f.write('TEXT %.12f %.12f %s\n'%(ra_d,dec_d,label));
  f.close()
