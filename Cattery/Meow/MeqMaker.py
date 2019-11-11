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

 # standard preamble
from Timba.TDL import *
from Timba.Meq import meq
import math
import inspect
import types
import re
import fnmatch

import Meow
from Meow import StdTrees,ParmGroup,Parallelization,MSUtils

DEG = math.pi/180.;

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

PER_SOURCE = "each source";
PER_ALL_SOURCES = "entire model";

class MeqMaker (object):
  def __init__ (self,namespace='me',solvable=False,
      use_correction=False,
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
    self.use_correction = use_correction;
    self._solvable = solvable;
    self._inspectors = [];
    self._skyjones_visualizer_mux = [];
    self._skyjones_visualizer_source = False;
    self._module_toggles = {};
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
      TDLMenu("Include time & bandwidth smearing",namespace=self,toggle='use_smearing',default=True,
        *self.smearing_options())
    );

    other_opt += Meow.IfrArray.compile_options();
    self._compile_options += [ TDLMenu("Measurement Equation options",*other_opt),TDLOptionSeparator("Image-plane components") ];
    # will be True once a uv-plane module has been added
    self._have_uvplane_modules = False;

  def compile_options (self):
    return self._compile_options;

  def smearing_options (self):
    return [ TDLOption('smearing_count',"Restrict to N brightest sources only",["all",10,100],more=int,namespace=self,
      doc="""<P>Smearing is somewhat expensive to calculate in this implementation, so you may choose to limit it
      to some number of brightest sources (e.g. 50)</P>""") ];

  def add_sky_models (self,modules,export_karma=False):
    if self._sky_models:
      raise RuntimeError("add_sky_models() may only be called once");
    self._compile_options.append(
        self._module_selector("Sky model","sky",modules,use_toggle=False,exclusive=False));
    self._sky_models = modules;
    global _annotation_label_doc;
    if export_karma:
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
    else:
      self.export_kvis = False;

  def add_uv_jones (self,label,name=None,modules=None,pointing=None,flaggable=False):
    return self._add_jones_modules(label,name,modules,is_uvplane=True,flaggable=flaggable,pointing=pointing);

  def add_vis_proc_module (self,label,name=None,modules=None):
    return self._add_vpm_modules(label,name,modules);

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
    FLAG_ABS  = "abs";
    FLAG_ABSDIAG  = "abs-diag";
    FLAG_NORM = "norm";
    def __init__ (self,label,name,modules,flaggable=False):
      self.label      = label;
      self.name       = name;
      self.modules    = modules;
      self.base_node  = None;
      self.solvable   = False;
      self.flaggable  = flaggable;

  class SkyJonesTerm (JonesTerm):
    """SkyJonesTerm represents a sky-Jones term.
    It adds a pointing_modules field, holding a list of poiting error modules.
    """;
    def __init__ (self,label,name,modules,pointing_modules=None,flaggable=False):
      MeqMaker.JonesTerm.__init__(self,label,name,modules,flaggable=flaggable);
      self.pointing_modules = pointing_modules;
      self.subset_selector = None;
      self.base_pe_node = None;
      self.pe_initialized = False;

  def _add_module_to_compile_options (self,option,is_uvplane):
    if is_uvplane:
      if not self._have_uvplane_modules:
        self._have_uvplane_modules = True;
        self._compile_options.append(TDLOptionSeparator("UV-plane components"));
    else:
      if self._have_uvplane_modules:
        raise RuntimeError("cannot add sky-plane terms after uv-plane terms");
    self._compile_options.append(option);

  def _add_vpm_modules (self,label,name,modules):
    if type(label) is types.ModuleType:
      modules = label;
      label = getattr(modules,'__default_label__',None);
      name = getattr(modules,'__default_name__','');
      if label is None:
        raise RuntimeError("""Module '%s' does not provide a __default_label__ attribute,
              so must be added with an explicit label"""%modules.__name__);
    elif not modules:
      raise RuntimeError("No modules specified for %s"%name);
    if not isinstance(modules,(list,tuple)):
      modules = [ modules ];
    # Add to internal list
    term = self.VPMTerm(label,name,modules);
    if not self._uv_vpm_list:
      self._compile_options.append(TDLOptionSeparator("Other visibility processing components"));
    self._uv_vpm_list.append(term);
    # make option menus for selecting a visiblity processor
    mainmenu = self._module_selector("Use %s"%name,label,modules);
    self._add_module_to_compile_options(mainmenu,True);

  def _add_jones_modules (self,label,name,modules,is_uvplane=True,flaggable=False,pointing=None):
    # see if only a module is specified instead of the full label,name,modules
    # syntax
    if type(label) is types.ModuleType:
      modules = label;
      label = getattr(modules,'__default_label__',None);
      name = getattr(modules,'__default_name__','');
      if label is None:
        raise RuntimeError("""Module '%s' does not provide a __default_label__ attribute,
              so must be added with an explicit label"""%modules.__name__);
      elif not modules:
        raise RuntimeError("No modules specified for %s Jones"%label);
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
    advanced_options = [ TDLOption(self._make_attr('all_stations',label),"Use same %s-Jones for all stations"%label,False,
        namespace=self,nonexclusive=True,
        doc="""If checked, then the same %s-Jones term will be used for all stations in the array. This may
        make your trees smaller and/or faster. Whether this is a valid approximation or not depends on your
        physics; typically this is valid if your array is small, and the effect arises far from the receivers."""%label) ];
    if is_uvplane:
      # force-set the advanced options to True, since we omit a control for them (unlike the sky-Jones case)
      setattr(self,self._make_attr('advanced',label),True);
#      setattr(self,self._make_attr('all_sources',label),False);
      jt = self.JonesTerm(label,name,modules,flaggable=flaggable);
      if flaggable:
        advanced_options += [
          TDLMenu("Allow flagging on out-of-bounds %s-Jones"%label,
            TDLOption(self._make_attr("flag_jones_type",label),
              "Flag on what quantities",
              { self.JonesTerm.FLAG_ABS:"all elements |%sij|"%label,
                self.JonesTerm.FLAG_ABSDIAG:"diagonal elements |%sii|"%label,
                self.JonesTerm.FLAG_NORM:"matrix norm ||%s||"%label},
              namespace=self,doc=
                """<P>This determines what quantitites are actually compared to the specified bounds.
                This can be either the absolute value of each matrix element
                <NOBR>(|%Jij|)</NOBR>, or else the diagonal elements only <NOBR>(|%Jii|)</NOBR>, or else the matrix norm
                <NOBR>(||%J||)</NOBR>. The matrix norm is defined as:</P>

                <P align="center"><BIG><I>&#x2225;%J&#x2225; = </I>tr<I>(%J%J<sup>H</sup>)<sup>&frac12;</sup>,</I></BIG></P>

                <P>where <I><BIG>%J<sup>H</sup></BIG></I> is the conjugate transpose, and </I>tr()</I> is the trace
                operator.</P>""".replace("%J",label)),
            TDLOption(self._make_attr("flag_jones_freqmean",label),"Take mean over frequency axis",False,namespace=self,doc=
                """<P>Enable this to take the mean over the frequency axis of the specified quantity before checking whether
                it is within bounds. This implies that entire timeslots will be flagged.
                </P>"""),
            TDLOption(self._make_attr("flag_jones_max",label),"Upper bound",[None,10,100],more=float,namespace=self),
            TDLOption(self._make_attr("flag_jones_min",label),"Lower bound",[None,.1,.01],more=float,namespace=self),
          toggle=self._make_attr("flag_jones",label),namespace=self,doc=
            """<P>Enable this to flag output visibilities when some metric characterizing the %J-Jones term goes out of bounds.
            (Note also that there's probably a top-level option that controls Jones-based flagging as a whole -- this must also
            be enabled.)</P>""".replace("%J",label)
        )];
    else:
      jt = self.SkyJonesTerm(label,name,modules,pointing);
      tdloption_namespace = "%s.%s"%(self.tdloption_namespace,self._make_attr('subset',label));
      jt.subset_selector = SourceSubsetSelector("Apply %s-Jones to subset of sources"%label,tdloption_namespace,annotate=False);
      advanced_options += [
        TDLOption(self._make_attr('per_source',label),"Use a unique %s-Jones term per"%label,[PER_SOURCE,PER_ALL_SOURCES],
                more=str,namespace=self,doc=
                """<P>If you specify '%s', each source will receive an independent %s-Jones term.
                If you specify '%s', the same %s-Jones term will be used for the entire sky (making it
                effectively direction-independent).</P>
                <P>Alternatively, you may specify a tag name to group sources by tag. For example,
                if you enter the name of "foo", and source A does not have the tag "foo",
                it will receive its own %s-Jones term. If source B then has the tag "foo=A",
                then it will use the Jones term associated with source A. This is commonly
                used with the "cluster" tag assigned by Tigger (see the tigger-convert script
                for details) to associate a single Jones term with a whole cluster of closely
                located sources. </P>"""%(PER_SOURCE,label,PER_ALL_SOURCES,label,label)) ];
      if self.use_correction:
        advanced_options += [
          TDLOption(self._make_attr('skip_correct',label),"Do not correct for %s-Jones"%label,False,
                  namespace=self,doc=
                  """<P>With this option you may exclude this particular Jones term from the overall
                  sky-Jones correction (if the latter is enabled in the top-level menu).</P>""") ];
      advanced_options += jt.subset_selector.options;
    # add advanced options to menu
    extra_options.append( TDLMenu("Advanced options",
                          toggle=self._make_attr('advanced',label),nonexclusive=True,
                          doc="""<P>This contains some advanced
                          options relevant to the %s-Jones term. Note that the
                          "Advanced options" checkbox itself must be checked for the options
                          within to have effect.</P>"""%label,
                          namespace=self,
                          *advanced_options) );
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

  def _make_attr (c1,label,*comps):
    """Forms up the name of an 'enabled' attribute for a given module group (e.g. a Jones term)"""
    return "_".join([label,c1]+list(comps));
  _make_attr = staticmethod(_make_attr);

  def _group_togglename (self,label):
    return self._make_attr("enable",label);

  def is_group_enabled (self,label):
    """Returns true if the given group is enabled"""
    return getattr(self,self._group_togglename(label),False);

  def is_module_enabled (self,label,mod):
    """Returns true if the given module is enabled for the given group label"""
    return self.is_group_enabled(label) and self._module_toggles.get(label,{}).get(_modname(mod).replace(".","_"));

  def are_advanced_options_enabled (self,label):
    return getattr(self,self._make_attr("advanced",label),False);

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
      exclusive = self._make_attr("module",label);
    else:
      exclusive = None;
    if len(modules) == 1:
      # make single entry
      doc = getattr(modules[0],'__doc__',None);
      modname = _modname(modules[0]);
      mainmenu = TDLMenu(menutext,toggle=toggle,namespace=self,name=modname,doc=doc,
                         *(_modopts(modules[0],'compile')+list(extra_opts)),**kw);
      # set the module's toggle to always on (will be controlled by outer toggle attribute)
      self._module_toggles[label] = { modname.replace(".","_"):True };
    else:
      # make menu of module options
      self._module_toggles[label] = dict(tdloption_namespace=self.tdloption_namespace+"."+label);
      submenus = [ TDLMenu("Use '%s' module"%_modname(mod),name=_modname(mod),
                            toggle=_modname(mod).replace(".","_"),
                            namespace=self._module_toggles[label],
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
    max_estimate = 0;
    for module in self._get_selected_modules('sky',self._sky_models):
      estimate = getattr(module,'estimate_image_size', 0);
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

  def make_bookmark_set (self,nodes,quals,title_inspector,title_folder=None,inspector_node=None,labels=None,freqmean=True,
   vells_labels=True):
    """Makes a standard set of bookmarks for a series of nodes: i.e. an inspector, and a folder of
    individual bookmarks""";
    labels = labels or [ ':'.join([getattr(q,'name',None) or str(q) for q in qq]) for qq in quals ];
    if freqmean:
      nodelist = [ nodes('freqmean',*qq) << Meq.Mean(nodes(*qq),reduction_axes="freq") for qq in quals ];
    else:
      nodelist = [ nodes(*qq) for qq in quals ];
    if vells_labels is True:
      vells_labels = Meow.Context.correlations or [];
    elif not vells_labels:
      vells_labels = [];
    # make inspector + bookmark
    insp = (inspector_node or nodes('inspector')) << Meq.Composer(dims=[0],plot_label=labels,vells_label=vells_labels,mt_polling=True,*nodelist);
    self.add_inspector_node(insp,title_inspector);
    # make folder of per-baseline plots
    if title_folder is not None:
      Meow.Bookmarks.make_node_folder(title_folder,[nodes(*qq) for qq in quals],sorted=True,ncol=2,nrow=2);

  def make_per_ifr_bookmarks (self,nodes,title,ifrs=None,freqmean=True,vells_labels=True):
    """Makes a standard set of bookmarks for a per-ifr node: i.e. an inspector, and a folder of
    per-baseline bookmarks""";
    ifrs = ifrs or Meow.Context.array.ifrs();
    return self.make_bookmark_set(nodes,ifrs,
        "%s: inspector plot"%title,"%s: by interferometer"%title,vells_labels=vells_labels,
        labels=[ "%s-%s"%(p,q) for p,q in ifrs ],freqmean=freqmean);

  def make_per_station_bookmarks (self,nodes,title,
        stations=None,vells_labels=True,freqmean=True):
    """Makes a standard set of bookmarks for a station node: i.e. an inspector, and a folder of
    per-station bookmarks""";
    stations = stations or Meow.Context.array.stations();
    return self.make_bookmark_set(nodes,[ (p,) for p in stations],
        "%s: inspector plot"%title,"%s: by station"%title,
        vells_labels=vells_labels,freqmean=freqmean);

  def make_per_source_per_station_bookmarks (self,nodes,title,sources,
        stations=None,vells_labels=True,freqmean=True):
    """Makes a standard set of bookmarks for a station node: i.e. an inspector, and a folder of
    per-station bookmarks""";
    stations = stations or Meow.Context.array.stations();
    return self.make_bookmark_set(nodes,[(src,p) for src in sources for p in stations],
        "%s: inspector plot"%title,"%s: by station"%title,
        vells_labels=vells_labels,freqmean=freqmean);

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
    raise KeyError("uv-Jones term %s not defined"%label);

  def get_skyjones_nodes (self,label):
    """Returns the Jones nodes associated with the given Jones label.
    Returns unqualified node that should be qualified with a source and station index.
    If this Jones term is not yet initialized, or is disabled via compile-time
    options, returns None. If Jones term is not found, raises KeyError.""";
    for jt in self._sky_jones_list:
      if jt.label == label:
        return jt.base_node;
    raise KeyError("sky-Jones term %s not defined"%label);

  def _get_pointing_error_nodes (self,ns,jt,stations):
    if not jt.pointing_modules:
      return None;
    if not jt.pe_initialized:
      jt.pe_initialized = True;
      pointing_module = self._get_selected_module(jt.label+"pe",jt.pointing_modules);
      if pointing_module:
        dlm = ns[jt.label]('dlm');
        inspectors = [];
        jt.base_pe_node = dlm = pointing_module.compute_pointings(dlm,stations=stations,tags=jt.label,
                                                label=jt.label+'pe',
                                                inspectors=inspectors,meqmaker=self);
        # add inspectors, if any were returned
        if inspectors:
          for node in inspectors:
            self.add_inspector_node(node);
        elif dlm and self.use_jones_inspectors:
          self.make_per_station_bookmarks(dlm,"%s pointing errors"%jt.label,stations,
                                           vells_labels=("dl","dm"),freqmean=False);
    return jt.base_pe_node;

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
      advanced_opt = self.are_advanced_options_enabled(jt.label);
      if isinstance(jt,self.SkyJonesTerm):
        # for sky-jones terms, apply subset selector to finalize list of sources
        if advanced_opt:
          sources = jt.subset_selector.filter(sources);
        if not sources:
          return None,False;
        # Depending on the per_sources option, make mapping between source names and the names of the Jones nodes
        # actually used.
        # sources_with_jones_nodes will be a list of sources for which unique Jones terms will actually be made.
        # jones_mapping[name] will give the name of a source from the above list, which is associated with this name
        per_source = advanced_opt and getattr(self,self._make_attr('per_source',jt.label),None);
        if per_source is PER_SOURCE or not per_source:
          sources_with_jones_node = sources;
          jones_mapping = dict();
        elif per_source is PER_ALL_SOURCES:
          sources_with_jones_node = [ sources[0] ];
          jones_mapping = dict([(src.name,sources[0].name) for src in sources[1:]]);
        else:
          uniqset = set([ (src.get_attr(per_source) or src.name) for src in sources ]);
          sources_with_jones_node = [ src for src in sources if src.name in uniqset ];
          srcdict = dict([(src.name,src) for src in sources_with_jones_node]);
          jones_mapping = dict();
          for src in sources:
            if src.name not in srcdict:
              tagval = src.get_attr(per_source);  # non-zero since otherwise it would be in srcdict already
              if not tagval in srcdict:
                sources_with_jones_node.append(src);
                srcdict[tagval] = src;
              else:
                jones_mapping[src.name] = srcdict[tagval];
      else:
        sources_with_jones_node = None;
      # are we doing a per-station Jones term?
      all_stations = stations or Context.array.stations();
      same_all_stations = advanced_opt and getattr(self,self._make_attr('all_stations',jt.label),False);
      if same_all_stations:
        stations = [ all_stations[0] ];
      # now get the module
      module = self._get_selected_module(jt.label,jt.modules);
      if not module:
        return None,False;
      jones_inspectors = [];
      jt.solvable = False;
      skyvis = None;
      # For sky-Jones terms, see if this module has pointing offsets enabled
      if isinstance(jt,self.SkyJonesTerm):
        sources0 = sources_with_jones_node;
        # is visualization enabled? insert extra source then for the visualization
        if self.use_skyjones_visualizers:
          skyvis = self._make_skyjones_visualizer_source(ns);
          sources_with_jones_node.insert(0,skyvis);
        dlm = self._get_pointing_error_nodes(ns,jt,stations);
      else:
        dlm = None;
      # now make the appropriate matrices
      # note that Jones terms are computed using the original source list.
      # this is to keep the extra corruption qualifiers from creeping into their names
      Jj = ns[jt.label];    # create name hint for nodes
      inspectors = [];
      solvable_sources_from_module = set();
      jt.base_node = Jj = module.compute_jones(Jj,sources=sources_with_jones_node,stations=stations,
                                          pointing_offsets=dlm,solvable_sources=solvable_sources_from_module,
                                          tags=jt.label,label=jt.label,
                                          meqmaker=self,
                                          inspectors=inspectors);
      # ignoring the pointing-solvable attribute for now
      jt.solvable = getattr(module,'solvable',True);
      if isinstance(jt,self.SkyJonesTerm) and Jj and jt.solvable:
        if solvable_sources_from_module:
          solvable_sources.update(solvable_sources_from_module);
        else:
          solvable_sources.update(sources_with_jones_node);
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
          sources_with_jones_node = sources0;
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
        # if same for all stations, loop over other stations and create identity nodes
        if same_all_stations:
          if sources_with_jones_node:
            for src in sources_with_jones_node:
              for p in all_stations[1:]:
                Jj(src,p) << Meq.Identity(Jj(src,all_stations[0]));
          else:
            for p in all_stations[1:]:
              Jj(p) << Meq.Identity(Jj(all_stations[0]));
        # if some sources need to map to these Jones nodes, do it here
        if sources_with_jones_node:
          for src,src0 in jones_mapping.items():
            for p in all_stations:
              Jj(src,p) << Meq.Identity(Jj(src0,p));
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
    Meow.Context.array.enable_uvw_derivatives(self.use_smearing);
    stations = Meow.Context.array.stations();
    ifrs = ifrs or Meow.Context.array.ifrs();
    # use sky model if no source list is supplied
    sources = sources if sources is not None else self.get_source_list(ns);
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
      # make list of (J,solvable) pairs per each enabled Jones term
      joneslist = [];
      for jt in self._sky_jones_list:
        Jj,solvable = self._get_jones_nodes(ns,jt,stations,sources=sources,solvable_sources=solvable_skyjones);
        # if this Jones is enabled (Jj not None), corrupt each source
        if Jj:
          joneslist.append((Jj,solvable));
      # now make multiplication node per each source
      for name,src in sourcelist:
        # get KJones and smear factor for source
        if not src.direction.is_phase_centre():
          Kj = src.direction.KJones();
          smear = src.is_smeared() and src.smear_factor();
        else:
          Kj = smear = None;
        # add solvables
        if src.get_solvables():
          solvable_skyjones.add(name);
        coh = src.coherency();
        # loop over baselines
        for p,q in ifrs:
          mulops = [ coh(p,q) ];
          # add KJones and smear factor to list of multiplication operands
          if Kj:
            mulops.insert(0,Kj(p));
            mulops.append(Kj(q)('conj') ** Meq.ConjTranspose(Kj(q)));
            if smear:
              mulops.insert(0,smear(p,q));
          # add other Jones terms to list
          for Jones,solvable in joneslist:
            J = Jones(name);
            if J(p).initialized():
              mulops.insert(0,J(p));
              mulops.append(J(q,'conj') ** Meq.ConjTranspose(J(q)));
              solvable and solvable_skyjones.add(name);
          # now make multiply node to compute sky visibility
          if len(mulops) > 1:
            ns.sky(name,p,q) << Meq.MatrixMultiply(*mulops);
          else:
            ns.sky(name,p,q) << Meq.Identity(*mulops);
    # now make two separate lists of sources with solvable corruptions, and sources without
    corrupted_sources   =  [ ns.sky(name) for name,src in sourcelist if name in solvable_skyjones ];
    uncorrupted_sources =  [ ns.sky(name) for name,src in sourcelist if name not in solvable_skyjones ];
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
    if len(sky_sources) > 1 or len(sky_sources) == 0:
      sky = Meow.Patch(ns,'sky1' if dec_sky else 'sky',Meow.Context.observation.phase_centre,components=sky_sources);
      skyvis = sky.visibilities();
    else:
      skyvis = sky_sources[0].visibilities() if isinstance(sky_sources[0],Meow.SkyComponent) else sky_sources[0];

    # add uv-plane effects
    # first make list of Jones basenodes
    joneslist = [];
    for jt in self._uv_jones_list:
      Jj,solvable = self._get_jones_nodes(ns,jt,stations);
      if Jj:
        joneslist.append(Jj);
    # now make matrix multiplication nodes
    if joneslist:
      for p,q in ifrs:
        mulops = [ skyvis(p,q) ];
        for Jones in joneslist:
          mulops.insert(0,Jones(p));
          mulops.append(Jones(q,'conj') ** Meq.ConjTranspose(Jones(q)));
        # now make multiply node to compute sky visibility
        if len(mulops) > 1:
          ns.visibility(p,q) << Meq.MatrixMultiply(*mulops);
        else:
          ns.visibility(p,q) << Meq.Identity(*mulops);
      skyvis = ns.visibility;

    # form up list of visibility contributions
    terms = [ skyvis ];
    if dec_sky:
      vis = ns.visibility1;
      for p,q in ifrs:
        vis(p,q) << Meq.Add(skyvis(p,q),dec_sky(p,q));
      skyvis = vis;

    ### OMS 24/10/2011: commenting this out. IFR errors will be applied directly to incoming spigots instead
    # now chain up any visibility processors
    # return self._apply_vpm_list(ns,skyvis,ifrs=ifrs);
    return skyvis;

  def apply_visibility_processing (self,ns,vis,ifrs=None):
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
    print("Your script uses the deprecated MeqMaker.make_tree() method. Please change it to use make_predict_tree()");
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

  def correct_uv_data (self,ns,inputs,outputs=None,sky_correct=None,inspect_ifrs=None,flag_jones=True):
    """makes subtrees for correcting the uv data given by 'inputs'.
    If 'outputs' is given, then it will be qualified by a jones label and by stations pairs
    to derive the output nodes. If it is None, then ns.correct(jones_label) is used as a base name.
    By default only uv-Jones corrections are applied, but if 'sky_correct' is set to
    a source object (or source name), then sky-Jones corrections for this particular source
    are also put in.
    If 'flag_jones' is given, then Jones terms will be flagged
    Returns an unqualified node that must be qualified with a station pair to get visibilities.
    """;
    stations = Meow.Context.array.stations();
    ifrs = Meow.Context.array.ifrs();
    inspect_ifrs = inspect_ifrs or ifrs;

    ### OMS 24/10/2011: commenting this out. IFR errors will be applied directly to incoming spigots instead
    ## apply vpm corrections, if any
    #for vpm in self._uv_vpm_list:
      #module = self._get_selected_module(vpm.label,vpm.modules);
      #if module and hasattr(module,'correct_visibilities'):
        #inspectors = [];
        #nodes = inputs(vpm.label);
        #if module.correct_visibilities(nodes,inputs,ns=getattr(ns,vpm.label).Subscope(),
             #tags=vpm.label,label=vpm.label,inspectors=inspectors,meqmaker=self) is not None:
          ## add inspectors to internal list
          #for insp in inspectors:
            #self.add_inspector_nodes(insp);
          #inputs = nodes;

    # now build up a correction chain for every station
    correction_chains = dict([(p,[]) for p in stations]);

    # if overall sky-Jones correction is enabled, collect all sky-Jones terms for which it hasn't been
    # individually disabled
    if sky_correct:
      for jt in self._sky_jones_list:
        # if using coherency decomposition, we will already have defined a "uvjones" node
        # containing a product of all the sky-Jones terms, so use that
        skyjones = ns.skyjones(sky_correct);
        if skyjones(stations[0]).initialized():
          for p in stations:
            correction_chains[p].insert(0,skyjones(p));
        else:
          skip = self.are_advanced_options_enabled(jt.label) and \
                    getattr(self,self._make_attr('skip_correct',jt.label),False);
          if not skip:
            Jj,solvable = self._get_jones_nodes(ns,jt,stations,sources=[sky_correct]);
            if Jj:
              Jj = Jj(sky_correct);
              if Jj(stations[0]).initialized():
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
          # if flagging is in effect, make flaggers
          if flag_jones and jt.flaggable and self.are_advanced_options_enabled(jt.label):
            enable,flagtype,flagmin,flagmax,freqmean = [ getattr(self,self._make_attr(attr,jt.label)) for attr in
                                                  ("flag_jones","flag_jones_type","flag_jones_min","flag_jones_max","flag_jones_freqmean") ];
            Jflag = Jj('flag');
            if enable and (flagmin is not None or flagmax is not None):
              if flagtype == self.JonesTerm.FLAG_ABS or flagtype == self.JonesTerm.FLAG_ABSDIAG:
                Jabs  = Jj('abs');
                for p in stations:
                  StdTrees.make_jones_abs_flagger(Jj(p),flagmin,flagmax,jflag=Jflag(p),jabs=Jabs(p),
                                                  diag=(flagtype==self.JonesTerm.FLAG_ABSDIAG),
                                                  freqmean=freqmean,
                                                  flagmask=MSUtils.FLAGMASK_OUTPUT);
                self.make_per_station_bookmarks(Jabs,"%s-Jones absolute value"%jt.label,stations=stations);
                self.make_per_station_bookmarks(Jflag,"%s-Jones flagged"%jt.label,stations=stations);
                Jj = Jflag;
              elif flagtype == self.JonesTerm.FLAG_NORM:
                Jnorm = Jj('norm');
                for p in stations:
                  StdTrees.make_jones_norm_flagger(Jj(p),flagmin,flagmax,jflag=Jflag(p),jnorm=Jnorm(p),
                                                   freqmean=freqmean,
                                                   flagmask=MSUtils.FLAGMASK_OUTPUT);
                self.make_per_station_bookmarks(Jnorm,"%s-Jones matrix norm"%jt.label,stations=stations);
                self.make_per_station_bookmarks(Jflag,"%s-Jones flagged"%jt.label,stations=stations);
                Jj = Jflag;
          # add to correction chain
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
        Jinv(p) << Meq.MatrixInvert22(jj);
        if p != stations[0]:
          Jtinv(p) << Meq.ConjTranspose(Jinv(p));
    elif correction_chains[stations[0]]:
      for p in stations:
        jj = correction_chains[p][0];
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
  docstring = """<P>Enter space- or comma-separated items. The following constructs are recognized:</P>
                  <P>"<I>name</I>" selects a source by name. The name may include shell-style
                  wildcards such as "A?" and "A*".</P>
                  <P>"all" or "*" selects all sources.</P>
                  <P>"=<i>tagname</i>" selects all sources that have a given tag.</P>
                  <P>"<i>tagname^^value</i>" or "=<i>tagname^^value</i>" selects sources by comparing the given
                  tag to the given (float) value. Here, "^^" represents a comparison operator, and may be any one of
                  = (or ==), !=, &lt;=, &gt;=, &lt;, &gt;, or .eq., .ne., .le., .ge., .lt., .gt. Note also that "value" may be
                  followed by "d", "m" or "s", in which case it is converted from degrees, minutes or seconds into
                  radians.</P>
                  <P>Successive items are treated as logical-OR (i.e. each item adds to the selection), unless the
                  item is prefixed by a modifier (with no space in between):</P>
                  <P>"&<I>selection</I>" means 'and', i.e. "I&ge;1 Q&ge;0 &r&lt;2d" is interpreted as "(I&ge;1 or Q&ge;0) and r&lt;2d".</P>
                  <P>"-<I>selection</I>" means 'except', i.e. "I&ge;1 Q&ge;0 -r&lt;2d" is interpreted as "(I&ge;1 or Q&ge;0) and r&ge;2d". Likewise, "-<i>name</i>" simply deselects sources by name.</P>"""

  def __init__ (self,title,tdloption_namespace=None,doc=None,annotate=True):
    if tdloption_namespace:
      self.tdloption_namespace = tdloption_namespace;
    self.subset_enabled = False;
    global _annotation_label_doc;
    opts = [ TDLOption('source_subset',"Sources",["all"],more=str,namespace=self,doc=self.docstring) ];
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

  # selection predicates (for the tag**value syntax)
  _select_predicates = {
    '=':lambda x,y:x==y,
    '==':lambda x,y:x==y,
    '!=':lambda x,y:x!=y,
    '>=':lambda x,y:x>=y,
    '<=':lambda x,y:x<=y,
    '>' :lambda x,y:x>y,
    '<' :lambda x,y:x<y,
    '.eq.':lambda x,y:x==y,
    '.ne.':lambda x,y:x!=y,
    '.ge.':lambda x,y:x>=y,
    '.le.':lambda x,y:x<=y,
    '.gt.' :lambda x,y:x>y,
    '.lt.' :lambda x,y:x<y
  };
  _units = dict(d=DEG,m=DEG/60.,s=DEG/3600.);

  # regex matching the tag**value[dms] operation
  # the initial "=" is allowed for backwards compatibility with =tag=value constructs
  _re_tagcomp = re.compile("^(?i)=?([^=<>!.]+)(%s)([^dms]+)([dms])?"%"|".join([key.replace('.','\.') for key in list(_select_predicates.keys())]));

  @staticmethod
  def _parse_float (strval):
    try:
      return float(strval);
    except:
      None;

  @staticmethod
  def filter_subset (subset,srclist0,tag_accessor=Meow.SkyComponent.get_attr):
    all = set([src.name for src in srclist0]);
    srcs = set();
    for ispec,spec0 in enumerate(re.split("[\s,]+",subset)):
      spec = spec0.strip();
      # "all" selects all sources
      if spec.lower() == "all":
        srcs = all;
        continue;
      if not spec:
        continue;
      # strip off modifier at beginning
      if spec[0] in "-&|":
        op = spec[0];
        spec = spec[1:];
      else:
        op = "|";
      # if first modifier is AND or EXCEPT, then implictly select all sources first
      if not ispec and op in "&-":
        srcs = all;
      # check for tag**value construct first
      match_tagcomp = SourceSubsetSelector._re_tagcomp.match(spec);
      if match_tagcomp:
        tag,oper,value,unit = match_tagcomp.groups();
        value = SourceSubsetSelector._parse_float(value);
        predicate = SourceSubsetSelector._select_predicates.get(oper.lower());
        scale = SourceSubsetSelector._units.get(unit,1);
        # ignore invalid selections
        if oper is None or value is None:
          print("Warning: invalid source subset selection '%s', ignoring"%spec0);
          continue;
        value *= scale;
        # do the selection
        print(predicate,tag,value);
        srctag = [ (src,tag_accessor(src,tag)) for src in srclist0 ];
        selection = [ src.name for src,tag in srctag if tag is not None and predicate(tag,value) ];
      # then, a =tag construct
      elif spec.startswith("="):
        spec = spec[1:];
        selection = [ src.name for src in srclist0 if src.get_attr(spec) ];
      # everything else treated as a source name (pattern)
      else:
        selection = fnmatch.filter(all,spec);
      # apply this selection to current source set
      if op == "-":
        srcs.difference_update(selection);
      elif op == "&":
        srcs.intersection_update(selection);
      else:
        srcs.update(selection);
      # print stats
      print("applied %s (involving %d sources), %d sources now selected"%(spec0,len(selection),len(srcs)));
    return [ src for src in srclist0 if src.name in srcs ];

  def filter (self,srclist0):
    if not self.subset_enabled or self.source_subset == "all":
      return srclist0;
    srclist = self.filter_subset(self.source_subset,srclist0);
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
  f = open(filename,'wt');
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
      attrkeys = list(attrs.keys());
      from past.builtins import cmp
      from functools import cmp_to_key
      attrkeys.sort(key=cmp_to_key(lambda a,b:cmp(len(b),len(a))));
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
