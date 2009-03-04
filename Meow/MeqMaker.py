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
import sets
import inspect

import Meow
from Meow import StdTrees
from Meow import ParmGroup
from Meow import Parallelization

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

class MeqMaker (object):
  def __init__ (self,namespace='me',solvable=False,use_decomposition=None,use_jones_inspectors=None):
    self.tdloption_namespace = namespace;
    self._uv_jones_list = [];
    self._sky_jones_list = [];
    self._uv_vpm_list = [];
    self._sky_vpm_list = [];
    self._compile_options = [];
    self._runtime_options = None;
    self._sky_models = None;
    self._source_list = None;
    self.export_kvis = False;
    self._solvable = solvable;
    self._inspectors = [];
    
    other_opt = [];
    if use_decomposition is None:
      self.use_decomposition = True;
      self.use_decomposition_opt = \
        TDLOption('use_decomposition',"Use source coherency decomposition, if available",True,namespace=self,
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
      
    other_opt.append(
      TDLMenu("Include time & bandwidth smearing",
        TDLOption('smearing_count',"Apply to N brightest sources only",["all",10,100],more=int,namespace=self),
        namespace=self,toggle='use_smearing',default=True,
      )
    );
    other_opt += Meow.IfrArray.compile_options();
    self._compile_options.append(TDLMenu("Measurement Equation options",*other_opt));

  def compile_options (self):
    return self._compile_options;

  def add_sky_models (self,modules):
    if self._sky_models:
      raise RuntimeError,"add_sky_models() may only be called once";
    self._compile_options.append(
        self._module_selector("Sky model","sky",modules,use_toggle=False));
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
                TDLOption("export_karma_incremental","Export each N sources as a separate file",
                          ["all",10,50],more=int,namespace=self)
              ))
    );

  def add_uv_jones (self,label,name,modules,pointing=None):
    return self._add_jones_modules(label,name,True,pointing,modules);
  
  def add_vis_proc_module (self,label,name,modules):
    return self._add_vpm_modules(label,name,True,modules);

  def add_sky_jones (self,label,name,modules,pointing=None):
    return self._add_jones_modules(label,name,False,pointing,modules);

  class VPMTerm (object):
    """Essentially a record representing one VPM in an ME.
    A VPM has the following fields:
      label:      a label (e.g. "G", "E", etc.)
      name:       a descriprive name
      modules:    a list of possible modules implementing the Jones set
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
      name:       a descriprive name
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

  def _add_vpm_modules (self,label,name,is_uvplane,modules):
    if not modules:
      raise RuntimeError,"No modules specified for %s"%name;
    if not isinstance(modules,(list,tuple)):
      modules = [ modules ];
    # make option menus for selecting a visiblity processor
    mainmenu = self._module_selector("Use %s"%name,label,modules);
    self._compile_options.append(mainmenu);
    # Add to internal list
    term = self.VPMTerm(label,name,modules);
    if is_uvplane:
      self._uv_vpm_list.append(term);
    else:
      self._sky_jones_list.append(term);

  def _add_jones_modules (self,label,name,is_uvplane,pointing,modules):
    if not modules:
      raise RuntimeError,"No modules specified for %s Jones"%label;
    if not isinstance(modules,(list,tuple)):
      modules = [ modules ];
    # extra options for pointing
    if pointing:
      if not isinstance(pointing,(list,tuple)):
        pointing = [ pointing ];
      pointing_menu = [ self._module_selector("Apply pointing errors to %s"%label,
                                              label+"pe",pointing,nonexclusive=True) ];
    else:
      pointing_menu = [];
    # make option menus for selecting a jones module
    mainmenu = self._module_selector("Use %s Jones (%s)"%(label,name),
                                     label,modules,extra_opts=pointing_menu);
    self._compile_options.append(mainmenu);
    # Add to internal list
    # each jones module set is represented by a list of:
    #   'label','full name',module_list,pointing_module_list,base_node
    # base_node is inintially None; when a module is ultimately invoked to
    # build the trees, the base node gets stashed here for later reuse
    if is_uvplane:
      self._uv_jones_list.append(self.JonesTerm(label,name,modules));
    else:
      self._sky_jones_list.append(self.SkyJonesTerm(label,name,modules,pointing));

  def runtime_options (self,nest=True):
    if self._runtime_options is None:
      self._runtime_options = [];
      # build list of currently selected modules
      mods = [];
      mods.append((self._get_selected_module('sky',self._sky_models),'Sky model options',None,None));
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
    return self._runtime_options;

  def _module_togglename (self,label,mod):
    return "%s_%s"%(label,_modname(mod));

  def _module_selector (self,menutext,label,modules,extra_opts=[],use_toggle=True,**kw):
    # Forms up a "module selector" submenu for the given module set.
    # for each module, we either take the options returned by module.compile_options(),
    # or pass in the module itself to let the menu "suck in" its options
    toggle = self._make_attr(label,"enable");
    exclusive = self._make_attr(label,"module");
    if not use_toggle:
      setattr(self,toggle,True);
      toggle = None;
    if len(modules) == 1:
      modname = _modname(modules[0]);
      mainmenu = TDLMenu(menutext,toggle=toggle,namespace=self,name=modname,
                         *(_modopts(modules[0],'compile')+list(extra_opts)),**kw);
      setattr(self,exclusive,modname);
    else:
      # note that toggle symbol is set to modname. The symbol is not really used,
      # since we make an exclusive parent menu, and the symbol will be assigned to
      # its option value
      submenus = [ TDLMenu("Use '%s' module"%_modname(mod),name=_modname(mod),
                            toggle=_modname(mod).replace('.','_'),namespace={},
                            *_modopts(mod,'compile'))
                    for mod in modules ];
      mainmenu = TDLMenu(menutext,toggle=toggle,exclusive=exclusive,namespace=self,
                          *(submenus+list(extra_opts)),**kw);
    return mainmenu;

  def _make_attr (*comps):
    """Forms up the name of an 'enabled' attribute for a given module group (e.g. a Jones term)"""
    return "_".join(comps);
  _make_attr = staticmethod(_make_attr);

  def is_group_enabled (self,label):
    """Returns true if the given module is enabled for the given Jones label"""
    return getattr(self,self._make_attr(label,"enable"));

  def get_inspectors (self):
    """Returns list of inspector nodes created by this MeqMaker""";
    return self._inspectors or None;

  def _get_selected_module (self,label,modules):
    # check global toggle for this module group
    if self.is_group_enabled(label):
      # figure out which module we'll be using. If only one is available, use it regardless.
      if len(modules) == 1:
        return modules[0];
      else:
        selname = getattr(self,self._make_attr(label,"module"));
        for mod in modules:
          if _modname(mod).replace('.','_') == selname:
            return mod;
    return None;

  def estimate_image_size (self):
    module = self._get_selected_module('sky',self._sky_models);
    if module:
      estimate = getattr(module,'estimate_image_size',None);
      if callable(estimate):
        return estimate();
    return None;

  def get_source_list (self,ns):
    if self._source_list is None:
      module = self._get_selected_module('sky',self._sky_models);
      if not module:
        raise RuntimeError,"No source list supplied and no sky model set up";
      self._source_list = module.source_list(ns);
      # add indices to source attrubutes (used by the %# directive when exporting  annotations)
      for isrc,src in enumerate(self._source_list):
        src.set_attr("#",isrc);
    return self._source_list;
  
  def _add_inspector (self,inspector_node,name=None):
    """adds an inspector node to internal list, and creates a bookmark page for it.""";
    self._inspectors.append(inspector_node);
    if not name:
      name = inspector_node.name.replace('_',' ');
    Meow.Bookmarks.Page(name).add(inspector_node,viewer="Collections Plotter");

  def _get_jones_nodes (self,ns,jt,stations,sources=None,solvable_sources=sets.Set()):
    """Returns the Jones nodes associated with the given JonesTerm ('jt'). If
    the term has been disabled (through compile-time options), returns None.
    'stations' is a list of stations.
    'sources' is a list of source (for sky-Jones only).
    Return value is tuple of (basenode,solvable). Basenode should be qualified
    with source (if sky-Jones) and station to obtain the Jones term. Solvable is True
    if the Jones set is to be treated as potentially solvable.
    If basenode==None, the given JonesTerm is not implemented and may be omitted from the ME.
    """;
    if jt.base_node is None:
      module = self._get_selected_module(jt.label,jt.modules);
      if not module:
        return None,False;
      prev_num_solvejobs = ParmGroup.num_solvejobs();  # to keep track of whether module creates its own
      jones_inspectors = [];
      jt.solvable = False;
      # For sky-Jones terms, see if this module has pointing offsets enabled
      if isinstance(jt,self.SkyJonesTerm):
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
              jones_inspectors += inspectors;
            elif dlm:
              jones_inspectors.append(
                  ns.inspector(jt.label)('dlm') << StdTrees.define_inspector(dlm,stations));
      else:
        dlm = None;
      # now make the appropriate matrices
      # note that Jones terms are computed using the original source list.
      # this is to keep the extra corruption qualifiers from creeping into
      # their names
      Jj = ns[jt.label];    # create name hint for nodes
      inspectors = [];
      jt.base_node = Jj = module.compute_jones(Jj,sources=sources,stations=stations,
                                          pointing_offsets=dlm,solvable_sources=solvable_sources,
                                          tags=jt.label,label=jt.label,
                                          meqmaker=self,
                                          inspectors=inspectors);
      # ignoring the pointing-solvable attribute for now
      jt.solvable = getattr(module,'solvable',True);
      # if module does not make its own inspectors, add automatic ones
      if Jj and self.use_jones_inspectors:
        if inspectors:
          jones_inspectors += inspectors;
        elif Jj:
          qual_list = [stations];
          if sources:
            qual_list.insert(0,[src for src in sources if Jj(src,stations[0]).initialized()]);
          jones_inspectors.append(
              ns.inspector(jt.label) << StdTrees.define_inspector(Jj,*qual_list));
        # add inspectors to internal list
        for insp in jones_inspectors:
          self._add_inspector(insp);
      ## NB: this is too messy, let's always use solvejobs instead
      ## see if module has created any solvejobs; create one automatically if not
      #if self._solvable:
        #if ParmGroup.num_solvejobs() == prev_num_solvejobs:
          #parms = jt.base_node.search(tags="solvable");
          #if parms:
            #pg = ParmGroup.ParmGroup(jt.label,parms);
            #ParmGroup.SolveJob("solve_%s"%jt.label,"Solve for %s"%jt.label,pg);
    return jt.base_node,jt.solvable;

  def make_predict_tree (self,ns,sources=None):
    """makes predict trees using the sky model and ME.
    'ns' is a node scope
    'sources' is a list of sources; the current sky model is used if None.
    Returns a base node which should be qualified with a station pair.
    """;
    stations = Meow.Context.array.stations();
    ifrs = Meow.Context.array.ifrs();
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
        if not sources:
          return self._apply_vpm_list(ns,dec_sky);
    
    # Now, proceed to build normal trees for non-decomposable sources.
    # An important optimization when solving for sky terms is to put all the sources containing
    # solvable corruptions in one patch, and the rest in another 
    # (this makes optimal use of cache when solving for the corruption).
    # This list will contain a list of (name,src) tuples, and will be updated as each sky-Jones term is applied
    sourcelist = [ (src.name,src) for src in sources ];
    # This set will contain names of sources to which solvable sky-Jones were been applied
    solvable_skyjones = sets.Set();

    # apply all sky Jones terms
    for jt in self._sky_jones_list:
      Jj,solvable = self._get_jones_nodes(ns,jt,stations,sources=sources,solvable_sources=solvable_skyjones);
      # if this Jones is enabled (Jj not None), corrupt each source
      if Jj:
        newlist = [];
        # go through sourcs, and make new list of corrupted sources.
        # note that not every source is necessarily corrupted
        # print jt.label,solvable,len(solvable_skyjones);
        for name,src in sourcelist:
          jones = Jj(name);
          # if Jones term is initialized, corrupt and append to corr_sources
          # if not initialized, then append to uncorr_sources 
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

    # now form up patch
    if dec_sky:
      patchname = 'sky2';
    else:
      patchname = 'sky';
    allsky = Meow.Patch(ns,patchname,Meow.Context.observation.phase_centre,components=sky_sources);

    # add uv-plane effects
    for jt in self._uv_jones_list:
      Jj,solvable = self._get_jones_nodes(ns,jt,stations);
      if Jj:
        allsky = allsky.corrupt(Jj);

    # now, if we also have a contribution from decomposable sources, add it here
    if dec_sky:
      vis = ns.visibility('sky');
      vis2 = allsky.visibilities();
      for p,q in ifrs:
        vis(p,q) << dec_sky(p,q) + vis2(p,q);
    else:
      vis = allsky.visibilities();
    
    # now chain up any visibility processors
    return self._apply_vpm_list(ns,vis);
          
  def _apply_vpm_list (self,ns,vis):
    # chains up any visibility processors, and applies them to the visibilities.
    # Returns new visibilities.
    for vpm in self._uv_vpm_list:
      module = self._get_selected_module(vpm.label,vpm.modules);
      if module:
        inspectors = [];
        nodes = vis(vpm.label);
        if module.process_visibilities(nodes,vis,ns=getattr(ns,vpm.label).Subscope(),
             tags=vpm.label,label=vpm.label,inspectors=inspectors) is not None:
          # add inspectors to internal list
          for insp in inspectors:
            self._add_inspector(insp);
          vis = nodes;
    return vis;

  make_tree = make_predict_tree; # alias for compatibility with older code

  def correct_uv_data (self,ns,inputs,outputs=None,sky_correct=None,inspect_ifrs=None):
    """makes subtrees for correcting the uv data given by 'inputs'.
    If 'outputs' is given, then it will be qualified by a jones label and by stations pairs
    to derive the output nodes. If it is None, then ns.correct(jones_label) is used as a base name.
    By default only uv-Jones corrections are applied, but if 'sky_correct' is set to
    a source object (or source name), then sky-Jones corrections for this particular source
    are also put in.
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
             tags=vpm.label,label=vpm.label,inspectors=inspectors) is not None:
          # add inspectors to internal list
          for insp in inspectors:
            self._add_inspector(insp);
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
        Jprod(p) << Meq.MatrixMultiply(*correction_chains[p]);
        Jinv(p) << Meq.MatrixInvert22(Jprod(p));
        if p != stations[0]:
          Jtinv(p) << Meq.ConjTranspose(Jinv(p));
    elif correction_chains[stations[0]]:
      for p in stations:
        Jinv(p) << Meq.MatrixInvert22(correction_chains[p][0]);
        if p != stations[0]:
          Jtinv(p) << Meq.ConjTranspose(Jinv(p));
    else:
      for p,q in ifrs:
        outputs(p,q) << Meq.Identity(inputs(p,q));
      # make an inspector for the results
      StdTrees.vis_inspector(ns.inspector('output'),outputs,ifrs=inspect_ifrs,bookmark=False);
      self._add_inspector(ns.inspector('output'),name='Inspect corrected data or residuals');
      return outputs;
    # now apply the correction matrices
    StdTrees.vis_inspector(ns.inspector('uncorr'),inputs,ifrs=inspect_ifrs,bookmark=False);
    self._add_inspector(ns.inspector('uncorr'),name='Inspect uncorrected data/residuals');
    for p,q in ifrs:
      outputs(p,q) << Meq.MatrixMultiply(Jinv(p),inputs(p,q),Jtinv(q));

    # make an inspector for the results
    StdTrees.vis_inspector(ns.inspector('output'),outputs,ifrs=inspect_ifrs,bookmark=False);
    self._add_inspector(ns.inspector('output'),name='Inspect corrected data/residuals');

    return outputs;
    
  def close (self):
    if self.export_kvis and self._source_list and self.export_karma_filename:
      if self.export_karma_incremental and self.export_karma_incremental != "all":
        filename,ext = os.path.splitext(self.export_karma_filename);
        for i0 in range(0,len(self._source_list),self.export_karma_incremental):
          i1 = min(i0+self.export_karma_incremental,len(self._source_list));
          export_karma_annotations(self._source_list[i0:i1],
            filename="%s-%d-%d%s"%(filename,i0,i1-1,ext),
            label_format=self.export_karma_label);
      else:
        export_karma_annotations(self._source_list,
            self.export_karma_filename,label_format=self.export_karma_label);
    
class SourceSubsetSelector (object):
  def __init__ (self,title,tdloption_namespace=None,doc=None):
    if tdloption_namespace:
      self.tdloption_namespace = tdloption_namespace;
    self.subset_enabled = False;
    global _annotation_label_doc;
    subset_opt = TDLMenu(title,toggle='subset_enabled',namespace=self,doc=doc,
      *( 
        TDLOption('source_subset',"Sources",["all"],more=str,namespace=self,
                  doc="""Enter source names separated by space."""),
        TDLMenu("Annotate selected sources",doc="""If you're exporting annotations for your sky model,
    you may choose to mark the selected sources using different colours and/or labels""",
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
        )
      ));
    self._src_set = None;
    self.options = [ subset_opt ];

  def filter (self,srclist):
    if not self.subset_enabled:
      return [];
    if not self._src_set:
      self._src_set = sets.Set(self.source_subset.split(" "));
    # make list
    srclist = [ src for src in srclist if src.name in self._src_set ];
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
      lbl_color='blue'):
  """Exports a list of sources as a Karma annotations file. 
  label_format is used to generate source labels.\n"""+_annotation_label_doc;
  f = file(filename,'wt');
  f.write('COORD W\nPA STANDARD\nCOLOR GREEN\nFONT hershey12\n');
  # calculate default size for crosses
  xcross=0.01;
  for src in sources:
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
      f.write('ELLIPSE %.12f %.12f %f %f %f\n'%(ra_d,dec_d,sx/math.pi*180.0,sy/math.pi*180.0,pa/math.pi*180.0));
    # else assume point source and do a cross
    else:
      f.write('CROSS %.12f %.12f %f %f\n'%(ra_d,dec_d,xcross,xcross));
    # now generate label
    label = src.get_attr('ANNOTATION_LABEL',label_format);
    if label:
      attrs = dict(src.attrs);
      # populate attribute dict
      attrs['%'] = '%';
      attrs['N'] = src.name;
      for st in "IQUV":
        attrs[st] = src.get_value(st) or 0;
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
