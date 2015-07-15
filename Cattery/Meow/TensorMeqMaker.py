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
from Meow.MeqMaker import *
import Meow
from Meow import StdTrees,ParmGroup,Parallelization,MSUtils

import itertools

TDLCompileOption("psv_class","Use which node for PSV tensors",["PSVTensor","CUDAPointSourceVisibility","ThrustPointSourceVisibility"]);
#psv_class = "PSVTensor";

class TensorMeqMaker (MeqMaker):
  def __init__ (self,**kw):
    kw['use_decomposition'] = False;
    MeqMaker.__init__(self,**kw);
    
  def smearing_options (self):
    return [ 
      TDLOption('fix_time_smearing',"Fix time interval for smearing calculation",[None],more=float,namespace=self),
      TDLOption('fix_freq_smearing',"Fix bandwidth for smearing calculation",[None],more=float,namespace=self),
      TDLOption('smearing_count',"Apply to N brightest non-point/gaussian sources only",["all",10,100],more=int,namespace=self,
        doc="""<P>Smearing is somewhat expensive to calculate in this implementation, so you may choose to limit it
        to some number of brightest sources (e.g. 50)</P>"""),
    ];


  def _get_skyjones_tensor (self,ns,jt,stations,source_lists,other_sources=[]):
    """Returns the per-source tensor node associated with the given JonesTerm ('jt').
    If the term has been disabled (through compile-time options), returns None.
      'stations' is a list of stations.
      'source_lists' is a list of source groups, each is a tuple of (srclist,lmntensor)
    Returns a tuple of (J1,J1conj),(J2,J2conj),... basenode pairs, one pair per source group.
    If a basenode pair is None,None, then that source group does not have that Jones term.
    Each basenode should be qualified with a station name to get the per-station tensor node.
    If the Jones module is completely disabled, just returns None.

    Note that 'other_sources', if passed in, is added to the list of sources for which Jones
    noes are computed. This is done for compatibility with mixed models (containing both
    tensor and non-tensor visiblities), since compute_jones() is called only once.
    """;
    # get the selected module
    attrname = 'tensor_basenodes_%s'%ns.name if ns.name else 'tensor_basenodes';
    basenodes = getattr(jt,attrname,None);
    if basenodes is None:
      module = self._get_selected_module(jt.label,jt.modules);
      if not module:
        return None;
      # get advanced options: is this Jones term per-source?
      advanced_opt = self.are_advanced_options_enabled(jt.label);
      per_source = getattr(self,self._make_attr('per_source',jt.label),None) if advanced_opt else PER_SOURCE;
      # if module has a tensor method, and we have one Jones term per source, use that
      # (if not per source, then fall back on calling compute_jones() and then assembling a tensor)
      if hasattr(module,"compute_jones_tensor") and per_source == PER_SOURCE:
        basenodes = [];
        # are we doing a per-station Jones term?
        same_all_stations = advanced_opt and getattr(self,self._make_attr('all_stations',jt.label),False);
        real_stations = stations[:1] if same_all_stations else stations;
        # do we have an associated pointing error module?
        dlm = self._get_pointing_error_nodes(ns,jt,stations);
        # loop over source groups and create a tensor for each
        for i,(srclist,lmnT) in enumerate(source_lists):
          # does this only apply to a subset of sources?
          if advanced_opt:
            srclist = jt.subset_selector.filter(srclist);
          if srclist:
            # call module to get the tensor
            inspectors = [];
            Tbase = module.compute_jones_tensor(ns["%sT%d"%(jt.label,i)],
                            srclist,real_stations,lmn=lmnT,pointing_offsets=dlm,inspectors=inspectors);
            if Tbase is None:
              # if no tensor available, fall back to Jones nodes
              basenodes = None;
              break;
            # make inspectors
            if inspectors:
              for node in inspectors:
                self.add_inspector_node(node);
            # make conjugate tensor
            Tbaseconj = ns["%sT%d^H"%(jt.label,i)];
            # has the user specified the same Jones term for all stations, or has the module
            # returned the same tensor for all stations? Make one conjuaget then
            p0 = stations[0];
            if same_all_stations or all([ Tbase(p) is Tbase(p0) for p in stations[1:] ]):
              tbase = Tbase(p0);
              tconj = Tbaseconj << Meq.ConjTranspose(Tbase(p0),tensor=True);
              Tbase = lambda p,T=tbase:T;
              Tbaseconj = lambda p,T=tconj:T;
            # else make conjugate tensor for each station
            else:
              for p in stations[1:]:
                Tbaseconj(p) << Meq.ConjTranspose(Tbase(p),tensor=True);
            basenodes.append((Tbase,Tbaseconj));
          else:
            basenodes.append((None,None));
      # if tensor method was not available, fall back to Jones matrices, and assemble a tensor from that
      if basenodes is None:
        basenodes = [];
        all_sources = list(itertools.chain(other_sources,*[src for src,lmnT in source_lists]));
        Jj,solvable = self._get_jones_nodes(ns,jt,stations,sources=all_sources);
        if Jj is None:
          return None;
        # compose tensors
        for i,(srclist,lmnT) in enumerate(source_lists):
          # get the nodes per each source and station, replace all unitialized  nodes with unity
          nodes = [ [ (Jj(src,p) if Jj(src,p).initialized() else 1) for src in srclist ] for p in stations ];
          # if none of them are initialized, skip this tensor entirely
          if all([all([n is 1 for n in nl]) for nl in nodes ]):
            basenodes.append((None,None));
          else:
            Tbase = ns["%sT%d"%(jt.label,i)];
            Tbaseconj = ns["%sT%dconj"%(jt.label,i)];
            for p in stations:
              # compose base tensor, using unity for uninitialized nodes
              nodes = [ Jj(src,p) for src in srclist ];
              Tbase(p) << Meq.Composer(dims=[0],*[ n if n.initialized() else 1 for n in nodes ]);
              # make conjugate tensor
              if p is not stations[0]:
                Tbaseconj(p) << Meq.ConjTranspose(Tbase(p),tensor=True);
            basenodes.append((Tbase,Tbaseconj));
      setattr(jt,attrname,basenodes);
    # return basenodes
    return basenodes;

  def make_predict_tree (self,ns,sources=None,uvdata=None,ifrs=None):
    """makes predict trees using the sky model and ME.
    'ns' is a node scope
    'sources' is a list of sources; the current sky model is used if None.
    'uvdata' is a basenode for additional model uv-data (which should be qualified with a station pair). If supplied, it will be
      added to the sky model, before applying any uv terms.
    Returns a base node which should be qualified with a station pair.
    """;
    # enable dUVW, if smearing is in effect
    Meow.Context.array.enable_uvw_derivatives(self.use_smearing);
    # setup station list and baselines
    stations = Meow.Context.array.stations();
    ifrs = ifrs or Meow.Context.array.ifrs();
    # use sky model if no source list is supplied
    all_sources = sources if sources is not None else self.get_source_list(ns);

    # split sky model into point/gaussian sources and others\
    # The first can be predicted via tensor nodes
    point_sources = [ src for src in all_sources if type(src) in (Meow.PointSource,Meow.GaussianSource) ];
    other_sources = [ src for src in all_sources if type(src) not in (Meow.PointSource,Meow.GaussianSource) ];

    ### build up the sky-Jones component for point/gaussian sources
    if point_sources:
      ### figure out how to split the source list into partitions with the same applicable tensors
      # initial list has one partition: all sources
      partitions = [ set([ src.name for src in point_sources ]) ];
      for jt in self._sky_jones_list:
        # if Jones term is enabled, and advanced options are enabled, and a subset
        # is specified, recalcaulte partitions, else no change
        if self.is_group_enabled(jt.label) and self.are_advanced_options_enabled(jt.label):
          subset = set([src.name for src in jt.subset_selector.filter(point_sources)]);
          if not subset:
            continue;
          new_partitions = [];
          # see what the overlap is with existing partitions
          for part in partitions:
            overlap = part&subset;
            # the overlap and the remainder of this partition, whatever is not empty,
            # become new parititons
            new_partitions += [ p1 for p1 in overlap,part-overlap if p1 ];
            # remove overlap from subset
            subset -= overlap;
          partitions = new_partitions;
      # Convert partitions to lists of sources
      sgroups = [
        [ src for src in point_sources if src.name in part ] for part in partitions
      ];

      ## create lmn tensor per each source group
      source_groups = [];
      for igrp,sources in enumerate(sgroups):
        lmn_static = [ src.direction.lmn_static() for src in sources ];
        lmnT = ns["lmnT%d"%igrp];
        # if all sources have static LMN coordinates, use a single constant node
        if all([ lmn is not None for lmn in lmn_static ]):
          lmnT << Meq.Constant(lmn_static);
        # else compose a tensor
        else:
          lmnT << Meq.Composer(dims=[0],*[src.direction.lmn() for src in sources]);
        # add to group list
        source_groups.append((sources,lmnT));

      ### collect sky-Jones tensors per each source group
      jones_tensors = [ list() for grp in source_groups ];
      for jt in self._sky_jones_list:
        # get the tensor basenodes
        basenodes = self._get_skyjones_tensor(ns,jt,stations,source_groups,other_sources=other_sources);
        if basenodes is not None:
          # now for every source grouping, add those basenodes that are not None to the jones_tensor list for that group
          for igrp,grp in enumerate(source_groups):
            if basenodes[igrp][0] is not None:
              jones_tensors[igrp].append(basenodes[igrp]);

      uvw = Meow.Context.array.uvw_ifr();
      ### now make the PSV tensor for esch source
      psvts = [];
      ns.null_shape << Meq.Constant([0,0,0]);
      psv_kwargs = dict();
      if self.fix_time_smearing is not None:
        psv_kwargs['fixed_time_smearing_interval'] = self.fix_time_smearing;
      if self.fix_freq_smearing is not None:
        psv_kwargs['fixed_freq_smearing_interval'] = self.fix_freq_smearing;
      for igrp,(sources,lmnT) in enumerate(source_groups):
        ### create brightness tensor
        any_pol = any([src.is_polarized() for src in sources])
        # see if we have constant brightnesses
        B_static = [ src.brightness_static() for src in sources ];
        if all([ b is not None for b in B_static ]):
          BT = ns["BT%d"%igrp] << Meq.Constant(B_static);
        else:
          BT = ns["BT%d"%igrp] << Meq.Composer(dims=[0],*[ src.brightness(always_matrix=any_pol) for src in sources ]);
        ### create shape tensor 
        shape_static = [ src.shape_static() if hasattr(src,'shape_static') else [0,0,0] for src in sources ];
        if all([ s is not None for s in shape_static ]):
          shapeT = ns["shapeT%d"%igrp] << Meq.Constant(shape_static);
        else:
          shapeT = ns["shapeT%d"%igrp] << Meq.Composer(dims=[0],*[Meq.Composer(*src.shape()) if hasattr(src,'shape') else ns.null_shape for src in sources]);
        ### now make the predict nodes
        psvt = ns["psvT%d"%igrp];
        psvts.append(psvt);
        for p,q in ifrs:
          jt = [];
          for Et,Etconj in jones_tensors[igrp]:
            jt += [ Et(p),Etconj(q) ];
          v = psvt(p,q) << Meq[psv_class](lmnT,BT,uvw(p,q),shapeT,*jt,**psv_kwargs);
      ### if uvdata is specified, add it to the summation terms
      if uvdata is not None:
        psvts.append(uvdata);
      ### generate output visibilities, as a sum of the terms
      if len(psvts) == 1:
        uvdata = psvts[0];
      else:
        uvdata = ns.psvTsum;
        for p,q in ifrs:
          uvdata(p,q) << Meq.Add(*[x(p,q) for x in psvts]);

    ### OK, at this stage the point source contributions (plus original uvdata, if supplied) are in uvdata.
    ### Invoke the base function to build trees for any other sources in the model, and to apply uv-plane terms.
    if uvdata is not None or other_sources:
      return MeqMaker.make_predict_tree(self,ns,other_sources,uvdata=uvdata,ifrs=ifrs);
    else:
      # no other sources and no data: return nulls
      for p,q in ifrs:
        ns.predict(p,q) << 0;
      return ns.predict;


