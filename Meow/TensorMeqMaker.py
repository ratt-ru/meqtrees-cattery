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

class TensorMeqMaker (MeqMaker):
  def __init__ (self,**kw):
    kw['use_decomposition'] = False;
    MeqMaker.__init__(self,**kw);

  def _get_skyjones_tensor (self,ns,jt,stations,source_lists,other_sources=[]):
    """Returns the per-source tensor node associated with the given JonesTerm ('jt'). 
    If the term has been disabled (through compile-time options), returns None.
      'stations' is a list of stations.
      'source_lists' is a list of source groups, each is a tuple of (srclist,lmntensor)
    Returns a tuple of (J1,J1conj),(J2,J2conj),... basenode pairs, one pair per source group.
    If a basenode pair is None,None, then that source group does not have that Jons term.
    Each basenode should be qualified with a station name to get the per-station tensor node. 
    If the Jones module is completely disabled, just returns None.
    
    Note that 'other_sources', if passed in, is added to the list of sources for which Jones
    noes are computed. This is done for compatibility with mixed models (containing both
    tensor and non-tensor visiblities), since compute_jones() is called only once.
    """;
    # get the selected module
    basenodes = getattr(jt,'tensor_basenodes',None);
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
            Tbase = ns["%sT%d"%(jt.label,i)];
            if module.compute_jones_tensor(Tbase,srclist,real_stations,lmn=lmnT,pointing_offsets=dlm) is None:
              # if no tensor available, fall back to Jones nodes
              basenodes = None;
              break;
            # make conjugate tensor
            Tbaseconj = ns["%sT%dconj"%(jt.label,i)];
            # propagate this to all stations
            if same_all_stations:
              tbase = Tbase(stations[0]);
              tconj = Tbaseconj(stations[0]) << Meq.ConjTranspose(Tbase(stations[0]),tensor=True);
              Tbase = lambda p:tbase;
              Tbaseconj = lambda p:tconj;
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
      jt.tensor_basenodes = basenodes;
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
    
    # split sky model into point sources and others
    point_sources = [ src for src in all_sources if isinstance(src,Meow.PointSource) ];
    other_sources = [ src for src in all_sources if not isinstance(src,Meow.PointSource) ];
    
    ### build up the sky-Jones component for point sources
    if point_sources:
      ### figure out if our source list needs to be partitioned into two subsets. 
      # We do this when e.g. one particular sky-Jones term is applied to only a subset of sources, or
      # if multiple sky-Jones terms have subsets, but all these subsets are the same
      fixed_subset = None;
      for jt in self._sky_jones_list:
        # if Jones term is enabled, and advanced options are enabled
        if self.is_group_enabled(jt.label) and self.are_advanced_options_enabled(jt.label):
          subset = set([src.name for src in jt.subset_selector.filter(point_sources)]);
          if len(subset) and len(subset) < len(point_sources):
            if fixed_subset is None:
              fixed_subset = subset;
            elif fixed_subset != subset:
              fixed_subset = None;
              break;
      # If fixed_subset is None, then we either have different subsets per Jones term, or no subsets at all.
      # Either way, this precludes partitioning
      if fixed_subset is None:
        sgroups = [ point_sources ];
      # else partition sources into two groups
      else:
        sgroups = [ 
          [ src for src in point_sources if src.name in fixed_subset ],
          [ src for src in point_sources if src.name not in fixed_subset ]
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
          lmnT << Meq.Composer(dims=[0],*[src.direction.lmn() for src in source_list]);
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
      for igrp,(sources,lmnT) in enumerate(source_groups):
        ### create brightness tensor
        # see if we have constant brightnesses
        B_static = [ src.brightness_static() for src in sources ];
        if all([ b is not None for b in B_static ]):
          BT = ns["BT%d"%igrp] << Meq.Constant(B_static);
        else:
          BT = ns["BT%d"%igrp] << Meq.Composer(dims=[0],*[src.brightness() for src in sources]);
        ### now make the predict nodes
        psvt = ns["psvT%d"%igrp];
        psvts.append(psvt);
        for p,q in ifrs:
          jt = [];
          for Et,Etconj in jones_tensors[igrp]:
            jt += [ Et(p),Etconj(q) ];
          v = psvt(p,q) << Meq.PSVTensor(lmnT,BT,uvw(p,q),*jt);
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

