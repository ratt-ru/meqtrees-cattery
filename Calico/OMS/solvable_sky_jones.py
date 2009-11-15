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
"""This implements generic solvable sky-Jones modules.
  DiagAmplPhase is a diagonal matrix, with separate solvable groups for
      amplitudes and phases
  FullRealImag is a full 2x2 matric, with separate solvable groups for
      diagonal and off-diagonal terms.
DiagAmplPhase is currently untested. In any case the solver shows poor convergence
for this type of matrix decomposition, so FullRealImag is preferred.
""";

from Timba.TDL import *
import Meow
from Meow import Context
from Meow import Jones,ParmGroup,Bookmarks
from Meow.Parameterization import resolve_parameter
import Meow.MeqMaker 

class DiagAmplPhase (object):
  def __init__ (self,label):
    self.tdloption_namespace = label+".diagamplphase";
    subset_opt = TDLOption('subset',"Apply this Jones term to a subset of sources",
        ["all"],more=str,namespace=self,doc="""Selects a subset of sources to which this 
        Jones term is applied. Enter soure names separated by space""");
    self.options = [ subset_opt ];

  def runtime_options (self):
    return self.options;

##  def compute_jones (self,nodes,stations=None,tags=None,label='',**kw):
##    stations = stations or Context.array.stations();
##    g_ampl_def = Meow.Parm(1);
##    g_phase_def = Meow.Parm(0);
##    nodes = Jones.gain_ap_matrix(nodes,g_ampl_def,g_phase_def,tags=tags,series=stations);

##    # make parmgroups for phases and gains
##    self.pg_phase = ParmGroup.ParmGroup(label+"_phase",
##                    nodes.search(tags="solvable phase"),
##                    table_name="%s_phase.fmep"%label,bookmark=4);
##    self.pg_ampl  = ParmGroup.ParmGroup(label+"_ampl",
##                    nodes.search(tags="solvable ampl"),
##                    table_name="%s_ampl.fmep"%label,bookmark=4);

##    # make solvejobs
##    ParmGroup.SolveJob("cal_"+label+"_phase","Calibrate %s phases"%label,self.pg_phase);
##    ParmGroup.SolveJob("cal_"+label+"_ampl","Calibrate %s amplitudes"%label,self.pg_ampl);

  def compute_jones (self,nodes,sources,stations=None,tags=None,label='',**kw):
    stations = stations or Context.array.stations();
    # figure out which sources to apply to
    if self.subset != "all":
      srcset = set(self.subset.split(" "));
      sources = [ src for src in sources if src.name in srcset ];
    if not sources:
      return None;

    g_ampl_def = Meow.Parm(1);
    g_phase_def = Meow.Parm(0);
    # loop over sources
    #print "tags",tags
    subgroups_phase = [];
    subgroups_ampl = [];
    for src in sources:
      #print  srci,src,self.solve_source,(srci==self.solve_source)
      jones = Jones.gain_ap_matrix(nodes(src),g_ampl_def,g_phase_def,tags=tags,series=stations);
      # add subgroup for this source
      subgroups_phase.append(ParmGroup.Subgroup(src.name,nodes(src).search(tags="solvable phase")));
      subgroups_ampl.append(ParmGroup.Subgroup(src.name,nodes(src).search(tags="solvable ampl")));

    # make parmgroups for phases and gains
    self.pg_phase = ParmGroup.ParmGroup(label+"_phase",
                    nodes.search(tags="solvable phase"),
                    subgroups = subgroups_phase,                    
                    table_name="%s_phase.mep"%label,bookmark=4);
    self.pg_ampl  = ParmGroup.ParmGroup(label+"_ampl",
                    nodes.search(tags="solvable ampl"),
                    subgroups = subgroups_ampl,                 
                    table_name="%s_ampl.mep"%label,bookmark=4);

    ParmGroup.SolveJob("cal_"+label+"_phase","Calibrate %s phases"%label,self.pg_phase);
    ParmGroup.SolveJob("cal_"+label+"_ampl","Calibrate %s amplitudes"%label,self.pg_ampl);


    return nodes;

class FullRealImag (object):
  def __init__ (self,label):
    self.tdloption_namespace = label+".fullrealimag";
    self.subset_selector = Meow.MeqMaker.SourceSubsetSelector("Apply this Jones term to a subset of sources",
                            tdloption_namespace=self.tdloption_namespace);
    self.options = self.subset_selector.options;

  def compile_options (self):
    return self.options;

  def compute_jones (self,jones,sources,stations=None,tags=None,label='',**kw):
    stations = stations or Context.array.stations();
    # figure out which sources to apply to
    sources = self.subset_selector.filter(sources);
    # set up qualifier labels depending on polarization definition
    if Context.observation.circular():
      x,y,X,Y = 'R','L','R','L';
    else:
      x,y,X,Y = 'x','y','X','Y';
    xx,xy,yx,yy = x+x,x+y,y+x,y+y;
    rxx,rxy,ryx,ryy = [ q+":r" for q in xx,xy,yx,yy ];
    ixx,ixy,iyx,iyy = [ q+":i" for q in xx,xy,yx,yy ];
    # create parm definitions for each jones element
    tags = NodeTags(tags) + "solvable";
    diag_real = Meq.Parm(1,tags=tags+"diag real");
    diag_imag = Meq.Parm(0,tags=tags+"diag imag");
    offdiag_real = Meq.Parm(0,tags=tags+"offdiag real");
    offdiag_imag = Meq.Parm(0,tags=tags+"offdiag imag");
    # loop over sources
    parms_diag = [];
    parms_offdiag = [];
    subgroups_diag = [];
    subgroups_offdiag = [];
    for src in sources:
      # now loop to create nodes
      sgdiag = [];
      sgoff = [];
      for p in stations:
        jj = jones(src,p);
        jj << Meq.Matrix22(
          jj(xx) << Meq.ToComplex(
              jj(rxx) << diag_real,
              jj(ixx) << diag_imag
          ),
          jj(xy) << Meq.ToComplex(
              jj(rxy) << offdiag_real,
              jj(ixy) << offdiag_imag
          ),
          jj(yx) << Meq.ToComplex(
              jj(ryx) << offdiag_real,
              jj(iyx) << offdiag_imag
          ),
          jj(yy) << Meq.ToComplex(
              jj(ryy) << diag_real,
              jj(iyy) << diag_imag
          )
        );
        # add to list of parms for groups and subgroups
        parms = [ jj(zz) for zz in rxx,ixx,ryy,iyy ]
        parms_diag += parms;
        sgdiag += parms;
        parms = [ jj(zz) for zz in rxy,ixy,ryx,iyx ];
        parms_offdiag += parms;
        sgoff  += parms;
      # add subgroup for this source
      subgroups_diag.append(ParmGroup.Subgroup(src.name,sgdiag));
      subgroups_offdiag.append(ParmGroup.Subgroup(src.name,sgoff));
    # re-sort by name
    subgroups_diag.sort(lambda a,b:cmp(a.name,b.name));
    subgroups_offdiag.sort(lambda a,b:cmp(a.name,b.name));
    
    # make parmgroups for diagonal and off-diagonal terms
    self.pg_diag  = ParmGroup.ParmGroup(label+"_diag",parms_diag,
                      subgroups=subgroups_diag,
                      table_name="%s_diag.fmep"%label,bookmark=False);
    self.pg_offdiag  = ParmGroup.ParmGroup(label+"_offdiag",parms_offdiag,
                      subgroups=subgroups_offdiag,
                      table_name="%s_offdiag.fmep"%label,bookmark=False);

    # make bookmarks
    Bookmarks.make_node_folder("%s diagonal terms"%label,
      [ jones(src,p,zz) for src in sources  
        for p in stations for zz in "xx","yy" ],sorted=True);
    Bookmarks.make_node_folder("%s off-diagonal terms"%label,
      [ jones(src,p,zz) for src in sources  
        for p in stations for zz in "xy","yx" ],sorted=True);

    # make solvejobs
    ParmGroup.SolveJob("cal_"+label+"_diag","Calibrate %s diagonal terms"%label,self.pg_diag);
    ParmGroup.SolveJob("cal_"+label+"_offdiag","Calibrate %s off-diagonal terms"%label,self.pg_offdiag);

    return jones;

      
class IntrinsicFR (object):
  def __init__ (self,label):
    self.tdloption_namespace = label+".intrinsic_fr";
    subset_opt = TDLOption('subset',"Apply this Jones term to a subset of sources",
        [None],more=str,namespace=self,doc="""Selects a subset of sources to which this 
        Jones term is applied. 'None' applies to all sources.
        You may specify individual indices (0-based) separated by commas or spaces, or ranges, e.g. "M:N" (M to N inclusive), or ":M" (0 to M), or "N:" (N to last).
        Example subset: ":3 5 8 10:12 16:".""");
    self._subset_parser = Meow.OptionTools.ListOptionParser(minval=0,name="sources");
    subset_opt.set_validator(self._subset_parser.validator);
    self.options = [ subset_opt ];
    self.subset = None;

  def compile_options (self):
    return self.options;

  def compute_jones (self,jones,sources,stations=None,tags=None,label='',**kw):
    stations = stations or Context.array.stations();
    # figure out which sources to apply to
    if self.subset:
      srclist = self._subset_parser.parse_list(self.subset);
      sources = [ sources[i] for i in srclist ];
    # create parm definitions for each jones element
    tags = NodeTags(tags) + "solvable";
    rm_parm = Meq.Parm(0,tags=tags+"rm");
    # loop over sources
    for src in sources:
      jj = jones(src);
      jj("RM") << rm_parm;
      fr = jj("fr") << jj("RM") / Meq.Sqr(Meq.Freq());
      cos_fr = jj("cos") << Meq.Cos(fr);
      sin_fr = jj("sin") << Meq.Sin(fr);
      jj << Meq.Matrix22(cos_fr,-sin_fr,sin_fr,cos_fr);
      for p in stations:
        jj(p) << Meq.Identity(jj);
    # make parmgroups for diagonal and off-diagonal terms
    self.pg_fr  = ParmGroup.ParmGroup(label+"_fr",
            [ jones(src,"RM") for src in sources ],
            table_name="%s_fr.fmep"%label,bookmark=False);

    # make bookmarks
    Bookmarks.make_node_folder("%s Faraday rotation"%label,
      [ jones(src,"fr") for src in sources  ],sorted=True);

    # make solvejobs
    ParmGroup.SolveJob("cal_"+label+"_fr","Calibrate %s"%label,self.pg_fr);

    return jones;
