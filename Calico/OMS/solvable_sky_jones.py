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
  """This implements a full 2x2 Jones matrix with solvable real and imaginary parts""";
  def __init__ (self,label):
    self.tdloption_namespace = label+".fullrealimag";
    self.subset_selector = Meow.MeqMaker.SourceSubsetSelector("Apply this Jones term to a subset of sources",
                            tdloption_namespace=self.tdloption_namespace);
    self.options = [
      TDLOption("matrix_type","Matrix type",["complex","real"],namespace=self),
      TDLOption("init_diag","Initial value, diagonal",[0,1],default=1,more=complex,namespace=self),
      TDLOption("init_offdiag","Initial value, off-diagonal",[0,1],default=0,more=complex,namespace=self),
    ] + self.subset_selector.options;
    self._offdiag = True;

  def compile_options (self):
    return self.options;

  def compute_jones (self,jones,sources,stations=None,tags=None,label='',**kw):
    stations = stations or Context.array.stations();
    # figure out which sources to apply to
    sources = self.subset_selector.filter(sources);
    is_complex = self.matrix_type != "real";
    # set up qualifier labels depending on polarization definition
    if Context.observation.circular():
      x,y,X,Y = 'R','L','R','L';
    else:
      x,y,X,Y = 'x','y','X','Y';
    xx,xy,yx,yy = x+x,x+y,y+x,y+y;
    rxx,rxy,ryx,ryy = [ q+":r" for q in xx,xy,yx,yy ];
    ixx,ixy,iyx,iyy = [ q+":i" for q in xx,xy,yx,yy ];
    # prepare parm definitions for real and imag parts of diagonal and off-diagonal elements
    tags = NodeTags(tags) + "solvable";
    diag_pdefs =     ( Meq.Parm(complex(self.init_diag).real,tags=tags+"diag real"),
                       Meq.Parm(complex(self.init_diag).imag,tags=tags+"diag imag") );
    offdiag_pdfefs = ( Meq.Parm(complex(self.init_offdiag).imag,tags=tags+"offdiag real"),
                       Meq.Parm(complex(self.init_offdiag).imag,tags=tags+"offdiag imag") );
    # loop over sources
    parms_diag = [];
    parms_offdiag = [];
    subgroups_diag = [];
    subgroups_offdiag = [];

    # Define a function to put together a matrix element, depending on whether we're in real or complex mode.
    # This is also responsible for appending parms to the appropriate parm and subgroup lists.
    # Note that for the purely-real case we still create parms called 'J:xx:r' (and not 'J:xx' directly),
    # this is to keep MEP tables mututally-compatible naming wise.
    if self.matrix_type == "real":
      def make_element (jj,parmdef,*parmlists):
        parm = jj('r') << parmdef[0];
        for plist in parmlists:
          plist.append(parm);
        return jj << Meq.Identity(parm);
    else:
      def make_element (jj,parmdef,*parmlists):
        parms = [ jj('r') << parmdef[0],jj('i') << parmdef[1] ];
        for plist in parmlists:
          plist += parms;
        return jj << Meq.ToComplex(*parms);
        
    for src in sources:
      # now loop to create nodes
      sgdiag = [];
      sgoff = [];
      for p in stations:
        jj = jones(src,p);
        jj << Meq.Matrix22(
          make_element(jj(xx),diag_pdefs,parms_diag,sgdiag),
          make_element(jj(xy),offdiag_pdefs,parms_offdiag,sgoff) if self._offdiag else 0,
          make_element(jj(yx),offdiag_pdefs,parms_offdiag,sgoff) if self._offdiag else 0,
          make_element(jj(yy),diag_pdefs,parms_diag,sgdiag)
        );
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
    if self._offdiag:
      self.pg_offdiag  = ParmGroup.ParmGroup(label+"_offdiag",parms_offdiag,
                        subgroups=subgroups_offdiag,
                        table_name="%s_offdiag.fmep"%label,bookmark=False);

    # make bookmarks
    Bookmarks.make_node_folder("%s diagonal terms"%label,
      [ jones(src,p,zz) for src in sources  
        for p in stations for zz in "xx","yy" ],sorted=True);
    if self._offdiag:
      Bookmarks.make_node_folder("%s off-diagonal terms"%label,
        [ jones(src,p,zz) for src in sources  
          for p in stations for zz in "xy","yx" ],sorted=True);

    # make solvejobs
    ParmGroup.SolveJob("cal_"+label+"_diag","Calibrate %s diagonal terms"%label,self.pg_diag);
    if self._offdiag:
      ParmGroup.SolveJob("cal_"+label+"_offdiag","Calibrate %s off-diagonal terms"%label,self.pg_offdiag);

    return jones;


class DiagRealImag (FullRealImag):
  """This implements a diagonal 2x2 Jones matrix with solvable real and imaginary parts""";
  # note that this is not really proper OO design, but it's convenient
  def __init__ (self,label):
    self.tdloption_namespace = label+".diagrealimag";
    self.subset_selector = Meow.MeqMaker.SourceSubsetSelector("Apply this Jones term to a subset of sources",
                            tdloption_namespace=self.tdloption_namespace);
    self.options = [
      TDLOption("matrix_type","Matrix type",["complex","real"],namespace=self),
      TDLOption("init_diag","Initial value",[0,1],default=1,more=float,namespace=self),
    ] + self.subset_selector.options;
    self._offdiag = False;
    self.init_offdiag = 0;


      
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
