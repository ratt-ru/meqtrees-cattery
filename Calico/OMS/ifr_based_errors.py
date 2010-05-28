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
"""This implements two Visibility Processing modules.
IfrGains implements solvable ifr-based gains
IfrBiases implements solvable ifr-based biases
""";

from Timba.TDL import *
import Meow
from Meow import Context
from Meow import ParmGroup,Bookmarks
from Meow.Parameterization import resolve_parameter

class IfrGains (object):
  def __init__ (self):
    self.options = [];

  def runtime_options (self):
    return self.options;

  def init_nodes (self,ns,tags=None,label=''):
    ifrs = Context.array.ifrs();
    G = ns.gain;
    if not G(*ifrs[0]).initialized():
      def1 = Meow.Parm(1,tags=tags);
      def0 = Meow.Parm(0,tags=tags);
      gains = [];
      for p,q in ifrs:
        gg = [];
        for corr in Context.correlations:
          if corr in Context.active_correlations:
            cc = corr.lower();
            gain_ri  = [ resolve_parameter(cc,G(p,q,cc,'r'),def1,tags="ifr gain real"),
                        resolve_parameter(cc,G(p,q,cc,'i'),def0,tags="ifr gain imag") ];
            gg.append(G(p,q,cc) << Meq.ToComplex(*gain_ri));
            gains += gain_ri;
          else:
            gg.append(1);
        G(p,q) << Meq.Matrix22(*gg);

      pg_ifr_ampl = ParmGroup.ParmGroup(label,gains,
                            individual_toggles=False,
                            table_name="%s.fmep"%label,bookmark=False);
      Bookmarks.make_node_folder("%s (interferometer-based gains)"%label,
        [ G(p,q) for p,q in ifrs ],sorted=True,nrow=2,ncol=2);
      ParmGroup.SolveJob("cal_%s"%label,
                        "Calibrate %s (interferometer-based gains)"%label,pg_ifr_ampl);
    return G;
    

  def process_visibilities (self,nodes,input_nodes,ns=None,
                             ifrs=None,tags=None,label='',**kw):
    ns = ns or nodes.Subscope();
    ifrs = ifrs or Context.array.ifrs();
    label = label or 'ifr_gains';
    G = self.init_nodes(ns,tags,label);

    # apply gains
    for p,q in ifrs:
      nodes(p,q) << G(p,q)*input_nodes(p,q);

    return nodes;

  def correct_visibilities (self,nodes,input_nodes,ns=None,
                            ifrs=None,tags=None,label='',**kw):
    ns = ns or nodes.Subscope();
    ifrs = ifrs or Context.array.ifrs();
    G = self.init_nodes(ns,tags,label);

    for p,q in ifrs:
      nodes(p,q) << input_nodes(p,q)/G(p,q);
    
    return nodes;

class IfrBiases (object):
  def __init__ (self):
    self.options = [];

  def runtime_options (self):
    return self.options;

  def init_nodes (self,ns,tags=None,label=''):
    ifrs = Context.array.ifrs();
    C = ns.bias;
    if not C(*ifrs[0]).initialized():
      def0 = Meow.Parm(0,tags=tags);
      biases = [];
      for p,q in ifrs:
        gg = [];
        for corr in Context.correlations:
          if corr in Context.active_correlations:
            cc = corr.lower();
            bias_ri  = [ resolve_parameter(cc,C(p,q,cc,'r'),def0,tags="ifr bias real"),
                        resolve_parameter(cc,C(p,q,cc,'i'),def0,tags="ifr bias imag") ];
            gg.append(C(p,q,cc) << Meq.ToComplex(*bias_ri));
            biases += bias_ri;
          else:
            gg.append(0);
        C(p,q) << Meq.Matrix22(*gg);

      pg_ifr_ampl = ParmGroup.ParmGroup(label,biases,
                            individual_toggles=False,
                            table_name="%s.fmep"%label,bookmark=False);
      Bookmarks.make_node_folder("%s (interferometer-based biases)"%label,
        [ C(p,q) for p,q in ifrs ],sorted=True,nrow=2,ncol=2);
      ParmGroup.SolveJob("cal_%s"%label,
                        "Calibrate %s (interferometer-based biases)"%label,pg_ifr_ampl);
    return C;

  def process_visibilities (self,nodes,input_nodes,ns=None,
                             ifrs=None,tags=None,label='',**kw):
    ns = ns or nodes.Subscope();
    ifrs = ifrs or Context.array.ifrs();
    label = label or 'ifr_biases';
    C = self.init_nodes(ns,tags,label);

    for p,q in ifrs:
      nodes(p,q) << input_nodes(p,q) + C(p,q);
    
    return nodes;

  def correct_visibilities (self,nodes,input_nodes,ns=None,
                            ifrs=None,tags=None,label='',**kw):
    ns = ns or nodes.Subscope();
    ifrs = ifrs or Context.array.ifrs();

    C = self.init_nodes(ns,tags,label);

    for p,q in ifrs:
      nodes(p,q) << input_nodes(p,q) - C(p,q);
    
    return nodes;

