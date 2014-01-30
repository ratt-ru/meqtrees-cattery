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

from Timba.TDL import *
import random
import math

class ErrorGenerator (object):
  def __init__ (self,name,nominal_value,typical_error,unit=None,label=None,**kw):
    label = label or "err-"+name.replace(' ','_');
    self.tdloption_namespace = label;
    self.name = name;
    self.value0 = nominal_value;
    self.opts = [];
    self.factor = 1;
    self.unit = None;
    if unit:
      if isinstance(unit,(list,tuple)):
        self.unit = unit[0];
        self.factor = unit[1];
      else:
        self.unit = unit;
        
  def make_label (self,label):
    label = label%self.name;
    if self.unit:
      label = "%s, %s"%(label,self.unit);
    return label;
    
  def make_node (self,node):
    raise TypeError,"make_node() must be implemented in subclass";


class NoError (ErrorGenerator):
  def make_node (self,node,**kw):
    node << self.value0;

class FixedOffset (ErrorGenerator):
  def __init__ (self,name,nominal_value=0,typical_error=0,**kw):
    ErrorGenerator.__init__(self,name,nominal_value,typical_error=[0],**kw);
    self.opts.append(TDLOption('offset',self.make_label("%s offset"),
                     typical_error,more=float,namespace=self));
  
  def make_node (self,node,**kw):
    node << (self.factor*self.offset + self.value0);

class ListOfValues (ErrorGenerator):
  def __init__ (self,name,nominal_value=0,typical_error=0,**kw):
    ErrorGenerator.__init__(self,name,nominal_value,typical_error=[0],**kw);
    self.opts.append(TDLOption('values_str',self.make_label("%s values"),
                     map(str,typical_error) if isinstance(typical_error,(list,tuple)) else [str(typical_error)],
                     more=str,
                     doc="Enter list of values, separated by spaces. Last value will be reused.",
                     namespace=self));
    self.ngen = 0;
    self.values = None;
    
  
  def make_node (self,node,**kw):
    if self.values is None:
      self.values = map(float,self.values_str.split(" "));
    if self.ngen >= len(self.values):
      node << self.values[-1]*self.factor;
    else:
      node << self.values[self.ngen]*self.factor;
    self.ngen+=1;

class RandomError (ErrorGenerator):
  def __init__ (self,name,nominal_value=0,typical_error=0,**kw):
    ErrorGenerator.__init__(self,name,nominal_value,typical_error,**kw);
    if isinstance(typical_error,(list,tuple)):
      self.opts.append(TDLOption('minval',self.make_label("Min %s"),
                            [typical_error[0]],more=float,namespace=self));
      self.opts.append(TDLOption('maxval',self.make_label("Max %s"),
                            [typical_error[1]],more=float,namespace=self));
      self._symmetric = False;
    else:
      self.opts.append(TDLOption('maxerr',self.make_label("Max %s error"),
                            [typical_error],more=float,namespace=self));
      self._symmetric = True;
    
  def get_range (self):  
    if self._symmetric:
      return self.factor*(self.value0-self.maxerr), \
             self.factor*(self.value0+self.maxerr);
    else:
      return self.factor*self.minval, \
             self.factor*self.maxval;
  
  def make_node (self,node,**kw):
    node << random.uniform(*self.get_range());

class SineError (RandomError):
  def __init__ (self,name,nominal_value=0,typical_error=0,**kw):
    RandomError.__init__(self,name,nominal_value,typical_error,**kw);
    self.opts.append(TDLOption("min_period","Min period of %s variation, hours"%name,
                          [1,2],more=float,namespace=self));
    self.opts.append(TDLOption("max_period","Max period of %s variation, hours"%name,
                          [2,4],more=float,namespace=self));
    
  def make_node (self,node,**kw):
    ns = node.Subscope();
    minval,maxval = self.get_range();
    ns.ampl << (maxval-minval)/2.;
    ns.offset << (maxval+minval)/2.;
    period = random.uniform(self.min_period,self.max_period)*3600;  # period in seconds
    # pick a random starting phase
    p0 = random.uniform(0,2*math.pi);
    node << ns.offset + ns.ampl*Meq.Sin(Meq.Time()*(2*math.pi/period)+p0);
    return node;

# This list shows the available generator classes
generator_classes = [
  (NoError,       'no error'),
  (FixedOffset,   'static offset'),
  (ListOfValues,  'list of values'),
  (RandomError,   'random error, constant in time'),
  (SineError,     'periodically varying error')
];



class Selector (object):
  def __init__ (self,name,nominal_value=0,typical_error=0,label=None,**kw):
    label = label or "err-"+name.replace(' ','_');
    self.tdloption_namespace = label;
    # create selector from list of generator classes
    model_option = TDLOption("error_model","Error model for %s"%name,\
                            dict(generator_classes),default=NoError,namespace=self);
    self.opts = [ model_option ];
    # create instances of each class
    self._generators = dict();
    for genclass,desc in generator_classes:
      gen = genclass(name,nominal_value,typical_error,label=label,**kw);
      self._generators[genclass] = gen;
      self.opts += gen.opts;
    # show/hide options interactively
    model_option.when_changed(self._show_gen_options);
  
  def has_errors (self):
    """Returns true if a non-trivial error model is selected""";
    return self.error_model is not NoError;
  
  def make_node (self,node,**kw):
    self._generators[self.error_model].make_node(node);
    
  def node_maker (self):
    return self._generators[self.error_model].make_node;
    
  def options (self):
    return self.opts;
  
  def _show_gen_options (self,model):
    genobj0 = None;
    # hide all "other" options
    for genclass,genobj in self._generators.iteritems():
      if genclass == model:
        genobj0 = genobj;
      else:
        for opt in genobj.opts:
          opt.hide();
    # show model options
    if genobj0:
      for opt in genobj0.opts:
        opt.show();
    
    