import ModelClasses
from qt import *
import math

# string used to indicate default value of an attribute
DefaultValue = "default";
# string used to indicate "none" value of an attribute
NoneValue = "none";

# definitive list of style attributes
StyleAttributes = [ "symbol","symbol_colour","symbol_size","symbol_linewidth","label","label_colour","label_size" ];

# dict of attribute labels (i.e. for menus, column headings and such)
StyleAttributeLabels = dict(symbol="symbol",symbol_colour="colour",symbol_size="size",symbol_linewidth="line width",
                                           label="label",label_colour="colour",label_size="size");
# dict of attribute types. Any attribute not in this dict is of type str.
StyleAttributeTypes = dict(symbol_size=int,symbol_linewidth=int,label_size=int);

# list of known colours
ColourList = [ "black","blue","green","cyan","red","purple","yellow","white" ];
DefaultColour = QColor("black");

# dict and method to pick a contrasting colour (i.e. suitable as background for specified colour)
ContrastColour = dict(white="#404040",yellow="#404040");
DefaultContrastColour = "#B0B0B0";

def getContrastColour (colour):
  return ContrastColour.get(colour,DefaultContrastColour);


# dict of possible user settings for each attribute
StyleAttributeOptions = dict(
  symbol = [ DefaultValue,NoneValue,"cross","plus","dot","circle","square","diamond" ],
  symbol_colour =  [ DefaultValue ] + ColourList,
  label = [ DefaultValue,NoneValue,"%N","%N Ia=%B","%N Ia=%B r=%R'" ],
  label_colour =  [ DefaultValue ] + ColourList,
  label_size = [ DefaultValue,"6","8","10","12","14" ],
);

DefaultPlotAttrs = dict(symbol=None,symbol_colour=DefaultColour,symbol_size=5,symbol_linewidth=0,
                                    label=None,label_colour=DefaultColour,label_size=10);

class PlotStyle (ModelClasses.ModelItem):
  optional_attrs = DefaultPlotAttrs;
  
  def copy (self):
    return PlotStyle(**dict([(attr,getattr(self,attr,default)) for attr,default in DefaultPlotAttrs.iteritems()]))
    
  def update (self,other):
    for attr in DefaultPlotAttrs.iterkeys():
      val = getattr(other,attr,None);
      if val is not None and val != DefaultValue:
        setattr(self,attr,val);
  
PlotStyle.registerClass();

BaselinePlotStyle = PlotStyle(symbol="plus",symbol_colour="yellow",symbol_size=2,symbol_linewidth=0,
                                              label="%N",label_colour="blue",label_size=10);

SelectionPlotStyle = PlotStyle(symbol=DefaultValue,symbol_colour="cyan",symbol_size=DefaultValue,symbol_linewidth=DefaultValue,
                                              label=DefaultValue,label_colour="green",label_size=12);
                                              
DefaultPlotStyle = PlotStyle(symbol=DefaultValue,symbol_colour=DefaultValue,symbol_size=DefaultValue,symbol_linewidth=DefaultValue,
                                              label=DefaultValue,label_colour=DefaultValue,label_size=DefaultValue);

# cache of precompiled labels
_compiled_labels = {};

# label replacements
_label_keys = {   "%N": lambda src:src.name,
                            "%B": lambda src:"%.2g"%src.brightness(),
                            "%R": lambda src:(hasattr(src,'r') and "%.2g"%(src.r/math.pi*180*60)) or "",
                            "%T": lambda src:src.typecode,
                            "%I": lambda src:"%.2g"%getattr(src.flux,'I',0),
                            "%Q": lambda src:"%.2g"%getattr(src.flux,'Q',0),
                            "%U": lambda src:"%.2g"%getattr(src.flux,'U',0),
                            "%V": lambda src:"%.2g"%getattr(src.flux,'V',0),
};
  
def makeSourceLabel (label,src):
  if label == NoneValue:
    return "";
  global _label_keys;
  lbl = label;
  for key,func in _label_keys.iteritems():
    if lbl.find(key) >= 0:
      lbl = lbl.replace(key,func(src));
  return lbl;
