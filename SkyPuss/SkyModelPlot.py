import sets
from qt import *
from qttable import *
from Qwt4 import *
import math

from Timba.utils import PersistentCurrier

import ModelClasses
import PlotStyles

class SourceMarker (object):
  QwtSymbolStyles = dict(none=QwtSymbol.None,cross=QwtSymbol.XCross,plus=QwtSymbol.Cross,dot=QwtSymbol.Ellipse,
                                        circle=QwtSymbol.Ellipse,square=QwtSymbol.Rect,diamond=QwtSymbol.Diamond);

  def __init__ (self,plotter,src,properties):
    self.src = src;
    self.plotter = plotter;
    self._x0,self._y0 = plotter.canvasCoordinates(src);
    self._size = plotter.getSymbolSize(src);
    self._marker = QwtPlotMarker(self.plotter.plot);
    self._marker.setXValue(self._x0);
    self._marker.setYValue(self._y0);
    self._symbol = QwtSymbol();
    self._font = QApplication.font();
    self._marker_key = None;
    self.setProperties(properties);
    
  def setProperties (self,properties):
    self._properties = properties;
    self.setStyle();
  
  def setStyle (self):
    # first property is always default, so init style from that
    self._show = self._properties[0].show_plot;
    self._style = self._properties[0].style.copy();
    # Then, override attributes in style object with non-default attributes found in each matching property
    # Go in reverse, so 'selected' overrides types overrides tags 
    for prop in self._properties[-1:0:-1]:
      if prop.func(self.src):
        # print self.src.name,"has property",prop.name,prop.apply,prop.show_plot;
        self._show = self._show or prop.show_plot;
        if prop.apply:
          self._style.update(prop.style);
    self._selected = getattr(self.src,'selected',False);
    # setup marker components
    self._setupMarker(self._style,self._show,self._selected);
    # insert marker into plot
    if self._marker_key is None:
      self._marker_key = self.plotter.plot.insertMarker(self._marker);
      
  def _setupMarker (self,style,shown,selected):
    """Sets up the plot marker based on style and "shown" and "selected" property""";
    if shown:
      self._symbol.setStyle(self.QwtSymbolStyles.get(style.symbol,QwtSymbol.Cross));
      self._font.setPointSize(style.label_size);
      symbol_colour = QColor(style.symbol_colour);
      label_colour = QColor(style.label_colour);
##      if selected:
##        self._symbol.setSize(self._size*2);
##        self._symbol.setPen(QPen(symbol_colour,style.symbol_linewidth*2));
##        lab_pen = QPen(label_colour,1);
##        lab_brush = QBrush(QColor(PlotStyles.getContrastColour(style.label_colour)));
      self._symbol.setSize(self._size);
      self._symbol.setPen(QPen(symbol_colour,style.symbol_linewidth));
      lab_pen = QPen(Qt.NoPen);
      lab_brush = QBrush(Qt.NoBrush);
      self._label = PlotStyles.makeSourceLabel(style.label,self.src);
      self._marker.setSymbol(self._symbol);
      self._marker.setLabel(self._label,self._font,label_colour,lab_pen,lab_brush);
      self._marker.setLabelAlignment(Qt.AlignBottom|Qt.AlignRight);
    else:
      self._symbol.setStyle(QwtSymbol.None);
      self._marker.setSymbol(self._symbol);
      self._marker.setLabel("");
    
  def checkSelected (self):
    """Checks the src.selected property, replots marker if changed.
    Returns True is something has changed.""";
    sel = getattr(self.src,'selected',False);
    if self._selected == sel:
      return False;
    self._selected = sel;
    self.setStyle();
    return True;
    
  def changeStyle (self,prop):
    if prop.func(self.src):
      self.setStyle();
      return True;
    return False;

class SkyModelPlotter (QWidget):
  def __init__ (self,parent,*args):
    QWidget.__init__(self,parent,*args);
    lo = QVBoxLayout(self);
    self.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding);
    self.plot  = QwtPlot(self);
    lo.addWidget(self.plot);
    # setup plot properties
    self._zoomer = None;
    self.plot.setCanvasBackground(QColor("#808080"));
    # other internal init
    self._markers = {};
    self._model = None;
    self.projection = None;
    
  def _initProjection (self,projection,lm_max):
    self.projection = projection;
    self.coordinates = "_lm_" + projection;
  
  def setModel (self,model):
    self.plot.clear();
    self._markers = {};
    self._model = model;
    # make sure proper LM coordinates are computed
    proj = self.projection = self.projection or 'ncp';
    lm_max = model.computeSourceLMs(self.projection);
    # compute min/max brightness
    b = [ src.brightness() for src in model.sources ];
    self._min_bright = min(b);
    self._max_bright = max(b);
    # reset canvas as needed
    self._initProjection(proj,lm_max);
    # make items for every object in the model
    for isrc,src in enumerate(model.sources):
      self._markers[src.name] = marker = SourceMarker(self,src,model.properties);
    # reset plot limits
    self.plot.enableAxis(QwtPlot.yLeft,False);
    self.plot.enableAxis(QwtPlot.xBottom,False);
    self.plot.setAxisScale(QwtPlot.yLeft,-lm_max*1.1,lm_max*1.1);
    self.plot.setAxisScale(QwtPlot.xBottom,-lm_max*1.1,lm_max*1.1);
    # update the plot
    self.plot.replot();
    # attach zoomer if None
    if self._zoomer is None:
      self._zoomer = QwtPlotZoomer(self.plot.canvas());
    self._zoomer.setZoomBase();
    
  def updateModelSelection (self):
    """This is callled when some other widget changes the set of selected model sources""";
    # call checkSelected() on all plot markers, replot if any return True
    if filter(lambda marker:marker.checkSelected(),self._markers.itervalues()):
      self.plot.replot();
      
  def changePropertyStyle (self,prop):
    if filter(lambda marker:marker.changeStyle(prop),self._markers.itervalues()):
      self.plot.replot();
      
  def canvasCoordinates (self,source):
    """Returns position of source in the canvas's coordinate system""";
    return getattr(source,self.coordinates);

  def getSymbolSize (self,src):
    return (math.log10(src.brightness()) - math.log10(self._min_bright))*2 + 3;
