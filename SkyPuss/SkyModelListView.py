import sets
import math
from qt import *
from qttable import *

import Timba.utils
from Timba.utils import PersistentCurrier

import ModelClasses
import PlotStyles

_verbosity = Timba.utils.verbosity(name="lsmlv");
dprint = _verbosity.dprint;
dprintf = _verbosity.dprintf;


ViewColumns = [ "name","RA","Dec","r","type","Iapp","I","Q","U","V","RM","spi","shape","tags" ];

for icol,col in enumerate(ViewColumns):
    globals()["Column%s"%col.capitalize()] = icol;
ColumnExtra = len(ViewColumns);

class SkyModelListView (QListView):
  """This implements a QListView for sky models""";
  def __init__ (self,*args):
    QListView.__init__(self,*args);
    self._currier = PersistentCurrier();
    # insert columns and init constants (self.ColumnName etc.) for column numbers
    for icol,col in enumerate(ViewColumns):
        self.addColumn((col == "Iapp" and "I(app)") or col);
        self.setColumnWidthMode(icol,QListView.Maximum);
    self.setColumnAlignment(ColumnR,Qt.AlignRight);
    self.setColumnAlignment(ColumnType,Qt.AlignHCenter);
    self.addColumn(""); # extra column
    self.setColumnWidthMode(ColumnExtra,QListView.Maximum);
    # _column_enabled[i] is True if column is available in the model. 
    # _column_show[i] is True if column is currently being shown (via a view control)
    self._column_enabled = [True]*(ColumnExtra+1);
    self._column_shown    = [True]*(ColumnExtra+1);
    # other listview init
    self.header().show();
    self.setAllColumnsShowFocus(True);
    self.setShowToolTips(True);
    # self.elv.setResizeMode(QListView.LastColumn);
    self.setSelectionMode(QListView.Extended);
    self._updating_selection = False;
    self.setRootIsDecorated(False);
    # connect signals to track selected sources
    QObject.connect(self,SIGNAL("selectionChanged()"),self._selectionChanged);
    # add "View" controls for different column categories
    self._column_views = [];
    self._column_widths = {};
    self.addColumnCategory("Position",[ColumnRa,ColumnDec]);
    self.addColumnCategory("Type",[ColumnType]);
    self.addColumnCategory("Flux",[ColumnIapp,ColumnI]);
    self.addColumnCategory("Polarization",[ColumnQ,ColumnU,ColumnV,ColumnRm]);
    self.addColumnCategory("Spectrum",[ColumnSpi]);
    self.addColumnCategory("Shape",[ColumnShape]);
    self.addColumnCategory("Tags",[ColumnTags]);
 
  def _showColumn (self,col,show=True):
    """Shows or hides the specified column. 
    (When hiding, saves width of column to internal array so that it can be restored properly.)"""
    if show:
      if not self.columnWidth(col):
        self.setColumnWidth(col,self._column_widths[col]);
        self.setColumnWidthMode(col,QListView.Maximum);
    else:
      if self.columnWidth(col):
        self._column_widths[col] = self.columnWidth(col);
        self.hideColumn(col);

  def _enableColumn (self,column,enable=True):
    self._column_enabled[column] = enable;
    self._showColumn(column,enable and self._column_shown[column]);
 
  def _showColumnCategory (self,columns,show):
    for col in columns:
      self._column_shown[col] = show;
      self._showColumn(col,self._column_enabled[col] and show);
          
  def _selectionChanged (self):
    if self._updating_selection:
      return;
    item = self.firstChild();
    while item:
      item._src.select(item.isSelected());
      item = item.nextSibling();
    self.emit(PYSIGNAL("modelSelectionChanged()"),(self,));
  
  def addColumnCategory (self,name,columns):
    qa = QAction(name,0,self);
    qa.setToggleAction(True);
    qa.setOn(True);
    QObject.connect(qa,SIGNAL("toggled(bool)"),self._currier.curry(self._showColumnCategory,columns));
    self._column_views.append((name,qa,columns));

  def setModel (self,model):
    self.clear();
    for src in model.sources:
      self.addSource(src);
    # show/hide columns based on tag availability
    self._enableColumn(ColumnIapp,'Iapp' in model.tagnames);
    self._enableColumn(ColumnR,'r' in model.tagnames);
    self.setSorting(('Iapp' in model.tagnames and ColumnIapp) or ColumnI,False);
      
  def addSource (self,src):
    return SkyModelListViewItem(src,self,self.lastItem());

  def addColumnViewActionsTo (self,menu):
    for name,qa,columns in self._column_views:
      qa.addTo(menu);
      
  def updateModelSelection (self):
    """This is called when some other widget changes the set of selected model sources""";
    self._updating_selection = True;
    self.clearSelection();
    item = self.firstChild();
    while item:
      self.setSelected(item,item._src.selected);
      item = item.nextSibling();
    self._updating_selection = False;
    
  def refreshSelectedItems (self):
    item = self.firstChild();
    while item:
      if item.isSelected():
        item.setSource(item._src);
      item = item.nextSibling();
      
  TagsWithOwnColumn = sets.Set(["Iapp","r"]);
    
class SkyModelListViewItem (QListViewItem):
  def __init__ (self,src,*args):
    QListViewItem.__init__(self,*args);
    self._src = src;
    # array of actual (i.e. numeric) column values
    self._values = [None]*(ColumnExtra+1);
    self.setSource(src);
    
  def setSource (self,src):
    # name
    self.setColumn(ColumnName,src.name);
    # coordinates
    self.setColumn(ColumnRa,src.pos.ra,"%2dh%02dm%05.2fs"%src.pos.ra_hms());
    self.setColumn(ColumnDec,src.pos.dec,("%+2d"+unichr(0xB0)+"%02d'%05.2f\"")%src.pos.dec_dms());
    if hasattr(src,'r'):
      self.setColumn(ColumnR,src.r,"%.1f'"%(src.r*180*60/math.pi));
    # type
    self.setColumn(ColumnType,src.typecode);
    # flux
    if hasattr(src,'Iapp'):
      self.setColumn(ColumnIapp,src.Iapp,"%.2g"%src.Iapp);
    self.setColumn(ColumnI,src.flux.I,"%.2g"%src.flux.I);
    # polarization
    if isinstance(src.flux,ModelClasses.Polarization):
      self.setColumn(ColumnQ,src.flux.Q,"%.2g"%src.flux.Q);
      self.setColumn(ColumnU,src.flux.U,"%.2g"%src.flux.U);
      self.setColumn(ColumnV,src.flux.V,"%.2g"%src.flux.V);
      if hasattr(src.flux,'rm'):
        self.setColumn(ColumnRm,src.flux.rm,"%.2f"%src.flux.rm);
    # spi
    if isinstance(src.spectrum,ModelClasses.SpectralIndex):
      self.setColumn(ColumnSpi,src.spectrum.spi,"%.2f"%getattr(src.spectrum,'spi',''));
    # shape
    shape = getattr(src,'shape',None);
    if isinstance(shape,ModelClasses.ModelItem):
      shapeval = [ val for attr,val in shape.getAttributes() ];
      shapestr = shape.strAttributes(label=False);
      self.setColumn(ColumnShape,shapeval,shapestr);
    # Tags. Tags are all extra attributes that do not have a dedicated column (i.e. not Iapp or r), and do not start
    # with "_" (which is reserved for internal attributes)
    truetags = [];
    falsetags = [];
    othertags = [];
    for attr,val in src.getExtraAttributes():
      if attr[0] != "_" and attr not in self.listView().TagsWithOwnColumn:
        if val is False:
          falsetags.append("-"+attr);
        elif val is True:
          truetags.append("+"+attr);
        else:
          othertags.append("%s=%s"%(attr,str(val)));
    for tags in truetags,falsetags,othertags:
      tags.sort();
    self.setColumn(ColumnTags,tags," ".join(truetags+falsetags+othertags));
    
  def setColumn (self,icol,value,text=None):
    """helper function to set the value of a column""";
    if text is None:
      text = str(value);
    self.setText(icol,text);
    self._values[icol] = value;

  def compare (self,other,icol,ascending):
    """Reimplement the compare method to compare true numeric values""";
    if isinstance(other,SkyModelListViewItem):
      return cmp(self._values[icol],other._values[icol]);
    else:
      return cmp(self.text(icol),other.text(icol));
  
class ModelPropsTable (QWidget):
  EditableAttrs = [ attr for attr in PlotStyles.StyleAttributes if attr in PlotStyles.StyleAttributeOptions ];
  ColList = 2;
  ColPlot = 3;
  ColApply = 4;
  AttrByCol = dict([(i+5,attr) for i,attr in enumerate(EditableAttrs)]);

  def __init__ (self,parent,*args):
    QWidget.__init__(self,parent,*args);
    self.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding);
    lo = QVBoxLayout(self);
    lo1 = QHBoxLayout(lo);
    # lo1.addLayout(lo1);
    lbl = QLabel(QString("<nobr><b>Source properties table:</b></nobr>"),self);
    lo1.addWidget(lbl,0);
    lo1.addStretch(1);
    # add show/hide button
    self._showattrbtn = QPushButton(self,"");
    self._showattrbtn.setMinimumWidth(256);
    lo1.addWidget(self._showattrbtn,0);
    lo1.addStretch();
    QObject.connect(self._showattrbtn,SIGNAL("clicked()"),self._hideShowAttributes);
    # add table
    self.table  = QTable(self);
    lo.addWidget(self.table);
    QObject.connect(self.table,SIGNAL("valueChanged(int,int)"),self._valueChanged);
    QObject.connect(self.table.verticalHeader(),SIGNAL("clicked(int)"),self._selectRow);
    self.table.setSelectionMode(QTable.NoSelection);
    # setup basic columns
    self.table.setNumCols(5+len(self.EditableAttrs));
    hdr = self.table.horizontalHeader();
    for i,label in enumerate(("property","total","","","")):
      hdr.setLabel(i,label);
    self.table.hideColumn(self.ColApply);
    # setup columns for editable property attributes
    for i,attr in self.AttrByCol.iteritems():
      hdr.setLabel(i,PlotStyles.StyleAttributeLabels[attr]);
      self.table.hideColumn(i);
    # other internal init
    self._attrs_shown = True;
    self._hideShowAttributes();
    self.model = None;
      
  def clear (self):
    self.table.setNumRows(0);
    self.model = None;
    
  def setModel (self,model):
    self.clear();
    self.model = model;
    self.table.setNumRows(len(model.properties));
    for i,prop in enumerate(model.properties):
      item = [None]*5;
      item[0] = QTableItem(self.table,QTableItem.Never,prop.name);
      item[1] = QTableItem(self.table,QTableItem.Never,str(prop.total));
      item[self.ColList] = QCheckTableItem(self.table,"list"); 
      item[self.ColList].setChecked(bool(prop.show_list));
      item[self.ColList].setEnabled(isinstance(prop.show_list,bool));
      item[self.ColPlot] = QCheckTableItem(self.table,"plot");
      item[self.ColPlot].setChecked(bool(prop.show_plot));
      item[self.ColApply] = QCheckTableItem(self.table,"apply");
      item[self.ColApply].setChecked(prop.apply);
      item[self.ColApply].setEnabled(i!=0);  # default always applies
      for j,it in enumerate(item):
        it.setReplaceable(False);
        self.table.setItem(i,j,it);
      for j,attr in self.AttrByCol.iteritems():
        # get list of options for this style attribute. If dealing with first property (i==0), which is
        # the "all sources" property, then remove the "default" option (which is always first in the list)
        options = PlotStyles.StyleAttributeOptions[attr];
        if i == 0:
          options = options[1:];
        # make item
        item = QComboTableItem(self.table,QStringList.split("|","|".join(options)));
        item.setReplaceable(False);
        # the "label" option is also editable
        if attr == "label":
          item.setEditable(True);
        item.setCurrentItem(str(getattr(prop.style,attr)));
        self.table.setItem(i,j,item);
    for col in range(self.table.numCols()):
      self.table.adjustColumn(col);
      
  def _hideShowAttributes (self):
    if self._attrs_shown:
      self._attrs_shown = False;
      self.table.hideColumn(self.ColApply);
      for col in self.AttrByCol.iterkeys():
        self.table.hideColumn(col);
      self._showattrbtn.setText("Show plot styles >>");
    else:
      self._attrs_shown = True;
      self.table.showColumn(self.ColApply);
      for col in self.AttrByCol.iterkeys():
        self.table.showColumn(col);
      self._showattrbtn.setText("<< Hide plot styles");
      
  def _valueChanged (self,row,col):
    """Called when a cell has been edited""";
    prop = self.model.properties[row];
    item = self.table.item(row,col);
    if col == self.ColList:
      prop.show_list = item.isChecked();
      self.emit(PYSIGNAL("listProperty()"),(prop,));
      return;
    elif col == self.ColPlot:
      prop.show_plot = item.isChecked();
      self.emit(PYSIGNAL("changePropertyStyle()"),(prop,));
    elif col == self.ColApply:
      prop.apply = item.isChecked();
      # enable/disable editable cells
      for j in self.AttrByCol.keys():
        self.table.item(row,j).setEnabled(prop.apply);
        self.table.updateCell(row,j);
      self.emit(PYSIGNAL("changePropertyStyle()"),(prop,));
    elif col in self.AttrByCol:
      attr = self.AttrByCol[col];
      setattr(prop.style,attr,PlotStyles.StyleAttributeTypes.get(attr,str)(item.currentText()));
      self.emit(PYSIGNAL("changePropertyStyle()"),(prop,));
      self.emit(PYSIGNAL("modelUpdated()"),());
      
  def _selectRow (self,row):
    # see which property the selected row corresponds to, ignore if "selected" property
    prop = self.model.properties[row];
    if prop is self.model.selprop:
      return;
    # select all sources with that property
    for src in self.model.sources:
      src.select(prop.func(src));
    self.emit(PYSIGNAL("modelSelectionChanged()"),(self,));
    self.model.selprop.computeTotal(self.model.sources);
    self.updateModelSelection();
      
  def updateModelSelection (self):
    """This is called when some other widget changes the set of selected model sources""";
    self.table.clearSelection();
    for i,prop in enumerate(self.model.properties):
      if prop is self.model.selprop:
        self.table.item(i,1).setText(str(prop.total));
        self.table.updateCell(i,1);
    
