import os.path
import time
import sys

from qt import *
from qttable import *

import SkyPuss
import ModelHTML
import ModelClasses
import Widgets
from SkyModelListView import *
from SkyModelPlot import *
from SkyPuss import Config,pixmaps,dprint,dprintf

class BusyIndicator (object):
  def __init__ (self):
    QApplication.setOverrideCursor(QCursor(Qt.WaitCursor));
  def __del__ (self):
    QApplication.restoreOverrideCursor();

class MainWindow (QMainWindow):
  ViewModelColumns = [ "name","RA","Dec","type","Iapp","I","Q","U","V","RM","spi","shape" ];
  def __init__ (self,parent,hide_on_close=False):
    QMainWindow.__init__(self,parent);
    # init column constants
    for icol,col in enumerate(self.ViewModelColumns):
        setattr(self,"Column%s"%col.capitalize(),icol);
    # init GUI
    self.setCaption("SkyPuss");
    # self.setIcon(pixmaps.purr_logo.pm());
    cw = QWidget(self);
    self.setCentralWidget(cw);
    cwlo = QVBoxLayout(cw);
    cwlo.setMargin(5);
    # make splitter
    spl1 = QSplitter(Qt.Vertical,cw);
    cwlo.addWidget(spl1);
    # Create listview of LSM entries
    self.lv = SkyModelListView(spl1);
    #self.connect(self.lv,SIGNAL("selectionChanged()"),self._entrySelectionChanged);
    #self.connect(self.lv,SIGNAL("doubleClicked(QListViewItem*,const QPoint&,int)"),self._viewEntryItem);
    #self.connect(self.lv,SIGNAL("returnPressed(QListViewItem*)"),self._viewEntryItem);
    #self.connect(self.lv,SIGNAL("spacePressed(QListViewItem*)"),self._viewEntryItem);
    #self.connect(self.lv,SIGNAL('contextMenuRequested(QListViewItem*,const QPoint &,int)'),
    #                 self._showItemContextMenu);
    
    # split bottom pane
    spl2 = QSplitter(Qt.Horizontal,spl1);
    # add plot
    self.skyplot = SkyModelPlotter(spl2);
    # add properties table
    self.proptab = ModelPropsTable(spl2);
    # connect signals
    QObject.connect(self.lv,PYSIGNAL("modelSelectionChanged()"),self._updateModelSelection);
    QObject.connect(self.skyplot,PYSIGNAL("modelSelectionChanged()"),self._updateModelSelection);
    QObject.connect(self.proptab,PYSIGNAL("modelSelectionChanged()"),self._updateModelSelection);
    QObject.connect(self.lv,PYSIGNAL("modelUpdated()"),self._modelUpdated);
    QObject.connect(self.skyplot,PYSIGNAL("modelUpdated()"),self._modelUpdated);
    QObject.connect(self.proptab,PYSIGNAL("modelUpdated()"),self._modelUpdated);
    QObject.connect(self.proptab,PYSIGNAL("changePropertyStyle()"),self.skyplot.changePropertyStyle);
    # enable status line
    self.statusBar().show();
    # Create and populate main menu
    menubar = self.menuBar();
    # File menu
    file_menu = QPopupMenu(self);
    menubar.insertItem("&File",file_menu);
    qa_open = QAction("&Open model...",Qt.CTRL+Qt.Key_O,self);   qa_open.addTo(file_menu);
    QObject.connect(qa_open,SIGNAL("activated()"),self.openFile);
    qa_merge = QAction("&Merge in model...",Qt.CTRL+Qt.SHIFT+Qt.Key_O,self);   qa_merge.addTo(file_menu);
    QObject.connect(qa_merge,SIGNAL("activated()"),self.mergeFile);
    QObject.connect(self,PYSIGNAL("hasContent()"),qa_merge.setEnabled);
    file_menu.insertSeparator();
    qa_save = QAction("&Save model",Qt.CTRL+Qt.Key_S,self);   qa_save.addTo(file_menu);    
    QObject.connect(qa_save,SIGNAL("activated()"),self._saveFile);
    QObject.connect(self,PYSIGNAL("isUpdated()"),qa_save.setEnabled);
    qa_save_as = QAction("Save model &as...",0,self);              qa_save_as.addTo(file_menu);    
    QObject.connect(qa_save_as,SIGNAL("activated()"),self.saveFileAs);
    QObject.connect(self,PYSIGNAL("hasContent()"),qa_save_as.setEnabled);
    qa_save_selection_as = QAction("Save selection as...",0,self); qa_save_selection_as.addTo(file_menu);
    QObject.connect(qa_save_selection_as,SIGNAL("activated()"),self._saveSelectionAs);
    QObject.connect(self,PYSIGNAL("hasSelection()"),qa_save_selection_as.setEnabled);
    file_menu.insertSeparator();
    qa_close = QAction("&Close model",Qt.CTRL+Qt.Key_W,self);   qa_close.addTo(file_menu);
    QObject.connect(qa_close,SIGNAL("activated()"),self._closeFile);
    QObject.connect(self,PYSIGNAL("hasContent()"),qa_close.setEnabled);
    qa_quit = QAction("Quit",Qt.CTRL+Qt.Key_Q,self);           qa_quit.addTo(file_menu);    
    QObject.connect(qa_quit,SIGNAL("activated()"),self.close);
    # Edit Menu
    self.edit_menu = QPopupMenu(self);
    menubar.insertItem("&Edit",self.edit_menu);
    qa_add_tag = QAction("&Tag selected..",Qt.CTRL+Qt.Key_T,self); qa_add_tag.addTo(self.edit_menu);
    QObject.connect(qa_add_tag,SIGNAL("activated()"),self._addTagToSelection);
    QObject.connect(self,PYSIGNAL("hasSelection()"),qa_add_tag.setEnabled);
    qa_del_tag = QAction("&Untag selected..",Qt.CTRL+Qt.Key_U,self); qa_del_tag.addTo(self.edit_menu);
    QObject.connect(qa_del_tag,SIGNAL("activated()"),self._removeTagsFromSelection);
    QObject.connect(self,PYSIGNAL("hasSelection()"),qa_del_tag.setEnabled);
    
    # View menu
    self.view_menu = QPopupMenu(self);
    menubar.insertItem("&View",self.view_menu);
    # Export menu
    export_menu = QPopupMenu(self);
    menubar.insertItem("&Export",export_menu);
    self.qa_export_karma = QAction("Export Karma annotations...",0,self); self.qa_export_karma.addTo(export_menu);
    # set initial state
    self.model = None;
    self.filename = None;
    self._open_file_dialog = self._save_as_dialog = None;
    self.emit(PYSIGNAL("isUpdated()"),(False,));
    self.emit(PYSIGNAL("hasContent()"),(False,));
    self.emit(PYSIGNAL("hasSelection()"),(False,));
    # resize selves
    width = Config.getint('main-window-width',512);
    height = Config.getint('main-window-height',512);
    self.resize(QSize(width,height));
    
  def resizeEvent (self,ev):
    QMainWindow.resizeEvent(self,ev);
    sz = ev.size();
    Config.set('main-window-width',sz.width());
    Config.set('main-window-height',sz.height());
    
  def _updateModelSelection (self,origin):
    """Called when a modelSelectionChanged() signal arrives fom one of its child widgets.
    Broadcasts this to other children.""";
    has_sel = bool(self.model and self.model.updateSelection());
    for widget in (self.lv,self.proptab,self.skyplot):
      if widget != origin:
        widget.updateModelSelection();
    self.emit(PYSIGNAL("hasSelection()"),(has_sel,));
    
  def mergeFile (self,filename=None,format=None):
    return self.openFile(filename,format=format,merge=True);
    
  def openFile (self,filename=None,format=None,merge=False):
    if filename is None:
        if not self._open_file_dialog:
            dialog = self._open_file_dialog = QFileDialog(".","*.html",self,None,True);
            dialog.setCaption("Open Sky Model");
            dialog.setMode(QFileDialog.ExistingFile);
            QObject.connect(dialog,SIGNAL("fileSelected(const QString&)"),self.openFile);
        self._open_file_dialog.rereadDir();
        self._open_file_dialog.exec_loop();
        return;
    # try to load the specified file
    filename = str(filename);
    busy = BusyIndicator();
    try:
        model = ModelHTML.loadModel(filename);
    except:
        busy = None;
        QMessageBox.warning(self,"Error loading sky model",
            """Error loading model file %s: %s"""%(filename,str(sys.exc_info()[1])),QMessageBox.Ok,0,0);
        return;
    # add to content
    if merge and self.model:
        self.model.append(model.sources);
        self.statusBar().message("""Merged in %d sources from model file %s"""%(len(model.sources),filename),3000);
    else:
        self.statusBar().message("""Loaded %d sources from model file %s"""%(len(model.sources),filename),3000);
        self.model = model;
        self.filename = filename;
        self.setCaption("SkyPuss - %s"%os.path.basename(self.filename));
        self._model_updated = False;
        self.emit(PYSIGNAL("isUpdated()"),(False,));
    self.emit(PYSIGNAL("hasContent()"),(True,));
    self.lv.setModel(model);
    self.proptab.setModel(model);
    self.skyplot.setModel(model);
    # add items to View menu
    self.view_menu.clear();
    self.lv.addColumnViewActionsTo(self.view_menu);
    
  def _closeFile (self):
    pass;
    
  def _saveFile (self,filename=None):
    filename = str(filename) or self.filename;
    if filename is None:
      self._saveFileAs();
    else:
      busy = BusyIndicator();
      try:
        ModelHTML.saveModel(filename,self.model);
      except:
          busy = None;
          QMessageBox.warning(self,"Error saving sky model",
              """Error saving model file %s: %s"""%(filename,str(sys.exc_info()[1])),QMessageBox.Ok,0,0);
          return False;
      self.setCaption("SkyPuss - %s"%os.path.basename(filename));
      self.filename = filename;
      self._model_updated = False;
      self.emit(PYSIGNAL("isUpdated()"),(False,));
    
  def saveFileAs (self,filename=None):
    if filename is None:
      if not self._save_as_dialog:
          dialog = self._save_as_dialog = QFileDialog(".","*.html",self,None,True);
          dialog.setMode(QFileDialog.AnyFile);
          dialog.setCaption("Save Sky Model");
          QObject.connect(dialog,SIGNAL("fileSelected(const QString&)"),self.saveFileAs);
      self._save_as_dialog.rereadDir();
      self._save_as_dialog.exec_loop();
      return;
    # filename supplied, so save
    self._saveFile(filename);
    
  def _saveSelectionAs (self):
    pass;
    
  def _exportKarmaAnnotations (self):
    pass;
    
  def _addTagToSelection (self):
    if not hasattr(self,'_add_tag_dialog'):
      self._add_tag_dialog = Widgets.AddTagDialog(self,modal=True);
    self._add_tag_dialog.setTags(self.model.tagnames);
    self._add_tag_dialog.setValue(True);
    if self._add_tag_dialog.exec_loop() != QDialog.Accepted:
      return;
    tagname,value = self._add_tag_dialog.getTag();
    if tagname is None or value is None:
      return None;
    # tag selected sources
    for src in self.model.sources:
      if src.selected:
        src.setAttribute(tagname,value);
    # add tag to model, and update property table if it is new
    if self.model.addTag(tagname):
      self.model.getTagProperty(tagname).computeTotal(self.model.sources);
      self.proptab.setModel(self.model);
      self.skyplot.setModel(self.model);
    # update listview
    self.lv.refreshSelectedItems();
    self._modelUpdated();
    
  def _removeTagsFromSelection (self):
    if not hasattr(self,'_remove_tag_dialog'):
      self._remove_tag_dialog = Widgets.SelectTagsDialog(self,modal=True,caption="Remove Tags",ok_button="Remove");
    # get set of all tags in selected sources
    tags = sets.Set();
    for src in self.model.sources:
      if src.selected:
        tags.update(src.getTagNames());
    if not tags:
      return;
    tags = list(tags);
    tags.sort();
    # show dialog
    self._remove_tag_dialog.setTags(tags);
    if self._remove_tag_dialog.exec_loop() != QDialog.Accepted:
      return;
    tags = self._remove_tag_dialog.getSelectedTags();
    if not tags:
      return;
    # ask for confirmation
    plural = (len(tags)>1 and "s") or "";
    if QMessageDialog.question(self,"Removing Tags","<P>Really remove the tag%s %s from selected sources?</P>"%(plural,", ".join(tags)),
        QMessageBox.Yes,QmessageBox.Cancel) != QMessageBox.Yes:
      return;
    # remove the tags
    for src in self.model.sources:
      if src.selected:
        for tag in tags:
          src.removeAttribute(tag);
    # reset model
    self.model.scanTags();
    self.model.initProperties();
    self.proptab.setModel(self.model);
    self.skyplot.setModel(self.model);
    # update listview
    self.lv.refreshSelectedItems();
    self._modelUpdated();
    
    
    
  def _modelUpdated (self):
    self._model_updated = True;
    self.emit(PYSIGNAL("isUpdated()"),(True,));
    pass;
