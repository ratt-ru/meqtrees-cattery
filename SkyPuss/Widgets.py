import sys
import math
import traceback

from qt import *

class ValueTypeEditor (QWidget):
  ValueTypes = (bool,int,float,complex,str);
  def __init__ (self,*args):
    QWidget.__init__(self,*args);
    lo = QHBoxLayout(self);
    lo.setSpacing(5);
    # type selector
    self.wtypesel = QComboBox(self);
    for i,tp in enumerate(self.ValueTypes):
      self.wtypesel.insertItem(tp.__name__,i);
    QObject.connect(self.wtypesel,SIGNAL("activated(int)"),self._selectTypeNum);
    typesel_lab = QLabel(self.wtypesel,"&Type:",self);
    lo.addWidget(typesel_lab,0);
    lo.addWidget(self.wtypesel,0);
    self.wvalue = QLineEdit(self);
    self.wvalue_lab = QLabel(self.wvalue,"&Value:",self);
    self.wbool = QCheckBox("set",self);
    self.wbool.setChecked(True);
    lo.addWidget(self.wvalue_lab,0);
    lo.addWidget(self.wvalue,1);
    lo.addWidget(self.wbool,1);
    self.wvalue_lab.hide();
    self.wvalue.hide();
    # make input validators
    self._validators = {int:QIntValidator(self),float:QDoubleValidator(self) };
    # select bool type initially
    self._selectTypeNum(0);

  def _selectTypeNum (self,index):
    tp = self.ValueTypes[index];
    self.wbool.setShown(tp is bool);
    self.wvalue.setShown(tp is not bool);
    self.wvalue_lab.setShown(tp is not bool);
    self.wvalue.setValidator(self._validators.get(tp,None));
    
  def setValue (self,value):
    """Sets current value""";
    for i,tp in enumerate(self.ValueTypes):
      if isinstance(value,tp):
        self.wtypesel.setCurrentItem(i);
        self._selectTypeNum(i);
        if tp is bool:
          self.wbool.setChecked(value);
        else:
          self.wvalue.setText(str(value));
        return;
    # unknown value: set bool
    self.setValue(True);
      
  def getValue (self):
    """Returns current value, or None if no legal value is set""";
    tp = self.ValueTypes[self.wtypesel.currentItem()];
    if tp is bool:
      return self.wbool.isChecked();
    else:
      try:
        return tp(self.wvalue.text());
      except:
        print "Error converting input to type ",tp.__name__;
        traceback.print_exc();
        return None;
  
class AddTagDialog (QDialog):
  def __init__ (self,parent,name=None,modal=True,flags=0):
    QDialog.__init__(self,parent,name,True,flags);
    self.setCaption("Add Tag");
    lo = QVBoxLayout(self);
    lo.setMargin(10);
    lo.setSpacing(5);
    # tag selector
    lo1 = QHBoxLayout(lo);
    lo1.setSpacing(5);
    self.wtagsel = QComboBox(True,self);
    self.wtagsel.setCurrentText("");
    wtagsel_lbl = QLabel(self.wtagsel,"&Tag:",self);
    lo1.addWidget(wtagsel_lbl,0);
    lo1.addWidget(self.wtagsel,1);
    QObject.connect(self.wtagsel,SIGNAL("activated(int)"),self._check_tag);
    QObject.connect(self.wtagsel,SIGNAL("textChanged(const QString &)"),self._check_tag_text);
    # value editor
    self.valedit = ValueTypeEditor(self);
    lo.addWidget(self.valedit);
    # buttons
    lo.addSpacing(10);
    lo2 = QHBoxLayout(lo);
    lo2.setMargin(5);
    self.wokbtn = QPushButton("OK",self);
    self.wokbtn.setMinimumWidth(128);
    QObject.connect(self.wokbtn,SIGNAL("clicked()"),self.accept);
    self.wokbtn.setEnabled(False);
    cancelbtn = QPushButton("Cancel",self);
    cancelbtn.setMinimumWidth(128);
    QObject.connect(cancelbtn,SIGNAL("clicked()"),self.reject);
    lo2.addWidget(self.wokbtn);
    lo2.addStretch(1);
    lo2.addWidget(cancelbtn);
    self.setMinimumWidth(384);
    
  def setTags (self,tagnames):
    while self.wtagsel.count():
      self.wtagsel.removeItem(0);
    for tag in tagnames:
      self.wtagsel.insertItem(tag);
    self.wtagsel.setCurrentText("");
  
  def setValue (self,value):
    self.valedit.setValue(value);
  
  def _check_tag (self,tag):
    self.wokbtn.setEnabled(True);
    
  def _check_tag_text (self,text):
    self.wokbtn.setEnabled(bool(str(text)!=""));
    
  def getTag (self):
    return str(self.wtagsel.currentText()),self.valedit.getValue();
    
class SelectTagsDialog (QDialog):
  def __init__ (self,parent,name=None,modal=True,flags=0,caption="Select Tags",ok_button="Select"):
    QDialog.__init__(self,parent,name,True,flags);
    self.setCaption(caption);
    lo = QVBoxLayout(self);
    lo.setMargin(10);
    lo.setSpacing(5);
    # tag selector
    self.wtagsel = QListBox(self);
    lo.addWidget(self.wtagsel);
    self.wtagsel.setColumnMode(QListBox.FitToWidth);
    self.wtagsel.setSelectionMode(QListBox.Multi);
    QObject.connect(self.wtagsel,SIGNAL("selectionChanged()"),self._check_tag);
    # buttons
    lo.addSpacing(10);
    lo2 = QHBoxLayout(lo);
    lo2.setMargin(5);
    self.wokbtn = QPushButton(ok_button,self);
    self.wokbtn.setMinimumWidth(128);
    QObject.connect(self.wokbtn,SIGNAL("clicked()"),self.accept);
    self.wokbtn.setEnabled(False);
    cancelbtn = QPushButton("Cancel",self);
    cancelbtn.setMinimumWidth(128);
    QObject.connect(cancelbtn,SIGNAL("clicked()"),self.reject);
    lo2.addWidget(self.wokbtn);
    lo2.addStretch(1);
    lo2.addWidget(cancelbtn);
    self.setMinimumWidth(384);
    self._tagnames = [];
    
  def setTags (self,tagnames):
    self._tagnames = tagnames;
    self.wtagsel.clear();
    for tag in tagnames:
      self.wtagsel.insertItem(tag);
      
  def _check_tag (self):
    for i in range(len(self._tagnames)):
      if self.wtagsel.isSelected(i):
        self.wokbtn.setEnabled(True);
        return;
    else:
      self.wokbtn.setEnabled(False);
    
  def getSelectedTags (self):
    return [ tag for i,tag in enumerate(self._tagnames) if self.wtagsel.isSelected(i) ];
    
