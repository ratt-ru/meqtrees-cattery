#!/usr/bin/python

###################################
### Class to select file export
### options
###################################


#% $Id$ 

#
# Copyright (C) 2002-2007
# ASTRON (Netherlands Foundation for Research in Astronomy)
# and The MeqTree Foundation
# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands, seg@astron.nl
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
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

from qt import *
import sys
from common_utils import *

class ExportDialog(QDialog):
    def __init__(self,parent = None,name = None,modal = 0,fl = 0,main_window=None):
        QDialog.__init__(self,parent,name,modal,fl)

        if not name:
            self.setName("ExportDialog")
        self.export_flag=EXPORT_IMG_KARMA
    
        self.main=main_window

        FormLayout = QVBoxLayout(self,11,6,"FormLayout")
#######################################
        self.selectMainBG = QButtonGroup(self,"selectMainBG")
        self.selectMainBG.setColumnLayout(0,Qt.Vertical)
        self.selectMainBG.layout().setSpacing(6)
        self.selectMainBG.layout().setMargin(11)
        selectMainBGLayout = QVBoxLayout(self.selectMainBG.layout())
        selectMainBGLayout.setAlignment(Qt.AlignTop)

        mainlayout = QVBoxLayout(None,0,6,"mainlayout")

        self.imageButton = QRadioButton(self.selectMainBG,"imageButton")
        mainlayout.addWidget(self.imageButton)

        self.ptableButton = QRadioButton(self.selectMainBG,"ptableButton")
        mainlayout.addWidget(self.ptableButton)

        self.stableButton = QRadioButton(self.selectMainBG,"stableButton")
        mainlayout.addWidget(self.stableButton)
        selectMainBGLayout.addLayout(mainlayout)
        FormLayout.addWidget(self.selectMainBG)

        self.imageButton.setChecked(1)
########################################
        self.imageBG = QButtonGroup(self,"imageBG")
        self.imageBG.setColumnLayout(0,Qt.Vertical)
        self.imageBG.layout().setSpacing(6)
        self.imageBG.layout().setMargin(11)
        imageBGLayout = QHBoxLayout(self.imageBG.layout())
        imageBGLayout.setAlignment(Qt.AlignTop)

        imlayout = QVBoxLayout(None,0,6,"imlayout")

        self.imKarmaButton = QRadioButton(self.imageBG,"imKarmaButton")
        imlayout.addWidget(self.imKarmaButton)
        self.imEPSButton = QRadioButton(self.imageBG,"imEPSButton")
        imlayout.addWidget(self.imEPSButton)
        self.imPNGButton = QRadioButton(self.imageBG,"imPNGButton")
        imlayout.addWidget(self.imPNGButton)
        self.imBMPButton = QRadioButton(self.imageBG,"imBMPButton")
        imlayout.addWidget(self.imBMPButton)
 
        imageBGLayout.addLayout(imlayout)
        FormLayout.addWidget(self.imageBG)

        self.imKarmaButton.setChecked(1)
        self.imageBG.show()
########################################
        self.ptableBG = QButtonGroup(self,"ptableBG")
        self.ptableBG.setColumnLayout(0,Qt.Vertical)
        self.ptableBG.layout().setSpacing(6)
        self.ptableBG.layout().setMargin(11)
        ptableBGLayout = QHBoxLayout(self.ptableBG.layout())
        ptableBGLayout.setAlignment(Qt.AlignTop)

        ptlayout = QVBoxLayout(None,0,6,"ptlayout")

        self.ptEPSButton = QRadioButton(self.ptableBG,"ptEPSButton")
        ptlayout.addWidget(self.ptEPSButton)

        self.ptTEXButton = QRadioButton(self.ptableBG,"ptTEXButton")
        ptlayout.addWidget(self.ptTEXButton)
        ptableBGLayout.addLayout(ptlayout)
        FormLayout.addWidget(self.ptableBG)

        self.ptEPSButton.setChecked(1)
        self.ptableBG.hide()
########################################
        self.stableBG = QButtonGroup(self,"stableBG")
        self.stableBG.setColumnLayout(0,Qt.Vertical)
        self.stableBG.layout().setSpacing(6)
        self.stableBG.layout().setMargin(11)
        stableBGLayout = QHBoxLayout(self.stableBG.layout())
        stableBGLayout.setAlignment(Qt.AlignTop)

        stlayout = QVBoxLayout(None,0,6,"stlayout")

        self.stEPSButton = QRadioButton(self.stableBG,"stEPSButton")
        stlayout.addWidget(self.stEPSButton)

        self.stTEXButton = QRadioButton(self.stableBG,"stTEXButton")
        stlayout.addWidget(self.stTEXButton)
        stableBGLayout.addLayout(stlayout)
        FormLayout.addWidget(self.stableBG)

        self.stEPSButton.setChecked(1)
        self.stableBG.hide()
############## Button Bar
        buttonLayout=QHBoxLayout(None,0,6,"buttonLayout")
        self.buttonHelp = QPushButton(self,"buttonHelp")
        self.buttonHelp.setAutoDefault(1)
        buttonLayout.addWidget(self.buttonHelp)

        Horizontal_Spacing = QSpacerItem(246,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
        buttonLayout.addItem(Horizontal_Spacing)

        self.buttonOk = QPushButton(self,"buttonOk")
        self.buttonOk.setAutoDefault(1)
        self.buttonOk.setDefault(1)
        buttonLayout.addWidget(self.buttonOk)

        self.buttonCancel = QPushButton(self,"buttonCancel")
        self.buttonCancel.setAutoDefault(1)
        buttonLayout.addWidget(self.buttonCancel)



        self.connect(self.buttonHelp,SIGNAL("clicked()"),self.help)
        self.connect(self.buttonOk,SIGNAL("clicked()"),self.accept)
        self.connect(self.buttonCancel,SIGNAL("clicked()"),self.reject)
        FormLayout.addLayout(buttonLayout)
########################################
        self.languageChange()

        #self.resize(QSize(600,480).expandedTo(self.minimumSizeHint()))
        self.clearWState(Qt.WState_Polished)

        self.connect(self.selectMainBG, SIGNAL("clicked(int)"), self.showGroups)
        self.connect(self.imageBG, SIGNAL("clicked(int)"), self.changeIMG)
        self.connect(self.ptableBG, SIGNAL("clicked(int)"), self.changePT)
        self.connect(self.stableBG, SIGNAL("clicked(int)"), self.changeST)

    def showGroups(self,id):
     if id==0:
        self.imageBG.show()
        self.ptableBG.hide()
        self.stableBG.hide()
     elif id==1:
        self.ptableBG.show()
        self.imageBG.hide()
        self.stableBG.hide()
     else:
        self.stableBG.show()
        self.imageBG.hide()
        self.ptableBG.hide()

    def changeIMG(self,id):
     if id==0:
      self.export_flag=EXPORT_IMG_KARMA
     elif id==1:
      self.export_flag=EXPORT_IMG_EPS
     elif id==2:
      self.export_flag=EXPORT_IMG_PNG
     elif id==3:
      self.export_flag=EXPORT_IMG_BMP

    def changePT(self,id):
     if id==0:
      self.export_flag=EXPORT_PT_EPS
     elif id==1:
      self.export_flag=EXPORT_PT_TEX
     elif id==2:
      self.export_flag=EXPORT_PT_TXT

    def changeST(self,id):
     if id==0:
      self.export_flag=EXPORT_ST_EPS
     elif id==1:
      self.export_flag=EXPORT_ST_TEX
     elif id==2:
      self.export_flag=EXPORT_ST_TXT


    def help(self):
     print "Not yet done"


    def accept(self):
     QDialog.accept(self)
     caption_str="Save File"
     if self.export_flag==EXPORT_IMG_KARMA:
      info_str="Choose Karma Annotation File"
      ext_str="*.ann"
     elif self.export_flag==EXPORT_IMG_EPS:
      info_str="Choose EPS File Name"
      ext_str="*.eps"
     elif self.export_flag==EXPORT_IMG_PNG:
      info_str="Choose PNG File Name"
      ext_str="*.png"
     elif self.export_flag==EXPORT_IMG_BMP:
      info_str="Choose BMP File Name"
      ext_str="*.bmp"
     elif self.export_flag==EXPORT_PT_EPS:
      info_str="Choose EPS File Name"
      ext_str="*.eps"
     elif self.export_flag==EXPORT_PT_TEX:
      info_str="Choose TEX File Name"
      ext_str="*.tex"
     elif self.export_flag==EXPORT_PT_TXT:
      info_str="Choose Text File Name"
      ext_str="*.txt"
     elif self.export_flag==EXPORT_ST_EPS:
      info_str="Choose EPS File Name"
      ext_str="*.eps"
     elif self.export_flag==EXPORT_ST_TEX:
      info_str="Choose TEX File Name"
      ext_str="*.tex"
     elif self.export_flag==EXPORT_ST_TXT:
      info_str="Choose Text File Name"
      ext_str="*.txt"
     else:
      return 
     s=None
     s=QFileDialog.getSaveFileName(".",ext_str,self,caption_str,info_str)
     if s==None:
      return
     if self.export_flag==EXPORT_IMG_KARMA:
      if not s.endsWith(".ann"):
       s.append(".ann")
      self.main.lsm.export_karma_annotations(s.ascii())
     elif self.export_flag==EXPORT_IMG_EPS:
      if not s.endsWith(".eps"):
       s.append(".eps")
      self.main.filePrintEPS(s.ascii())
     elif self.export_flag==EXPORT_IMG_PNG:
      if not s.endsWith(".png"):
       s.append(".png")
      self.main.fileExportImage(s.ascii(),"PNG")
     elif self.export_flag==EXPORT_IMG_BMP:
      if not s.endsWith(".bmp"):
       s.append(".bmp")
      self.main.fileExportImage(s.ascii(),"BMP")
     elif self.export_flag==EXPORT_PT_EPS:
      if not s.endsWith(".eps"):
       s.append(".eps")
     elif self.export_flag==EXPORT_PT_TEX:
      if not s.endsWith(".tex"):
       s.append(".tex")
      self.main.exportPUTableTEX(s.ascii())
     elif self.export_flag==EXPORT_PT_TXT:
      if not s.endsWith(".txt"):
       s.append(".txt")
     elif self.export_flag==EXPORT_ST_EPS:
      if not s.endsWith(".eps"):
       s.append(".eps")
     elif self.export_flag==EXPORT_ST_TEX:
      if not s.endsWith(".tex"):
       s.append(".tex")
      self.main.exportSTableTEX(s.ascii())
     elif self.export_flag==EXPORT_ST_TXT:
      if not s.endsWith(".txt"):
       s.append(".txt")
 

    def reject(self):
     print "reject options"
     QDialog.reject(self)

    def languageChange(self):
        self.setCaption(self.__tr("Export As"))
        self.selectMainBG.setTitle(self.__tr("Select..."))
        self.imageButton.setText(self.__tr("Export Image"))
        self.ptableButton.setText(self.__tr("Export PUnit Table"))
        self.stableButton.setText(self.__tr("Export Source Table"))
        self.imageBG.setTitle(self.__tr("Export Image as:"))
        self.imKarmaButton.setText(self.__tr("Export Karma Annotations"))
        self.imEPSButton.setText(self.__tr("Export Image as EPS Image"))
        self.imPNGButton.setText(self.__tr("Export Image as PNG Image"))
        self.imBMPButton.setText(self.__tr("Export Image as BMP Image"))
        self.ptableBG.setTitle(self.__tr("Export PUnit Table as:"))
        self.ptEPSButton.setText(self.__tr("Export PUnit Table as EPS File"))
        self.ptTEXButton.setText(self.__tr("Export PUnit Table as LaTex Source"))
        self.stableBG.setTitle(self.__tr("Export Source Table as:"))
        self.stEPSButton.setText(self.__tr("Export Source Table as EPS File"))
        self.stTEXButton.setText(self.__tr("Export Source Table as LaTex Source"))
        self.buttonHelp.setText(self.__tr("&Help"))
        self.buttonHelp.setAccel(self.__tr("F1"))
        self.buttonOk.setText(self.__tr("&OK"))
        self.buttonOk.setAccel(QString.null)
        self.buttonCancel.setText(self.__tr("&Cancel"))
        self.buttonCancel.setAccel(QString.null)



    def __tr(self,s,c = None):
        return qApp.translate("ExportDialog",s,c)

#############################################################################
def main(args):
  app=QApplication(args)
  win=ExportDialog(None)
  #win.setTitle("Testing...")
  win.show()
  app.connect(app,SIGNAL("lastWindowClosed()"),
               app,SLOT("quit()"))
  app.exec_loop()

if __name__=="__main__":
   main(sys.argv)

