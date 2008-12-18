
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

import sys
from qt import *


class PatchOptionsDialog(QDialog):
    def __init__(self,parent = None,name = None,modal = 0,fl = 0,canvas_view=None,lsm=None):
        QDialog.__init__(self,parent,name,modal,fl)

        if not name:
            self.setName("PatchOptionsDialog")

        self.minS=2 # there should be at least 2 sources for a patch
        if lsm!=None:
         self.minBval=lsm.getMinBrightness('A')
         self.maxBval=lsm.getMaxBrightness('A')
         self.maxS=len(lsm.s_table)
        else:
         self.minBval=0.0
         self.maxBval=10.0
         self.maxS=100
   
        self.cview=canvas_view

        self.setSizeGripEnabled(1)

        PatchOptionsDialogLayout = QHBoxLayout(self,11,6,"PatchOptionsDialogLayout")

        layout4 = QVBoxLayout(None,0,6,"layout4")

        layout2 = QHBoxLayout(None,0,6,"layout2")

        self.minSLabel = QLabel(self,"minSLabel")
        layout2.addWidget(self.minSLabel)

        self.spinBox_minSources = QSpinBox(self,"spinBox_min")
        self.spinBox_minSources.setMaxValue(self.maxS)
        self.spinBox_minSources.setLineStep(1)

        layout2.addWidget(self.spinBox_minSources)
        layout4.addLayout(layout2)

        layout3 = QHBoxLayout(None,0,6,"layout3")

        self.BrightnessLabel = QLabel(self,"BrightnessLabel")
        layout3.addWidget(self.BrightnessLabel)
        spacer3 = QSpacerItem(0,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
        layout3.addItem(spacer3)

        self.minBLabel = QLabel(self,"minBLabel")
        layout3.addWidget(self.minBLabel)

        self.ledit_minB=QLineEdit(self,"lineEdit_minB")
        self.ledit_minB.setText(str(self.minBval))
        layout3.addWidget(self.ledit_minB)

        self.maxBLabel = QLabel(self,"maxBLabel")
        layout3.addWidget(self.maxBLabel)

        self.ledit_maxB=QLineEdit(self,"lineEdit_maxB")
        self.ledit_maxB.setText(str(self.maxBval))
        layout3.addWidget(self.ledit_maxB)
        layout4.addLayout(layout3)

        Layout1 = QHBoxLayout(None,0,6,"Layout1")

        self.buttonHelp = QPushButton(self,"buttonHelp")
        self.buttonHelp.setAutoDefault(1)
        Layout1.addWidget(self.buttonHelp)
        Horizontal_Spacing2 = QSpacerItem(0,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
        Layout1.addItem(Horizontal_Spacing2)

        self.buttonOk = QPushButton(self,"buttonOk")
        self.buttonOk.setAutoDefault(1)
        self.buttonOk.setDefault(1)
        Layout1.addWidget(self.buttonOk)

        self.buttonCancel = QPushButton(self,"buttonCancel")
        self.buttonCancel.setAutoDefault(1)
        Layout1.addWidget(self.buttonCancel)
        layout4.addLayout(Layout1)
        PatchOptionsDialogLayout.addLayout(layout4)

        self.languageChange()

        #self.resize(QSize(511,282).expandedTo(self.minimumSizeHint()))
        self.clearWState(Qt.WState_Polished)

        self.connect(self.buttonOk,SIGNAL("clicked()"),self.accept)
        self.connect(self.buttonCancel,SIGNAL("clicked()"),self.reject)
        self.connect(self.ledit_minB,SIGNAL("textChanged(const QString &)"),self.minBright)
        self.connect(self.ledit_maxB,SIGNAL("textChanged(const QString &)"),self.maxBright)

        self.connect(self.spinBox_minSources,SIGNAL("valueChanged(int)"),self.spinBox_min_valueChanged)


    def languageChange(self):
        self.setCaption(self.__tr("Options For Patches"))
        self.minSLabel.setText(self.__tr("Minimum Number of Sources "))
        self.BrightnessLabel.setText(self.__tr("Apparent Brightness :"))
        self.minBLabel.setText(self.__tr("Min :"))
        self.spinBox_minSources.setSpecialValueText(self.__tr(str(self.minS),"All"))
        QToolTip.add(self.spinBox_minSources,self.__tr("Minimum Number of Sources"))
        self.maxBLabel.setText(self.__tr("Max :"))
        self.buttonHelp.setText(self.__tr("&Help"))
        self.buttonHelp.setAccel(self.__tr("F1"))
        self.buttonOk.setText(self.__tr("&OK"))
        self.buttonOk.setAccel(QString.null)
        self.buttonCancel.setText(self.__tr("&Cancel"))
        self.buttonCancel.setAccel(QString.null)


    def minBright(self,newText):
     if newText.length()>0:
      try:
       self.minBval=float(newText.ascii())
      except (TypeError,ValueError):
       self.ledit_minB.setText(str(self.minBval))

    def maxBright(self,newText):
     if newText.length()>0:
      try:
       self.maxBval=float(newText.ascii())
      except (TypeError,ValueError):
       self.ledit_maxB.setText(str(self.maxBval))




    def spinBox_min_valueChanged(self,a0):
      self.minS=a0

    def __tr(self,s,c = None):
        return qApp.translate("PatchOptionsDialog",s,c)

    def accept(self):
     if self.cview !=None:
      self.cview.createPatchesFromGrid(self.minBval,self.maxBval,self.minS)
     QDialog.accept(self)

    def reject(self):
     QDialog.reject(self)

############################################################
def main(args):
  app=QApplication(args)
  win=PatchOptionsDialog(None)
  win.show()
  app.connect(app,SIGNAL("lastWindowClosed()"),
               app,SLOT("quit()"))
  app.exec_loop()

if __name__=="__main__":
   main(sys.argv)
