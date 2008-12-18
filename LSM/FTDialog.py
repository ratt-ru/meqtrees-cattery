#!/usr/bin/python

########################################
#### Change Freq. and Time
########################################

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

from common_utils import *
import sys

class FTDialog(QDialog):
    def __init__(self,parent = None,name = None,modal = 0,fl = 0,**kw):
        QDialog.__init__(self,parent,name,modal,fl)

        if not name:
            self.setName("Form1")
        ## initialize the range controls
        if kw.has_key('f0'):
         self.f0=kw['f0']
        else:
         self.f0=0.0
        if kw.has_key('f1'):
         self.f1=kw['f1']
        else:
         self.f1=1.0
        if kw.has_key('fticks'):
         self.fticks=kw['fticks']
        else:
         self.fticks=2
        if kw.has_key('f'):
         self.f=kw['f']
        else:
         self.f=0


        if kw.has_key('t0'):
         self.t0=kw['t0']
        else:
         self.t0=0.0
        if kw.has_key('t1'):
         self.t1=kw['t1']
        else:
         self.t1=1.0
        if kw.has_key('tticks'):
         self.tticks=kw['tticks']
        else:
         self.tticks=2
        if kw.has_key('t'):
         self.t=kw['t']
        else:
         self.t=0



        FormLayout = QHBoxLayout(self,11,6,"FormLayout")

        layout8 = QVBoxLayout(None,0,6,"layout8")

        layout3 = QHBoxLayout(None,0,6,"layout3")

        self.textLabelFreq = QLabel(self,"textLabelFreq")
        self.textLabelFreq.setFixedSize(100,20)
        layout3.addWidget(self.textLabelFreq)


        self.sliderF = QSlider(0,self.fticks,1,self.f, QSlider.Horizontal,self,"sliderF")
        self.sliderF.setFixedSize(200,20)
        self.sliderF.setTickmarks(1)
        self.connect( self.sliderF, SIGNAL("valueChanged( int )"), self.updateFreq )
        layout3.addWidget(self.sliderF)
        spacer1 = QSpacerItem(11,11,QSizePolicy.Expanding,QSizePolicy.Minimum)
        layout3.addItem(spacer1)


        self.textLabelF = QLabel(self,"textLabelF")
        layout3.addWidget(self.textLabelF)
        layout8.addLayout(layout3)

        layout4 = QHBoxLayout(None,0,6,"layout4")

        self.textLabelTime = QLabel(self,"textLabelTime")
        self.textLabelTime.setFixedSize(100,20)
        layout4.addWidget(self.textLabelTime)


        self.sliderT = QSlider(0,self.tticks,1,self.t, QSlider.Horizontal,self,"sliderT")
        self.sliderT.setFixedSize(200,20)
        self.sliderT.setTickmarks(1)
        self.connect( self.sliderT, SIGNAL("valueChanged( int )"), self.updateTime )
        layout4.addWidget(self.sliderT)
        spacer2 = QSpacerItem(11,11,QSizePolicy.Expanding,QSizePolicy.Minimum)
        layout4.addItem(spacer2)


        self.textLabelT = QLabel(self,"textLabelT")
        layout4.addWidget(self.textLabelT)
        layout8.addLayout(layout4)

        layout7 = QHBoxLayout(None,0,6,"layout7")

        self.pushButton1 = QPushButton(self,"pushButton1")
        layout7.addWidget(self.pushButton1)
        self.connect(self.pushButton1,SIGNAL("clicked()"),self.accept)

        self.pushButton2 = QPushButton(self,"pushButton2")
        layout7.addWidget(self.pushButton2)
        self.connect(self.pushButton2,SIGNAL("clicked()"),self.reject)
        layout8.addLayout(layout7)
        FormLayout.addLayout(layout8)

        self.languageChange()

        #self.resize(QSize(600,480).expandedTo(self.minimumSizeHint()))
        self.clearWState(Qt.WState_Polished)


    def languageChange(self):
        self.setCaption(self.__tr("Change F,T"))
        self.textLabelFreq.setText(self.__tr("Frequency :"))
        fval=self.f0+(self.f1-self.f0)/self.fticks*self.f
        tmpval=stdForm(fval,'%3.4f')
        self.textLabelF.setText(tmpval[0]+tmpval[1]+'Hz') 

        self.textLabelTime.setText(self.__tr("Time         :"))
        tval=self.t0+(self.t1-self.t0)/self.tticks*self.t
        tmpval=stdForm(tval,'%3.4f')
        #self.textLabelT.setText('%5.4f'%tval) 
        self.textLabelT.setText(tmpval[0]+tmpval[1]+'s') 

        self.pushButton1.setText(self.__tr("OK"))
        self.pushButton2.setText(self.__tr("Cancel"))


    def __tr(self,s,c = None):
        return qApp.translate("FTDialog",s,c)

    def updateFreq(self,intval):
      fval=self.f0+(self.f1-self.f0)/self.fticks*intval
      tmpval=stdForm(fval,'%3.4f')
      self.textLabelF.setText(tmpval[0]+tmpval[1]+'Hz') 

    def updateTime(self,intval):
      tval=self.t0+(self.t1-self.t0)/self.tticks*intval
      tmpval=stdForm(tval,'%3.4f')
      self.textLabelT.setText(tmpval[0]+tmpval[1]+'s') 

    def accept(self):
     self.parentWidget().cview.updateDisplay(self.parentWidget().cview.default_mode,self.sliderF.value(),self.sliderT.value())
     QDialog.accept(self)

    def reject(self):
     QDialog.accept(self)

def main(args):
  app=QApplication(args)
  win=FTDialog(None,f0=100,f1=1000,fticks=5,f=3,t0=0,t1=120,tticks=4,t=2)
  win.show()
  app.connect(app,SIGNAL("lastWindowClosed()"),
               app,SLOT("quit()"))
  app.exec_loop()

if __name__=="__main__":
   main(sys.argv)

