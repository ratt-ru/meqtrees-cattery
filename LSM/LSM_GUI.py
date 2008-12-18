#!/usr/bin/python

#########################################################
### The main Window and the Tabs
#########################################################


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
import math
from copy import deepcopy
from qt import *
from qttable import QTable
from qtcanvas import *

from LSM import *
from common_utils import *

from OptionsDialog import *
from FTDialog import *
from MyCanvasView import *
from TreeDisp import *
from ExportDialog import *
from PatchOptionsDialog import *
from TransDialog import *

image0_data = \
    "\x89\x50\x4e\x47\x0d\x0a\x1a\x0a\x00\x00\x00\x0d" \
    "\x49\x48\x44\x52\x00\x00\x00\x16\x00\x00\x00\x16" \
    "\x08\x06\x00\x00\x00\xc4\xb4\x6c\x3b\x00\x00\x00" \
    "\x77\x49\x44\x41\x54\x38\x8d\x63\x60\xa0\x11\x60" \
    "\xc4\x22\xf6\x9f\x02\xbd\x70\xc0\x84\x4d\xf0\xff" \
    "\xff\xff\x38\x31\x03\x03\x03\x83\xa6\x8e\x26\x41" \
    "\x07\x60\x35\x98\x10\x30\x34\x34\x24\x68\x38\x0b" \
    "\x39\x06\x6b\x6a\x68\xc2\xd9\xd7\xaf\x5c\xff\xcf" \
    "\x80\x25\x58\xc8\x32\xd8\xc9\xd6\x09\x85\x8f\xcd" \
    "\x70\xb2\x0c\xb6\xb2\xb5\xc2\x10\x43\x37\x9c\x2c" \
    "\x83\x19\x19\xf1\x26\x08\xf2\x0c\x86\xa5\x0c\x42" \
    "\x96\x91\x95\x2a\x88\x01\xa3\x06\x8f\x1a\x3c\x6a" \
    "\x30\x1e\x40\x49\x0d\x42\xac\x79\x43\x04\x00\x00" \
    "\x25\x19\x2a\x1c\x1c\xcc\x28\x4a\x00\x00\x00\x00" \
    "\x49\x45\x4e\x44\xae\x42\x60\x82"
image1_data = \
    "\x89\x50\x4e\x47\x0d\x0a\x1a\x0a\x00\x00\x00\x0d" \
    "\x49\x48\x44\x52\x00\x00\x00\x16\x00\x00\x00\x16" \
    "\x08\x06\x00\x00\x00\xc4\xb4\x6c\x3b\x00\x00\x00" \
    "\x99\x49\x44\x41\x54\x38\x8d\xed\x94\x41\x0e\x85" \
    "\x20\x0c\x44\x5f\x89\xd7\x36\x3f\xc6\x83\xd7\x85" \
    "\x82\x44\x8a\x02\xf6\xef\x9c\x4d\xb5\xe8\xeb\x38" \
    "\x41\xe0\xd3\x21\x19\x78\x47\x1f\x18\x3a\xc2\xbd" \
    "\x42\x63\x4f\xaf\xd7\x23\x8e\x5b\x86\x4a\xf8\x03" \
    "\x14\x80\xa9\xb2\xd0\xf3\x25\x56\x3c\x08\xa0\xaa" \
    "\x71\x5d\x00\x45\xa4\x99\x5b\x7d\x70\x3a\x87\x4a" \
    "\xaa\xfb\x24\x29\xfa\x3d\xc3\x43\x66\xbc\xb3\xde" \
    "\x2b\x58\x8e\xdb\xea\xbd\xcc\x8c\x07\xb2\x2e\x64" \
    "\x66\x9c\x43\xd7\xa5\x0f\x38\xff\xf6\x6a\x66\xfc" \
    "\x16\xca\xf9\x83\xf8\x42\x0b\xc7\x5e\xd0\x0c\xec" \
    "\x0b\x8d\x37\x69\xef\x78\x41\x21\x39\xf6\x85\xc6" \
    "\xe6\xf3\x6e\xaf\xcb\xf3\xd8\x6d\xd3\x06\xc2\x20" \
    "\x45\x2c\xdf\x60\xc8\xbd\x00\x00\x00\x00\x49\x45" \
    "\x4e\x44\xae\x42\x60\x82"
image2_data = \
    "\x89\x50\x4e\x47\x0d\x0a\x1a\x0a\x00\x00\x00\x0d" \
    "\x49\x48\x44\x52\x00\x00\x00\x16\x00\x00\x00\x16" \
    "\x08\x06\x00\x00\x00\xc4\xb4\x6c\x3b\x00\x00\x00" \
    "\xa2\x49\x44\x41\x54\x38\x8d\x63\x60\xa0\x11\x60" \
    "\x84\xd2\xff\xa9\x6d\x2e\x0b\x8c\xd5\xdc\x08\xa1" \
    "\x6b\xeb\x19\x18\x8e\x1e\x3a\x4a\x92\x29\xd6\x76" \
    "\xd6\x70\x3d\xd6\x76\xd6\x0c\x0c\x0c\x0c\x0c\x4c" \
    "\x54\x72\x21\x06\x60\x21\xac\x04\x3b\xe0\xf8\x22" \
    "\xc0\xf0\x83\xe7\x03\x9c\x0f\x73\x29\x0c\x10\x74" \
    "\x31\xc7\x17\x01\x06\x8e\x2f\x02\x18\xe2\xc8\x86" \
    "\x32\x30\x40\x82\x12\x16\x9c\x44\x19\x4c\x2e\x20" \
    "\x18\x14\xe8\x2e\x23\x16\xd0\xcc\xc5\x43\xcf\x60" \
    "\x8c\x30\x6e\x6e\xc4\x4c\x3a\x84\x00\x72\x6a\xc0" \
    "\x69\x30\x2e\x85\xa4\x02\xfa\x26\xb7\xda\x7a\xd2" \
    "\x0d\x42\xf7\x25\xd9\x59\x1a\xb9\xa0\xc2\x16\x27" \
    "\x43\x2f\xb9\x0d\x3d\x83\x19\x91\xd8\xd4\xaa\x9e" \
    "\x18\x09\x2b\xa1\x00\x00\x00\x1a\x6d\x1e\xb1\xec" \
    "\x51\x25\x55\x00\x00\x00\x00\x49\x45\x4e\x44\xae" \
    "\x42\x60\x82"
image3_data = \
    "\x89\x50\x4e\x47\x0d\x0a\x1a\x0a\x00\x00\x00\x0d" \
    "\x49\x48\x44\x52\x00\x00\x00\x16\x00\x00\x00\x16" \
    "\x08\x06\x00\x00\x00\xc4\xb4\x6c\x3b\x00\x00\x02" \
    "\xb9\x49\x44\x41\x54\x38\x8d\x8d\x94\xbd\x6f\x13" \
    "\x4d\x10\xc6\x7f\x46\x2e\x66\x25\x8a\xbd\xee\x22" \
    "\x51\xd8\x1d\xee\xe2\x54\x90\x8e\xd0\x59\xa1\x09" \
    "\x69\x00\xd1\xbc\xe1\x7d\x25\x3e\xfe\x01\x24\xf3" \
    "\x17\xd0\x52\x92\x0a\x10\x05\x21\x0d\x28\x14\x08" \
    "\x28\x10\x4e\x45\x68\x5e\xe2\x74\x9b\x02\x61\x57" \
    "\xdc\x76\x37\xdd\x52\xec\xdd\xe5\xce\x1f\x82\x69" \
    "\xf6\x6e\x76\xf6\x99\x67\x66\x9e\xdd\x16\x4b\xec" \
    "\xed\xc1\xdb\x60\xc4\x60\xad\x6d\xf8\x55\xb5\xb1" \
    "\x66\x3e\x23\xb1\x09\x57\x37\xae\xb6\xea\x71\xed" \
    "\x65\xc0\x46\x0c\xfd\xd5\x7e\x05\x50\x37\x11\xa9" \
    "\xbe\xbd\xf7\xb8\x53\x37\x17\xb3\x14\x18\x60\x7c" \
    "\x32\x9e\x03\x9b\x4d\xa4\xaa\x8d\x44\x7f\x04\xce" \
    "\x35\x9f\x03\xc8\x35\xc7\x88\x69\xec\xab\x2a\x89" \
    "\x4d\xe6\xce\x9f\x5b\x04\xfa\x65\xf4\x25\xd4\x99" \
    "\x88\x08\x22\xd2\x00\x30\x62\x48\x6c\x82\x88\xcc" \
    "\xf5\x77\x29\xe3\x92\x45\xe6\x33\x8c\x98\x3f\xf6" \
    "\x79\x91\x2d\x64\x0c\xb1\xd4\x45\x25\x96\xa0\xb3" \
    "\xea\xf8\x2b\xc6\xcb\x06\x55\x5a\x1d\x74\x19\xf3" \
    "\xb9\xde\x00\xec\xed\xef\x85\x6e\xa7\xcb\x4a\xba" \
    "\xb2\x14\xbc\x1c\xe4\xe8\x70\x14\xc9\x18\xcb\xf6" \
    "\xd6\x66\x85\xd7\x60\xfc\xe1\xe3\x87\xe0\x9c\xe3" \
    "\xf0\xfd\x21\x72\x4d\xaa\x92\x17\xad\xa5\x3a\xfc" \
    "\xd4\x33\xfa\x3c\x62\xe2\x27\xdc\xbb\x7b\x2f\x0c" \
    "\xae\x0d\xd8\x1c\x6c\xb6\x5a\x75\xd0\xf1\xd8\x91" \
    "\xfd\x98\x10\x49\x2e\x66\x0a\x20\x62\x11\x0c\x9c" \
    "\x07\xf7\xc3\xf1\xed\xdb\x11\x93\xe9\x14\xd5\x9c" \
    "\xde\xc5\x2e\x3b\xff\xee\xc4\x56\x94\x4c\x37\x2e" \
    "\x6f\x20\xd6\xc6\xf2\x55\xc9\x67\x11\x35\x07\x0c" \
    "\x87\x8c\xf8\xc8\x01\xeb\xba\x8e\x1d\x5b\x76\x9f" \
    "\xed\xe2\x4e\xe2\xed\x53\x55\xba\x17\x7b\x11\xf4" \
    "\xe9\xee\xf3\xf0\xeb\xe7\xaf\xf0\x37\xf6\x3d\x7c" \
    "\x0f\xbd\xd0\x0b\x12\x24\xf4\x42\x2f\x7c\x3d\xfa" \
    "\x1a\x1e\x3f\x7c\x12\xfa\xab\xfd\x90\xa6\x9d\x60" \
    "\xad\x0d\x9d\xb4\x13\xda\xd9\x24\xe3\xe6\x8d\x2d" \
    "\x8c\x31\x4b\x4b\xaf\x5b\x42\xc2\x15\x2e\x23\x08" \
    "\x03\x06\xac\xf5\xd7\x48\x64\x05\xaf\x13\x76\x5f" \
    "\xbe\xc0\x7b\x8f\x4d\x2d\x6d\x05\xdc\x69\xec\x0f" \
    "\xaa\x20\x42\x54\x50\x4c\x64\x04\x04\x89\x1d\x2f" \
    "\xf6\xef\xeb\x90\xeb\x38\x52\x4d\x39\xe2\x18\x54" \
    "\xb1\x17\x56\xa0\x94\x21\xd0\xf6\xd3\x29\xee\xe4" \
    "\x38\x1e\xcc\x95\x6a\x6e\x52\x00\x51\xe8\x54\x14" \
    "\x54\xaa\x4d\x01\x9c\x4e\x50\x94\x5c\x15\xf7\xff" \
    "\x71\x35\x6e\xa1\x90\x5b\xe9\x38\xc3\xd1\xca\xa9" \
    "\xe8\x59\x8e\xb3\xc8\xea\x5f\x22\x9f\xe6\x8c\x51" \
    "\xda\x08\x68\xd6\x3c\x20\x34\xc5\x16\x13\xd6\x19" \
    "\x97\x26\xa8\x7a\x4a\xfd\x94\x77\x70\x6d\x75\x9d" \
    "\x73\xd4\x59\x56\x64\x0b\xc6\x5a\xac\x33\x6d\x28" \
    "\xa3\x55\x3d\x79\xe9\x11\x83\x02\x69\x9a\xf2\xcf" \
    "\x7f\x37\x23\x63\x77\xea\x8a\x70\x03\x92\x83\x37" \
    "\x40\x06\x62\x0a\xba\x79\xfc\x06\x4c\xc1\xcf\x54" \
    "\x1c\x35\xfe\x69\x8e\x08\x0c\x1f\x0d\xd9\xde\xda" \
    "\x6e\xb5\x3b\xdd\x0e\xee\xcd\x27\x5e\xbf\x3b\x40" \
    "\x8b\x4d\x30\x34\xc4\x27\x52\x95\x49\xf9\x46\x57" \
    "\xcd\xa3\x52\xcb\xf0\xd1\x90\x07\xf7\x1f\xb4\xa0" \
    "\xf6\x08\xed\xed\xef\x85\xfd\x57\x07\x4c\x9c\xab" \
    "\x06\x30\x7b\xb5\xcb\x07\x29\x57\x50\xcd\xa2\xcf" \
    "\x2b\xfd\x4b\x7d\x6e\xdd\xbe\xc5\x9d\x9d\x3b\x15" \
    "\xde\x6f\xf6\x3c\x9e\x4e\x2a\x0d\x25\x44\x00\x00" \
    "\x00\x00\x49\x45\x4e\x44\xae\x42\x60\x82"

######################################################

#######################################################

class LSMWindow(QMainWindow):
    def __init__(self,lsm_object,parent = None,name = None,fl = 0):
        QMainWindow.__init__(self,parent,name,fl)
        self.statusBar()

        #self.image0 = QPixmap()  # FILE - New
        #self.image0.loadFromData(image0_data,"PNG")
        self.image1 = QPixmap()
        self.image1.loadFromData(image1_data,"PNG")
        self.image2 = QPixmap()
        self.image2.loadFromData(image2_data,"PNG")
        self.image3 = QPixmap()
        self.image3.loadFromData(image3_data,"PNG")
        if not name:
            self.setName("LSM")

        self.lsm=lsm_object
        self.savefile=None # file to save the LSM
        if self.lsm==None:
          raise NameError, "LSM object not defined"

        self.setCentralWidget(QWidget(self,"qt_central_widget"))
        MainLayout = QVBoxLayout(self.centralWidget(),11,6,"MainLayout")

        self.tabWidget = QTabWidget(self.centralWidget(),"tabWidget")
###################### Tab 1
        self.sourceTab = QWidget(self.tabWidget,"sourceTab")

        # layout for table1
        sourceLayout=QVBoxLayout(self.sourceTab)
        
        self.table1 = QTable(self.sourceTab,"table1")
        self.table1.setGeometry(QRect(-3,7,410,460))
        self.table1.setNumRows(self.lsm.getSources())
        self.table1.setNumCols(self.lsm.getSourceColumns())
        self.table1.horizontalHeader().setLabel(0,self.tr("Source Name"))
        self.table1.horizontalHeader().setLabel(1,self.tr("Template Tree Type"))
        self.table1.horizontalHeader().setLabel(2,self.tr("MeqParm Table"))
        row=0
        for sname in self.lsm.s_table.keys():
         source=self.lsm.s_table[sname]
         self.table1.setText(row,0,QString(source.name))
         if source.treeType==POINT_TYPE:
          self.table1.setText(row,1,QString("Point"))
         elif source.treeType==GAUSS_TYPE:
          self.table1.setText(row,1,QString("Gaussian"))
         else:
          self.table1.setText(row,1,QString("Undefined"))
         self.table1.setText(row,2,QString(source.tableName))
         row+=1
        self.table1.adjustColumn(0)
        self.table1.adjustColumn(1)
        self.table1.adjustColumn(2)
        self.table1.setReadOnly(1)
        self.table1.setSorting(0)

        self.tabWidget.insertTab(self.sourceTab,QString.fromLatin1(""))
        sourceLayout.addWidget(self.table1)


##################### Tab 2
        self.punitTab = QWidget(self.tabWidget,"punitTab")

        # layout for table2
        sourceLayout=QVBoxLayout(self.punitTab)

        self.table2 = QTable(self.punitTab,"table2")
        self.table2.setGeometry(QRect(-3,-3,401,471))
        self.table2.setNumRows(self.lsm.getPUnits())
        self.table2.setNumCols(self.lsm.getPUnitColumns())
        self.table2.setSorting(1)
        self.tabWidget.insertTab(self.punitTab,QString.fromLatin1(""))

        self.table2.horizontalHeader().setLabel(PCOL_NAME,self.tr("PUnit Name"))
        self.table2.horizontalHeader().setLabel(PCOL_TYPE,self.tr("Type"))
        self.table2.horizontalHeader().setLabel(PCOL_SLIST,self.tr("Source List"))
        self.table2.horizontalHeader().setLabel(PCOL_CAT,self.tr("Category"))
        self.table2.horizontalHeader().setLabel(PCOL_BRIGHT,self.tr("Brightness"))
        self.table2.horizontalHeader().setLabel(PCOL_FOV,self.tr("FOV Distance"))
        self.table2.horizontalHeader().setLabel(PCOL_I,self.tr("I"))
        self.table2.horizontalHeader().setLabel(PCOL_Q,self.tr("Q"))
        self.table2.horizontalHeader().setLabel(PCOL_U,self.tr("U"))
        self.table2.horizontalHeader().setLabel(PCOL_V,self.tr("V"))
        self.table2.horizontalHeader().setLabel(PCOL_RA,self.tr("RA"))
        self.table2.horizontalHeader().setLabel(PCOL_DEC,self.tr("Dec"))
        row=0
        # use a hash table to match row number to name
        self.table2_names={}
        for sname in self.lsm.p_table.keys():
         punit=self.lsm.p_table[sname]
         if punit._patch_name==None:
          self.table2.setText(row,PCOL_NAME,QString(punit.name))
          mytype=punit.getType()
          if mytype==POINT_TYPE:
           self.table2.setText(row,PCOL_TYPE,self.tr("Point"))
          elif mytype==GAUSS_TYPE:
           self.table2.setText(row,PCOL_TYPE,self.tr("Ext"))
          else:
           self.table2.setText(row,PCOL_TYPE,self.tr("Patch"))
          # do not print all the source names in case of a patch
          if mytype==POINT_TYPE or mytype==GAUSS_TYPE:
           self.table2.setText(row,PCOL_SLIST,QString(str(punit.getSources())))
          else: #patch
           srclist=punit.getSources()
           self.table2.setText(row,PCOL_SLIST,QString(str(srclist[0])+"...."))

          self.table2.setText(row,PCOL_CAT,QString(str(punit.getCat())))
          self.table2.setText(row,PCOL_BRIGHT,QString(str(punit.getBrightness())))
          self.table2.setText(row,PCOL_FOV,QString(str(punit.getFOVDist())))
         # self.table2.setText(row,6,QString(str(punit.sp.getI())))
         # self.table2.setText(row,7,QString(str(punit.sp.getQ())))
         # self.table2.setText(row,8,QString(str(punit.sp.getU())))
         # self.table2.setText(row,9,QString(str(punit.sp.getV())))
         # self.table2.setText(row,10,QString(str(punit.sp.getRA())))
         # self.table2.setText(row,11,QString(str(punit.sp.getDec())))

          if mytype==POINT_TYPE or mytype==GAUSS_TYPE:
           self.table2.setText(row,PCOL_I,QString("MeqTree"))
           # see if this source is polarized
           [qq,uu,vv]=extract_polarization_parms(punit.getSixpack(),self.lsm.getNodeScope())
           if qq!=0:
             self.table2.setText(row,PCOL_Q,QString("Polarized"))
           else:
             self.table2.setText(row,PCOL_Q,QString("MeqTree"))
           if uu!=0:
             self.table2.setText(row,PCOL_U,QString("Polarized"))
           else:
             self.table2.setText(row,PCOL_U,QString("MeqTree"))
           if vv!=0:
             self.table2.setText(row,PCOL_V,QString("Polarized"))
           else:
             self.table2.setText(row,PCOL_V,QString("MeqTree"))
           self.table2.setText(row,PCOL_RA,QString("MeqTree"))
           self.table2.setText(row,PCOL_DEC,QString("MeqTree"))
          else: # patch
           self.table2.setText(row,PCOL_I,QString("N/A"))
           self.table2.setText(row,PCOL_Q,QString("N/A"))
           self.table2.setText(row,PCOL_U,QString("N/A"))
           self.table2.setText(row,PCOL_V,QString("N/A"))
           self.table2.setText(row,PCOL_RA,QString("N/A"))
           self.table2.setText(row,PCOL_DEC,QString("N/A"))

          # rememeber the name
          self.table2_names[sname]=row
          row+=1
         else:
          #print "This belongs to a patch"
          pass

        for i in range(self.table2.numCols()):
         self.table2.adjustColumn(i)
        sourceLayout.addWidget(self.table2)
        self.table2.setSorting(0)
        self.table2.setReadOnly(1)
        ### now some signals
        self.connect( self.table2, SIGNAL("clicked( int, int, int,const QPoint&)"),self.putable_getcell)


####### Tab 3 ############################
        self.imageTab= QWidget(self.tabWidget,"imageTab")
        sourceLayout=QVBoxLayout(self.imageTab)
        self.canvas=QCanvas(500,500)

        self.xlabel=QLabel(self.imageTab,"xlabel")
        self.xlabel.setText("0.0000000")
        self.xlabel.setMinimumWidth(self.xlabel.sizeHint().width())
        self.ylabel=QLabel(self.imageTab,"ylabel")
        self.ylabel.setText("0.0000000")
        self.ylabel.setMinimumWidth(self.ylabel.sizeHint().width())
        self.zlabel=QLabel(self.imageTab,"zlabel")
        self.zlabel.setText("f=000000000 T=0000000000")
        self.zlabel.setMinimumWidth(self.zlabel.sizeHint().width())
        self.sliderLabel= QLabel(self.imageTab,"slider")
        self.sliderLabel.setText("0.000")

        self.cview=MyCanvasView(self.canvas,self.imageTab,"canvas",self.xlabel,\
         self.ylabel,self.zlabel,self.sliderLabel,self.lsm,self)
        sourceLayout.addWidget(self.cview)

        layout1 = QHBoxLayout(None,0,6,"layout1")
        layout1.addWidget(self.xlabel)
        layout1.addWidget(self.ylabel)
        layout1.addWidget(self.zlabel)

        spacer1 = QSpacerItem(11,11,QSizePolicy.Expanding,QSizePolicy.Minimum)
        layout1.addItem(spacer1)

        ## add two buttons for enlarge(shrinking sources
        layout3 = QVBoxLayout(None,0,0,"layout3")
        # buttons
        self.enlarge_button=QPushButton(self.imageTab,"+")
        self.shrink_button=QPushButton(self.imageTab,"-")
        self.enlarge_button.setFlat(0)
        self.enlarge_button.setText("+")
        self.shrink_button.setFlat(0)
        self.shrink_button.setText("-")
        layout3.addWidget(self.enlarge_button)
        layout3.addWidget(self.shrink_button)
        layout1.addLayout(layout3)
        self.connect( self.enlarge_button, SIGNAL("clicked()"), self.cview.enlarge_display )
        self.connect( self.shrink_button, SIGNAL("clicked()"), self.cview.shrink_display )


        # slider
        layout2 = QVBoxLayout(None,0,6,"layout2")
        layout2.addWidget(self.sliderLabel)

        self.sliderCut = QSlider(0,100,1,0, QSlider.Horizontal,self.imageTab,"sliderCut")
        self.connect( self.sliderCut, SIGNAL("valueChanged( int )"), 
           self.cview.display)

        layout2.addWidget(self.sliderCut)

        layout1.addLayout(layout2)
        sourceLayout.addLayout(layout1)

        self.tabWidget.insertTab(self.imageTab,QString.fromLatin1(""))
        # show the image by default
        self.tabWidget.setCurrentPage(2)
######## End of Tabs #####################

        MainLayout.addWidget(self.tabWidget)

        self.fileOpenAction = QAction(self,"fileOpenAction")
        self.fileOpenAction.setIconSet(QIconSet(self.image1))
        self.fileSaveAction = QAction(self,"fileSaveAction")
        self.fileSaveAction.setIconSet(QIconSet(self.image2))
        self.fileSaveAsAction = QAction(self,"fileSaveAsAction")
        self.filePrintAction = QAction(self,"filePrintAction")
        self.filePrintAction.setIconSet(QIconSet(self.image3))
        self.fileExportAction = QAction(self,"fileExportAction")
        self.helpAboutAction = QAction(self,"helpAboutAction")
        self.viewZoom_WindowAction = QAction(self,"viewZoom_WindowAction")
        self.viewZoom_MoreAction = QAction(self,"viewZoom_MoreAction")
        self.viewZoom_LessAction = QAction(self,"viewZoom_LessAction")
        self.viewZoom_AllAction = QAction(self,"viewZoom_AllAction")
        self.view_selectAction = QAction(self,"view_selectAction")
        self.viewZoom_OptionsAction = QAction(self,"viewZoom_OptionsAction")
        self.view_nextAction= QAction(self,"view_nextAction")
        self.view_patchAction= QAction(self,"view_patchAction")
        self.view_refreshAction = QAction(self,"view_refreshAction")
        self.edit_moveAction = QAction(self,"edit_moveAction")
        self.edit_transformAction = QAction(self,"edit_transformAction")



        self.MenuBar = QMenuBar(self,"MenuBar")


        self.fileMenu = QPopupMenu(self)
        self.fileOpenAction.addTo(self.fileMenu)
        self.fileSaveAction.addTo(self.fileMenu)
        self.fileSaveAsAction.addTo(self.fileMenu)
        self.fileMenu.insertSeparator()
        self.fileExportAction.addTo(self.fileMenu)
        self.filePrintAction.addTo(self.fileMenu)
        self.fileMenu.insertSeparator()
        self.fileMenu.insertItem('&Close',self,SLOT('close()'),Qt.CTRL + Qt.Key_C)
        self.fileMenu.insertItem('&Quit',qApp,SLOT('closeAllWindows()'),Qt.CTRL + Qt.Key_Q)

        self.MenuBar.insertItem(QString(""),self.fileMenu,1)

        self.viewMenu = QPopupMenu(self)
        self.viewZoom_WindowAction.addTo(self.viewMenu)
        self.viewZoom_AllAction.addTo(self.viewMenu)
        #self.viewZoom_CancelAction.addTo(self.viewMenu)
        self.viewMenu.insertSeparator()
        self.viewZoom_MoreAction.addTo(self.viewMenu)
        self.viewZoom_LessAction.addTo(self.viewMenu)
        self.viewMenu.insertSeparator()
        self.view_refreshAction.addTo(self.viewMenu)
        self.viewMenu.insertSeparator()
        self.view_nextAction.addTo(self.viewMenu)
        self.viewZoom_OptionsAction.addTo(self.viewMenu)
        self.MenuBar.insertItem(QString(""),self.viewMenu,2)


        self.editMenu = QPopupMenu(self)
        self.edit_moveAction.addTo(self.editMenu)
        self.edit_transformAction.addTo(self.editMenu)
        self.editMenu.insertSeparator()
        self.view_selectAction.addTo(self.editMenu)
        self.view_patchAction.addTo(self.editMenu)
        self.MenuBar.insertItem(QString(""),self.editMenu,3)

        self.MenuBar.insertSeparator(4)

        self.helpMenu = QPopupMenu(self)
        self.helpAboutAction.addTo(self.helpMenu)
        self.MenuBar.insertItem(QString(""),self.helpMenu,5)


        self.languageChange()

        self.resize(QSize(530,670).expandedTo(self.minimumSizeHint()))
        self.clearWState(Qt.WState_Polished)

        self.connect(self.fileOpenAction,SIGNAL("activated()"),self.fileOpen)
        self.connect(self.fileSaveAction,SIGNAL("activated()"),self.fileSave)
        self.connect(self.fileSaveAsAction,SIGNAL("activated()"),self.fileSaveAs)
        self.connect(self.filePrintAction,SIGNAL("activated()"),self.filePrint)
        self.connect(self.fileExportAction,SIGNAL("activated()"),self.fileExport)

        self.connect(self.viewZoom_WindowAction,SIGNAL("activated()"),self.zoomStart)
        self.connect(self.viewZoom_MoreAction,SIGNAL("activated()"),self.zoomIn)
        self.connect(self.viewZoom_LessAction,SIGNAL("activated()"),self.zoomOut)
        self.connect(self.viewZoom_AllAction,SIGNAL("activated()"),self.zoomAll)
        self.connect(self.view_refreshAction,SIGNAL("activated()"),self.refreshView)
        self.connect(self.viewZoom_OptionsAction,SIGNAL("activated()"),self.changeOptions)
        self.connect(self.view_nextAction,SIGNAL("activated()"),self.viewNextMode)
        self.connect(self.edit_moveAction,SIGNAL("activated()"),self.moveItem)
        self.connect(self.edit_transformAction,SIGNAL("activated()"),self.linearTransform)
        self.connect(self.view_selectAction,SIGNAL("activated()"),self.viewSelectWindow)
        self.connect(self.view_patchAction,SIGNAL("activated()"),self.viewCreatePatches)

        self.connect(self.helpAboutAction,SIGNAL("activated()"),self.helpAbout)

    # event handelr fro PUnit table
    def putable_getcell(self,cellx,celly,button,point):
     #print cellx,celly
     #print button
     #print point
     if button==1:
       puname=self.table2.text(cellx,PCOL_NAME).ascii()
       plist=self.lsm.queryLSM(name=puname)
       punit=plist[0]
       sp=plist[0].getSP()
       win=None
       if punit.getType()==POINT_TYPE or punit.getType()==GAUSS_TYPE: 
        if (celly==PCOL_I):
         # get I node stub corresponding to this punit
         win=TreeDisp(self,puname,0,0,sp.stokesI())
        elif (celly==PCOL_Q):
         win=TreeDisp(self,puname,0,0,sp.stokesQ())
        elif (celly==PCOL_U):
         win=TreeDisp(self,puname,0,0,sp.stokesU())
        elif (celly==PCOL_V):
         win=TreeDisp(self,puname,0,0,sp.stokesV())
        elif (celly==PCOL_RA):
         win=TreeDisp(self,puname,0,0,sp.ra())
        elif (celly==PCOL_DEC):
         win=TreeDisp(self,puname,0,0,sp.dec())
        elif (celly==PCOL_SLIST):
         win=TreeDisp(self,puname,0,0,sp.sixpack())
       else: # this is a patch
        if (celly==PCOL_SLIST) or (celly==PCOL_I) or\
          (celly==PCOL_Q) or (celly==PCOL_U) or (celly==PCOL_V) or\
          (celly==PCOL_RA) or (celly==PCOL_DEC):
         win=TreeDisp(self,puname,0,0,sp.root())
       if win!=None:
         win.setTitle(puname)
         win.show()
     else: # sort column in descending order, with whole row 
       self.table2.sortColumn(celly,0,1)

    def exportPUTableTEX(self,filename="tab.tex"):
      mytable=self.table2
      self.exportTableTEX(mytable,filename)

    def exportSTableTEX(self,filename="tab.tex"):
      mytable=self.table1
      self.exportTableTEX(mytable,filename)

    def exportTableTEX(self,mytable,filename):
     # open file
     try:
      f=open(filename, 'w')
     except:
      print "Cannont open file %s",filename
      return
     nrows=mytable.numRows()
     ncols=mytable.numCols()
     # write preamble
     f.write("% Table Automatically generated by LSM browser\n")
     # header 
     hdr="||c||"+(ncols-1)*'c|'+'|'
     f.write("\\begin{tabular}{"+hdr+"}\\hline\n")
     thead=mytable.horizontalHeader()
     hdr=thead.label(0).ascii()
     for cj in range(ncols-1):
      hdr+="&"+thead.label(cj+1).ascii()
     f.write(hdr+"\\\\\\hline\\hline\n")
     # process each row
     for ci in range(nrows):
      hdr=mytable.text(ci,0).ascii()
      for cj in range(ncols-1):
       hdr+="&"+mytable.text(ci,cj+1).ascii()
      f.write(hdr+"\\\\\\hline\n")
     # close table
     f.write("\\end{tabular}\n")
     f.close()


    def languageChange(self):
        self.setCaption(self.__tr("File: "+self.lsm.getFileName()))
        self.tabWidget.changeTab(self.sourceTab,self.__tr("Source Table"))
        self.tabWidget.changeTab(self.punitTab,self.__tr("P-Unit Table"))
        self.tabWidget.changeTab(self.imageTab,self.__tr("Image"))
        self.fileOpenAction.setText(self.__tr("Open"))
        self.fileOpenAction.setMenuText(self.__tr("&Open..."))
        self.fileOpenAction.setAccel(self.__tr("Ctrl+O"))
        self.fileSaveAction.setText(self.__tr("Save"))
        self.fileSaveAction.setMenuText(self.__tr("&Save"))
        self.fileSaveAction.setAccel(self.__tr("Ctrl+S"))
        self.fileSaveAsAction.setText(self.__tr("Save As"))
        self.fileSaveAsAction.setMenuText(self.__tr("Save &As..."))
        self.fileSaveAsAction.setAccel(QString.null)
        self.filePrintAction.setText(self.__tr("Print"))
        self.filePrintAction.setMenuText(self.__tr("&Print..."))
        self.filePrintAction.setAccel(self.__tr("Ctrl+P"))
        self.fileExportAction.setText(self.__tr("Export"))
        self.fileExportAction.setMenuText(self.__tr("&Export..."))
        self.fileExportAction.setAccel(self.__tr("Ctrl+E"))

        self.helpAboutAction.setText(self.__tr("About"))
        self.helpAboutAction.setMenuText(self.__tr("&About"))
        self.helpAboutAction.setAccel(QString.null)
        self.viewZoom_WindowAction.setText(self.__tr("Enable Zoom"))
        self.viewZoom_WindowAction.setMenuText(self.__tr("Enable &Zoom"))
        self.viewZoom_WindowAction.setAccel(self.__tr("Ctrl+Z"))
        self.viewZoom_MoreAction.setText(self.__tr("Zoom In"))
        self.viewZoom_MoreAction.setMenuText(self.__tr("Zoom &In"))
        self.viewZoom_MoreAction.setAccel(self.__tr("Ctrl+I"))
        self.viewZoom_LessAction.setText(self.__tr("Zoom Out"))
        self.viewZoom_LessAction.setMenuText(self.__tr("Zoom O&ut"))
        self.viewZoom_LessAction.setAccel(self.__tr("Ctrl+U"))


        self.viewZoom_AllAction.setText(self.__tr("Zoom All"))
        self.viewZoom_AllAction.setMenuText(self.__tr("Zoom &All"))
        self.viewZoom_AllAction.setAccel(self.__tr("Ctrl+A"))
        self.view_refreshAction.setText(self.__tr("Refresh"))
        self.view_refreshAction.setMenuText(self.__tr("&Refresh"))
        self.view_refreshAction.setAccel(self.__tr("Ctrl+R"))
        self.viewZoom_OptionsAction.setText(self.__tr("Change Options for Plotting"))
        self.viewZoom_OptionsAction.setMenuText(self.__tr("Change &Options..."))
        self.viewZoom_OptionsAction.setAccel(self.__tr("Ctrl+O"))
        self.view_nextAction.setText(self.__tr("Next"))
        self.view_nextAction.setMenuText(self.__tr("Next &Mode..."))
        self.view_nextAction.setAccel(self.__tr("Ctrl+N"))

        self.edit_moveAction.setText(self.__tr("Move"))
        self.edit_moveAction.setMenuText(self.__tr("&Move"))
        self.edit_moveAction.setAccel(self.__tr("Ctrl+M"))
        self.edit_transformAction.setText(self.__tr("Transform"))
        self.edit_transformAction.setMenuText(self.__tr("&Transform..."))
        self.edit_transformAction.setAccel(self.__tr("Ctrl+T"))
        self.view_patchAction.setText(self.__tr("Create Patches"))
        self.view_patchAction.setMenuText(self.__tr("Create Patches (&Grid)..."))
        self.view_patchAction.setAccel(self.__tr("Ctrl+G"))
        self.view_selectAction.setText(self.__tr("Create Patches"))
        self.view_selectAction.setMenuText(self.__tr("Create Patches (&Window)"))
        self.view_selectAction.setAccel(self.__tr("Ctrl+W"))


        if self.MenuBar.findItem(1):
            self.MenuBar.findItem(1).setText(self.__tr("&File"))
        if self.MenuBar.findItem(2):
            self.MenuBar.findItem(2).setText(self.__tr("&View"))
        if self.MenuBar.findItem(3):
            self.MenuBar.findItem(3).setText(self.__tr("&Edit"))
        if self.MenuBar.findItem(5):
            self.MenuBar.findItem(5).setText(self.__tr("&Help"))




    def fileOpen(self):
        s=QFileDialog.getOpenFileName(".","*.lsm",self,"Open File","Choose LSM File to Load")
        #ry:
        self.lsm.load(s.ascii())
        self.setCaption("File: "+s.ascii())
        #except Exception:
        #print "file %s is not a valid LSM file" % s.ascii()

    def fileSave(self):
        if self.savefile==None:
         s=QFileDialog.getSaveFileName(".","*.lsm",self,"Save File","Choose File Name to Save LSM")
        else:
         s=QString(self.savefile)
        if s!=None and len(s)>0:
         if not s.endsWith(".lsm"):
           s.append(".lsm")
         self.savefile=s.ascii()
         self.lsm.save(self.savefile)
         self.setCaption("File: "+self.savefile)
         #print "Error saving to file %s" % self.savefile


    def fileSaveAs(self):
        s=QFileDialog.getSaveFileName(".","*.*",self,"Save File","Choose File Name to Save LSM")
        if s==None or len(s)==0:
         return
        self.savefile=s.ascii()
        try:
         self.lsm.save(self.savefile)
         self.setCaption("File: "+self.savefile)
        except Exception:
         print "Error saving to file %s" % self.savefile


    def filePrint(self):
        margin=10
        tab_margin=10
        if not self.cview.printer:
            self.cview.printer = QPrinter()
        if  self.cview.printer.setup(self.cview):
            pp=QPainter(self.cview.printer)
            # set font
            pp.setFont(self.cview.fonts['default'])
            # create special font as well
            ufont=QFont(self.cview.fonts['default'])
            ufont.setBold(1)
            ufont.setItalic(1)
            yPos=margin
            fm=pp.fontMetrics()
            metrics=QPaintDeviceMetrics(self.cview.printer)
            #self.canvas.drawArea(QRect(0,0,self.canvas.width(),self.canvas.height()),pp,False)
            # draw some text
            pWidth=metrics.width()
            pHeight=metrics.height()
            # get filename
            ttext=QString(self.lsm.getBaseFileName())
            # center
            pp.drawText((pWidth-ttext.length())/2,margin+yPos,metrics.width(),fm.lineSpacing(),Qt.ExpandTabs|Qt.DontClip,ttext)
            yPos=yPos+fm.lineSpacing()
            # create a Pixmap of the canvas
            pm=QPixmap(self.canvas.width(),self.canvas.height())
            pn=QPainter(pm)
            self.canvas.drawArea(self.canvas.rect(),pn)
            # now print the Pixmap
            pp.drawPixmap((pWidth-pm.width())/2,margin+yPos+1,pm)
            # destroy Pixmap painter
            pn.end()
 
            yPos=yPos+fm.lineSpacing()+self.canvas.height()
            ttext=QString("Additional Information :")
            # to draw this, use underlined font
            pp.setFont(ufont)
            ufm=pp.fontMetrics()
            pp.drawText(margin,margin+yPos,metrics.width(),ufm.lineSpacing(),Qt.ExpandTabs|Qt.DontClip,ttext)
            yPos=yPos+ufm.lineSpacing()
            pp.setFont(self.cview.fonts['default'])
            ttext=QString("Filename : "+self.lsm.getFileName())
            pp.drawText(margin,margin+yPos,metrics.width(),fm.lineSpacing(),Qt.ExpandTabs|Qt.DontClip,ttext)
            yPos=yPos+fm.lineSpacing()
            ttext=QString("PUnits : "+str(self.lsm.getPUnits())+" Sources : "+str(self.lsm.getSources()))
            pp.drawText(margin,margin+yPos,metrics.width(),fm.lineSpacing(),Qt.ExpandTabs|Qt.DontClip,ttext)
            yPos=yPos+fm.lineSpacing()
            tmpval0='%3.4f'%self.lsm.getMaxBrightness()
            tmpval1='%3.4f'%self.lsm.getMinBrightness()
            ttext=QString("Apparent Brightness : Max= "+tmpval0+" Min= "+ tmpval1+ "Jy")
            pp.drawText(margin,margin+yPos,metrics.width(),fm.lineSpacing(),Qt.ExpandTabs|Qt.DontClip,ttext)
            yPos=yPos+fm.lineSpacing()
            # print some info on the projection
            proj_info=self.cview.proj.info()
            if proj_info['state']==1:
             ttext=QString("Projection : phase centre = ("+str(proj_info['ra0'])+","+str(proj_info['dec0'])+"), rotation = "+str(proj_info['rot']))
            else:
             ttext=QString("Projection : None")
            pp.drawText(margin,margin+yPos,metrics.width(),fm.lineSpacing(),Qt.ExpandTabs|Qt.DontClip,ttext)
            yPos=yPos+fm.lineSpacing()
            n_punits=5
            if self.lsm.getPUnits()<n_punits:
              n_punits=self.lsm.getPUnits()
            ttext=QString("First "+str(n_punits)+" PUnits are:")
            pp.setFont(ufont)
            pp.drawText(margin,margin+yPos,metrics.width(),ufm.lineSpacing(),Qt.ExpandTabs|Qt.DontClip,ttext)
            yPos=yPos+ufm.lineSpacing()
            pp.setFont(self.cview.fonts['default'])
            # get the PUnits
            l_punits=self.lsm.queryLSM(count=n_punits)
            for punit in l_punits:
              tmpval='%3.4f'%punit.getBrightness()
              ttext=QString(punit.name+": brightness: "+tmpval+"Jy"+", type: ")
              if punit.getType()==POINT_TYPE:
                ttext=QString(ttext.ascii()+"Point")
              elif punit.getType()==GAUSS_TYPE:
                ttext=QString(ttext.ascii()+"Ext")
              else:
                ttext=QString(ttext.ascii()+"Patch")
              # polarization
              [qq,uu,vv]=extract_polarization_parms(punit.getSixpack(),self.lsm.getNodeScope())
              if qq!=0 or uu!=0 or vv!=0:
                tmpvalQ='%3.4e'%qq
                tmpvalU='%3.4e'%uu
                tmpvalV='%3.4e'%vv
                ttext=QString(ttext.ascii()+" % ["+tmpvalQ+", "+tmpvalU+", "+tmpvalV+"]")
              pp.drawText(margin+tab_margin,margin+yPos,metrics.width(),fm.lineSpacing(),Qt.ExpandTabs|Qt.DontClip,ttext)
              yPos=yPos+fm.lineSpacing()
            # print each PUnit on a line
            pp.end()
            if pp.isActive():
              pp.flush()

    def fileExport(self):
     win=ExportDialog(self,"Export Dialog",1,0,self)
     #win.setTitle("Testing...")
     win.show()

    # print to EPS
    def filePrintEPS(self,filename='./print.eps'):
      # write EPS
      if not self.cview.printer:
       self.cview.printer = QPrinter()
      self.cview.printer.setColorMode(QPrinter.Color)
      self.cview.printer.setOrientation(QPrinter.Portrait)
      self.cview.printer.setOutputToFile(True)
      self.cview.printer.setOutputFileName(filename)
      self.cview.printer.setPageSize(QPrinter.A4)
      if  self.cview.printer.setup(self.cview):
       pp=QPainter(self.cview.printer)
       self.canvas.drawArea(QRect(0,0,self.canvas.width(),self.canvas.height()),pp,False)
       pp.end()
       while pp.isActive():
        pp.flush()

       if self.cview.printer.outputToFile():
         self.filePStoEPS(self.cview.printer.outputFileName().ascii())
            
    # convert a PS file to an EPS file by changing the 
    # bounding box
    def filePStoEPS(self,filename):
     fn=open(filename,"r")   
     text=fn.read()
     fn.close()
     fileContent=QString(text)
     rx=QRegExp("%%BoundingBox:\\s*(-?[\\d\\.]+)\\s*(-?[\\d\\.]+)\\s*(-?[\\d\\.]+)\\s*(-?[\\d\\.]+)")
     pos = rx.search(fileContent)
     left = rx.cap(1).toFloat()
     top = rx.cap(4).toFloat()
     # parsed [left,bottom] [right,top]
     # replace  the margins
     # note padding is arbitrary
     padding=100 
     newstr="%%BoundingBox: "+str(left[0])+" "+str(top[0]-(self.canvas.height()-padding))+" "+str(left[0]+self.canvas.width()-padding)+" "+str(top[0])
     #print "Writing "+newstr
     # FIXME: need to check if canvas is scaled (zoomed)
     fileContent.replace(pos,rx.cap(0).length(),QString(newstr))
     #print "File modified",fileContent.ascii()

     newfilename=filename
     try:
       f = open(newfilename,'w+')
     except:
      print "Print EPS: file not writable"
      return
     f.write(fileContent.ascii())
     f.close()

    def fileExportImage(self,filename='./image',filetype="PNG"):
      pm=QPixmap(self.canvas.width(),self.canvas.height())
      pn=QPainter(pm)
      self.canvas.drawArea(self.canvas.rect(),pn)
      pm.save(filename,filetype)
      pn.end()

    def zoomStart(self):
       if (self.cview.zoom_status!=GUI_ZOOM_WINDOW):
        self.cview.zoom_status=GUI_ZOOM_WINDOW
        self.viewZoom_WindowAction.setMenuText(self.__tr("Disable &Zoom"))
        self.viewZoom_WindowAction.setText(self.__tr("Disable Zoom"))
        self.viewZoom_WindowAction.setAccel(self.__tr("Ctrl+Z"))
       else:
        self.cview.zoom_status=GUI_ZOOM_NONE
        self.viewZoom_WindowAction.setMenuText(self.__tr("Enable &Zoom"))
        self.viewZoom_WindowAction.setText(self.__tr("Enable Zoom"))
        self.viewZoom_WindowAction.setAccel(self.__tr("Ctrl+Z"))

    def zoomIn(self):
       m = self.cview.worldMatrix()
       dx=m.dx()
       dy=m.dy()
       # first scale
       m.scale( 2.0, 2.0 )
       # then move the origin by the scaled amount
       dx1=m.dx()
       dy1=m.dy()
       m.translate(2*(dx1-dx),2*(dy1-dy))

       g=QWMatrix()
       g.translate(2*(dx-dx1),2*(dy-dy1))
       g.scale(0.5,0.5)
       if (self.cview.tmstack) == None:
         self.cview.tmstack=g
       else:
         # calculate the total
         self.cview.tmstack*=g

       self.cview.setWorldMatrix( m )

    def zoomOut(self):
       m = self.cview.worldMatrix()
       dx=m.dx()
       dy=m.dy()
       g=QWMatrix()
       g.translate(0.5*dx,0.5*dy)
       g.scale(2,2)
       print m
       if (self.cview.tmstack) == None:
         self.cview.tmstack=g
       else:
         # calculate the total
         self.cview.tmstack*=g
       # first scale
       m.scale( 0.5, 0.5 )
       # then move the origin by the scaled amount
       m.translate(-dx*0.5,-dy*0.5)
       self.cview.setWorldMatrix( m )


    def zoomAll(self):
       # no Zooming if image tab is not seen
       if self.tabWidget.currentPageIndex() == 2:
        if (self.cview.tmstack) != None:
         m = self.cview.tmstack
         m *=self.cview.worldMatrix()
         self.cview.setWorldMatrix( m )
         self.cview.tmstack=None

    def refreshView(self):
      self.cview.updateCanvas()
 
    def helpAbout(self):
      tmp_str="<font color=\"blue\">LSM Browser</font><br/>"
      tmp_str+="<p>For more information please visit Timba MeqWiki Page at<br/>"
      tmp_str+="<span style=\"font-style: italic;\">http://lofar9.astron.nl/meqwiki/</span><br>"
      tmp_str+="</p>"
      dialog=SDialog(self)
      dialog.setInfoText(tmp_str)
      dialog.setTitle("Help")
      dialog.show()

    def changeOptions(self):
     win=OptionsDialog(self)
     win.show()

    def viewNextMode(self):
      win=FTDialog(self,f0=self.lsm.cells.grid.freq[self.lsm.cells.segments.freq.start_index],\
       f1=self.lsm.cells.grid.freq[self.lsm.cells.segments.freq.end_index],\
       fticks=self.lsm.cells.segments.freq.end_index,\
        f=self.cview.default_freq_index,\
        t0=self.lsm.cells.grid.time[self.lsm.cells.segments.time.start_index],\
       t1=self.lsm.cells.grid.time[self.lsm.cells.segments.time.end_index],\
       tticks=self.lsm.cells.segments.time.end_index,\
        t=self.cview.default_time_index)
      win.show()

    def viewSelectWindow(self):
      self.cview.zoom_status=GUI_SELECT_WINDOW

    def moveItem(self):
      self.cview.zoom_status=GUI_MOVE_START

    def linearTransform(self):
      win=TransDialog(self,"Linear Transform",1,0,self.cview,self.lsm) 
      win.show()


    # remove rows from PUnit table (table2)
    # with given name
    def removePUnitRows(self,pname_list):
        # create a list of row numbers corresponding
        # to the given names
        row_list=[]
        for pname in pname_list:
         if self.table2_names.has_key(pname):
          row_list.append(self.table2_names[pname])
        row_list.sort()
 
        # remove these rows
        if len(row_list)>0:
         self.table2.removeRows(row_list)
         # recreate name map
         self.table2_names={} 
         for ci in range(self.table2.numRows()):
          self.table2_names[self.table2.text(ci,0).ascii()]=ci



    # remove one row from PUnit table (table2)
    # with given name
    def removePUnitRow(self,pname):
        if self.table2_names.has_key(pname):
         self.table2.removeRow(self.table2_names[pname])
         # recreate name map
         self.table2_names={} 
         for ci in range(self.table2.numRows()):
           self.table2_names[self.table2.text(ci,0).ascii()]=ci

    # insert row to PUnit table (table2)
    # with given name
    def insertPUnitRow(self,pname):
        if self.table2_names.has_key(pname):
          print "Key alreay exits"
        else:
         punit=self.lsm.p_table[pname]
         if punit._patch_name==None:
          row=self.table2.numRows()
          self.table2.insertRows(row,1)
          self.table2.setText(row,PCOL_NAME,QString(punit.name))
          mytype=punit.getType()
          if mytype==POINT_TYPE:
           self.table2.setText(row,PCOL_TYPE,self.tr("Point"))
          elif mytype==GAUSS_TYPE:
           self.table2.setText(row,PCOL_TYPE,self.tr("Ext"))
          else:
           self.table2.setText(row,PCOL_TYPE,self.tr("Patch"))
          # do not print all the source names in case of a patch
          if mytype==POINT_TYPE:
           self.table2.setText(row,PCOL_SLIST,QString(str(punit.getSources())))
          elif mytype==GAUSS_TYPE:
           self.table2.setText(row,PCOL_SLIST,QString(str(punit.getSources())))
          else: #patch
           srclist=punit.getSources()
           self.table2.setText(row,PCOL_SLIST,QString(str(srclist[0])+"...."))

          self.table2.setText(row,PCOL_CAT,QString(str(punit.getCat())))
          self.table2.setText(row,PCOL_BRIGHT,QString(str(punit.getBrightness())))
          self.table2.setText(row,PCOL_FOV,QString(str(punit.getFOVDist())))

          if mytype==POINT_TYPE or mytype==GAUSS_TYPE  :
           self.table2.setText(row,PCOL_I,QString("MeqTree"))
           self.table2.setText(row,PCOL_Q,QString("MeqTree"))
           self.table2.setText(row,PCOL_U,QString("MeqTree"))
           self.table2.setText(row,PCOL_V,QString("MeqTree"))
           self.table2.setText(row,PCOL_RA,QString("MeqTree"))
           self.table2.setText(row,PCOL_DEC,QString("MeqTree"))
          else: # Patch
           self.table2.setText(row,PCOL_I,QString("N/A"))
           self.table2.setText(row,PCOL_Q,QString("N/A"))
           self.table2.setText(row,PCOL_U,QString("N/A"))
           self.table2.setText(row,PCOL_V,QString("N/A"))
           self.table2.setText(row,PCOL_RA,QString("N/A"))
           self.table2.setText(row,PCOL_DEC,QString("N/A"))

          # rememeber the name
          self.table2_names[pname]=row

    def viewCreatePatches(self):
      dialog=PatchOptionsDialog(self,"Patch Options",1,0,self.cview,self.lsm)
      dialog.show()


    def __tr(self,s,c = None):
        return qApp.translate("LSMWindow",s,c)

##########################################################################

