#!/usr/bin/env python

######################################################################
############# Dialog to change user options on plotting the LSM
######################################################################


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
from Timba.Meq import meq
import sys,string,math

from SDialog import *
import common_utils

class OptionsDialog(QDialog):
    def __init__(self,parent = None,name = "Change Options",modal = 1,fl = 0):
        QDialog.__init__(self,parent,name,modal,fl)
        self.setName("Change Options")
        self.setSizeGripEnabled(1)

        LayoutWidget = QVBoxLayout(self,11,6,"OuterLayout")
        self.tabWidget = QTabWidget(self,"tabWidget")

        # remember what has changed
        self.axes_changed=0
        self.turn_gridon=-1
        self.turn_legendon=-1
        self.display_source_type=-1
        self.coord_radians=-1
        self.plot_z_type=-1 # 0:brightness, 1,2,3,4,=IQUV

        # karma annotations
        self.karma_font0='hershey12';

        self.patch_center=-1# 0: geometric, 1: centroid 
        self.patch_method=-1# 0,1

        self.proj_changed=0
        self.proj_id=-1
        self.proj_ra0=self.parentWidget().cview.ra0
        self.proj_dec0=self.parentWidget().cview.dec0
################# Tab 1
        self.axisTab=QWidget(self.tabWidget,"axisTab")

        axistabLayout=QVBoxLayout(self.axisTab,11,6,"axistabLayout")

        layoutH=QHBoxLayout(None,0,6,"layoutH")
        axistabLayout.addLayout(layoutH)
        ################
        self.markerBG=QButtonGroup(self.axisTab,"markerBG")
        self.markerBG.setColumnLayout(0,Qt.Vertical)
        self.markerBG.layout().setSpacing(6)
        self.markerBG.layout().setMargin(6)
        bgLayout=QVBoxLayout(self.markerBG.layout())
        bgLayout.setAlignment(Qt.AlignCenter)

        layoutbgV=QVBoxLayout(None,0,6,"layoutbgV")
        self.marker_rad=QRadioButton(self.markerBG,"radian")
        layoutbgV.addWidget(self.marker_rad)
        self.marker_deg=QRadioButton(self.markerBG,"degrees")
        layoutbgV.addWidget(self.marker_deg)
        bgLayout.addLayout(layoutbgV)

        layoutH.addWidget(self.markerBG)

        if self.parentWidget().cview.default_coords=='rad':
         self.marker_rad.setChecked(1)
        else:
         self.marker_deg.setChecked(1)
        self.connect(self.markerBG, SIGNAL("clicked(int)"), self.markerBGradioClick)

        ###############
        self.gridBG=QButtonGroup(self.axisTab,"gridBG")
        self.gridBG.setColumnLayout(0,Qt.Vertical)
        self.gridBG.layout().setSpacing(6)
        self.gridBG.layout().setMargin(6)
        bgLayout=QVBoxLayout(self.gridBG.layout())
        bgLayout.setAlignment(Qt.AlignCenter)

        layoutbgV=QVBoxLayout(None,0,6,"layoutbgV")
        self.grid_ON=QRadioButton(self.gridBG,"gridon")
        layoutbgV.addWidget(self.grid_ON)
        self.grid_OFF=QRadioButton(self.gridBG,"gridoff")
        layoutbgV.addWidget(self.grid_OFF)
        bgLayout.addLayout(layoutbgV)


        layoutH.addWidget(self.gridBG)
        if(self.parentWidget().cview.grid_on==1):
         self.grid_ON.setChecked(1)
        else:
         self.grid_OFF.setChecked(1)

        self.connect(self.gridBG, SIGNAL("clicked(int)"), self.gridBGradioClick)
        #################
        self.xticksBG=QGroupBox(self.axisTab,"xticksBG")
        self.xticksBG.setColumnLayout(0,Qt.Vertical)
        self.xticksBG.layout().setSpacing(6)
        self.xticksBG.layout().setMargin(11)
        xticksLayout= QHBoxLayout(self.xticksBG.layout())
        xticksLayout.setAlignment(Qt.AlignCenter)

        self.textLabel_xspace = QLabel(self.xticksBG,"textLabel_xspace")
        xticksLayout.addWidget(self.textLabel_xspace)

        self.lineEdit_xspace = QLineEdit(self.xticksBG,"lineEdit_xspace")
        tempstr='%9.7f'%((self.parentWidget().cview.x_max-self.parentWidget().cview.x_min)/self.parentWidget().cview.xdivs)
        self.lineEdit_xspace.setText(tempstr)
        self.lineEdit_xspace.setReadOnly(1)
        xticksLayout.addWidget(self.lineEdit_xspace)
        spacer = QSpacerItem(31,31,QSizePolicy.Expanding,QSizePolicy.Minimum)
        xticksLayout.addItem(spacer)
        self.connect( self.lineEdit_xspace, SIGNAL("textChanged(const QString &)"), self.lineEditXspace)

        self.textLabel_xticks = QLabel(self.xticksBG,"textLabel_xticks")
        xticksLayout.addWidget(self.textLabel_xticks)

        self.lineEdit_xticks = QLineEdit(self.xticksBG,"lineEdit_xticks")
        self.lineEdit_xticks.setText(str(self.parentWidget().cview.xdivs))
        xticksLayout.addWidget(self.lineEdit_xticks)
        self.connect( self.lineEdit_xticks, SIGNAL("textChanged(const QString &)"), self.lineEditXticks )

        axistabLayout.addWidget(self.xticksBG)
        ################################### 
        self.yticksBG=QGroupBox(self.axisTab,"yticksBG")
        self.yticksBG.setColumnLayout(0,Qt.Vertical)
        self.yticksBG.layout().setSpacing(6)
        self.yticksBG.layout().setMargin(11)
        yticksLayout= QHBoxLayout(self.yticksBG.layout())
        yticksLayout.setAlignment(Qt.AlignCenter)

        self.textLabel_yspace = QLabel(self.yticksBG,"textLabel_yspace")
        yticksLayout.addWidget(self.textLabel_yspace)

        self.lineEdit_yspace = QLineEdit(self.yticksBG,"lineEdit_yspace")
        tempstr='%9.7f'%((self.parentWidget().cview.y_max-self.parentWidget().cview.y_min)/self.parentWidget().cview.ydivs)
        self.lineEdit_yspace.setText(tempstr)
        self.lineEdit_yspace.setReadOnly(1)
        yticksLayout.addWidget(self.lineEdit_yspace)
        self.connect( self.lineEdit_yspace, SIGNAL("textChanged(const QString &)"), self.lineEditYspace)
        spacer = QSpacerItem(31,31,QSizePolicy.Expanding,QSizePolicy.Minimum)
        yticksLayout.addItem(spacer)

        self.textLabel_yticks = QLabel(self.yticksBG,"textLabel_yticks")
        yticksLayout.addWidget(self.textLabel_yticks)

        self.lineEdit_yticks = QLineEdit(self.yticksBG,"lineEdit_yticks")
        self.lineEdit_yticks.setText(str(self.parentWidget().cview.ydivs))
        self.connect( self.lineEdit_yticks, SIGNAL("textChanged(const QString &)"), self.lineEditYticks)
        yticksLayout.addWidget(self.lineEdit_yticks)

        axistabLayout.addWidget(self.yticksBG)
        ################################### 
        self.fontBG= QGroupBox(self.axisTab,"fontBG")
        self.fontBG.setColumnLayout(0,Qt.Vertical)
        self.fontBG.layout().setSpacing(6)
        self.fontBG.layout().setMargin(11)
        fontBGLayout = QVBoxLayout(self.fontBG.layout())
        fontBGLayout.setAlignment(Qt.AlignCenter)

        self.fontButton= QPushButton(self.fontBG,"fontToolButton")
        fontBGLayout.addWidget(self.fontButton)
 
        self.connect( self.fontButton, SIGNAL("clicked()"), self.parentWidget().cview.chooseFont)

        axistabLayout.addWidget(self.fontBG)
        ####################################
        self.tabWidget.insertTab(self.axisTab,QString.fromLatin1(""))

################### Tab 2
        self.displayTab=QWidget(self.tabWidget,"displayTab")
        displaytabLayout=QVBoxLayout(self.displayTab,11,6,"displaytabLayout")
  
        ######## Group 1
        self.sourceBG=QButtonGroup(self.displayTab,"sourceBG")
        self.sourceBG.setColumnLayout(0,Qt.Vertical)
        self.sourceBG.layout().setSpacing(6)
        self.sourceBG.layout().setMargin(6)
        bgLayout=QVBoxLayout(self.sourceBG.layout())
        bgLayout.setAlignment(Qt.AlignCenter)

        layoutbgV=QVBoxLayout(None,0,6,"layoutbgV")
        self.source_point=QRadioButton(self.sourceBG,"points")
        layoutbgV.addWidget(self.source_point)
        self.source_cross=QRadioButton(self.sourceBG,"crosses")
        layoutbgV.addWidget(self.source_cross)
        self.source_ext=QRadioButton(self.sourceBG,"extended")
        layoutbgV.addWidget(self.source_ext)
        bgLayout.addLayout(layoutbgV)
        displaytabLayout.addWidget(self.sourceBG)

        if self.parentWidget().cview.display_point_sources=='cross':
         self.source_cross.setChecked(1)
        elif self.parentWidget().cview.display_point_sources=='point':
         self.source_point.setChecked(1)
        else:
         self.source_ext.setChecked(1)

        self.connect(self.sourceBG, SIGNAL("clicked(int)"), self.sourceBGradioClick)
        ######## Group 2
        self.legendBG=QButtonGroup(self.displayTab,"legendBG")
        self.legendBG.setColumnLayout(0,Qt.Vertical)
        self.legendBG.layout().setSpacing(6)
        self.legendBG.layout().setMargin(6)
        bgLayout=QVBoxLayout(self.legendBG.layout())
        bgLayout.setAlignment(Qt.AlignCenter)

        layoutbgV=QVBoxLayout(None,0,6,"layoutbgV")
        self.legend_ON=QRadioButton(self.legendBG,"legendON")
        layoutbgV.addWidget(self.legend_ON)
        self.legend_OFF=QRadioButton(self.legendBG,"legendOFF")
        layoutbgV.addWidget(self.legend_OFF)
        bgLayout.addLayout(layoutbgV)
        displaytabLayout.addWidget(self.legendBG)

       
        if self.parentWidget().cview.legend_on==0:
         self.legend_OFF.setChecked(1)
        else:
         self.legend_ON.setChecked(1)
        self.connect(self.legendBG, SIGNAL("clicked(int)"), self.legendBGradioClick)
        ###############################



        self.tabWidget.insertTab(self.displayTab,QString.fromLatin1(""))
################### Tab 3
        self.zTab=QWidget(self.tabWidget,"zTab")
        displaytabLayout=QVBoxLayout(self.zTab,11,6,"ztabLayout")
  
        ######## Group 1
        self.plotBG=QButtonGroup(self.zTab,"sourceBG")
        self.plotBG.setColumnLayout(0,Qt.Vertical)
        self.plotBG.layout().setSpacing(6)
        self.plotBG.layout().setMargin(6)
        bgLayout=QVBoxLayout(self.plotBG.layout())
        bgLayout.setAlignment(Qt.AlignCenter)

        layoutbgV=QVBoxLayout(None,0,6,"layoutbgV")
        self.z_A=QRadioButton(self.plotBG,"zA")
        layoutbgV.addWidget(self.z_A)
        self.z_I=QRadioButton(self.plotBG,"zI")
        layoutbgV.addWidget(self.z_I)
        self.z_Q=QRadioButton(self.plotBG,"zQ")
        layoutbgV.addWidget(self.z_Q)
        self.z_U=QRadioButton(self.plotBG,"zU")
        layoutbgV.addWidget(self.z_U)
        self.z_V=QRadioButton(self.plotBG,"zV")
        layoutbgV.addWidget(self.z_V)
        bgLayout.addLayout(layoutbgV)
        displaytabLayout.addWidget(self.plotBG)

        if self.parentWidget().cview.default_mode=='A':
         self.z_A.setChecked(1)
        elif self.parentWidget().cview.default_mode=='I':
         self.z_I.setChecked(1)
        elif self.parentWidget().cview.default_mode=='Q':
         self.z_Q.setChecked(1)
        elif self.parentWidget().cview.default_mode=='U':
         self.z_U.setChecked(1)
        elif self.parentWidget().cview.default_mode=='V':
         self.z_V.setChecked(1)


        self.connect(self.plotBG, SIGNAL("clicked(int)"), self.plotBGradioClick)
        ######## Group 2
        self.scaleBG=QButtonGroup(self.zTab,"scaleBG")
        self.scaleBG.setColumnLayout(0,Qt.Vertical)
        self.scaleBG.layout().setSpacing(6)
        self.scaleBG.layout().setMargin(6)
        bgLayout=QVBoxLayout(self.scaleBG.layout())
        bgLayout.setAlignment(Qt.AlignCenter)

        layoutbgV=QVBoxLayout(None,0,6,"layoutbgV")
        self.scale_lin=QRadioButton(self.scaleBG,"linear")
        layoutbgV.addWidget(self.scale_lin)
        self.scale_log=QRadioButton(self.scaleBG,"log")
        layoutbgV.addWidget(self.scale_log)
        bgLayout.addLayout(layoutbgV)
        displaytabLayout.addWidget(self.scaleBG)

       
        self.scale_lin.setChecked(1)
        self.connect(self.scaleBG, SIGNAL("clicked(int)"), self.scaleBGradioClick)
        ###############################



        self.tabWidget.insertTab(self.zTab,QString.fromLatin1(""))
################### Tab 4 - Karma Annotations
        self.karmaTab=QWidget(self.tabWidget,"karmaTab")
        karmaTabLayout=QVBoxLayout(self.karmaTab,11,6,"karmaTabLayout")


        ################# Box 1 - Font
        self.karmaFontBG=QGroupBox(self.karmaTab,"karmaFontBG")
        self.karmaFontBG.setColumnLayout(0,Qt.Vertical)
        self.karmaFontBG.layout().setSpacing(6)
        self.karmaFontBG.layout().setMargin(11)
        FLayout= QHBoxLayout(self.karmaFontBG.layout())
        FLayout.setAlignment(Qt.AlignCenter)

        self.textLabel_karmaf0= QLabel(self.karmaFontBG,"textLabel_karmaf0")
        FLayout.addWidget(self.textLabel_karmaf0)

        self.karmafontButton0= QPushButton(self.karmaFontBG,"fontToolButton")
        FLayout.addWidget(self.karmafontButton0)
 
        self.connect( self.karmafontButton0, SIGNAL("clicked()"), self.chooseFont_karma0 )


        karmaTabLayout.addWidget(self.karmaFontBG)
 
        ################# Box 2
        self.karmaColBG=QGroupBox(self.karmaTab,"karmaColBG")
        self.karmaColBG.setColumnLayout(0,Qt.Vertical)
        self.karmaColBG.layout().setSpacing(6)
        self.karmaColBG.layout().setMargin(11)
        TLayout= QHBoxLayout(self.karmaColBG.layout())
        TLayout.setAlignment(Qt.AlignCenter)

        self.textLabel_karmaC0= QLabel(self.karmaColBG,"textLabel_karmaC0")
        TLayout.addWidget(self.textLabel_karmaC0)

        karmaTabLayout.addWidget(self.karmaColBG)
 

        self.tabWidget.insertTab(self.karmaTab,QString.fromLatin1(""))
################### Tab 5
        self.patchTab=QWidget(self.tabWidget,"patchTab")
        patchtabLayout=QVBoxLayout(self.patchTab,11,6,"patchtabLayout")
  
        ######## Group 1
        self.patchBG=QButtonGroup(self.patchTab,"patchBG")
        self.patchBG.setColumnLayout(0,Qt.Vertical)
        self.patchBG.layout().setSpacing(6)
        self.patchBG.layout().setMargin(6)
        bgLayout=QVBoxLayout(self.patchBG.layout())
        bgLayout.setAlignment(Qt.AlignCenter)

        layoutbgV=QVBoxLayout(None,0,6,"layoutbgV")
        self.patch_geom=QRadioButton(self.patchBG,"geom")
        layoutbgV.addWidget(self.patch_geom)
        self.patch_centroid=QRadioButton(self.patchBG,"centroid")
        layoutbgV.addWidget(self.patch_centroid)
        bgLayout.addLayout(layoutbgV)
        patchtabLayout.addWidget(self.patchBG)

        if self.parentWidget().cview.lsm.default_patch_center=='G':
         self.patch_geom.setChecked(1)
        elif self.parentWidget().cview.lsm.default_patch_center=='C':
         self.patch_centroid.setChecked(1)

        ######## Group 2
        self.patchMBG=QButtonGroup(self.patchTab,"patchMBG")
        self.patchMBG.setColumnLayout(0,Qt.Vertical)
        self.patchMBG.layout().setSpacing(6)
        self.patchMBG.layout().setMargin(6)
        bgLayout=QVBoxLayout(self.patchMBG.layout())
        bgLayout.setAlignment(Qt.AlignCenter)

        layoutbgV=QVBoxLayout(None,0,6,"layoutbgV")
        self.patch_method1=QRadioButton(self.patchMBG,"method1")
        layoutbgV.addWidget(self.patch_method1)
        self.patch_method2=QRadioButton(self.patchMBG,"method2")
        layoutbgV.addWidget(self.patch_method2)
        bgLayout.addLayout(layoutbgV)
        patchtabLayout.addWidget(self.patchMBG)

        if self.parentWidget().cview.lsm.default_patch_method==1:
         self.patch_method1.setChecked(1)
        elif self.parentWidget().cview.lsm.default_patch_method==2:
         self.patch_method2.setChecked(1)


        self.connect(self.patchMBG, SIGNAL("clicked(int)"), self.patchMBGradioClick)
        ##########################################

        self.tabWidget.insertTab(self.patchTab,QString.fromLatin1(""))
################### Tab 6
        self.projTab=QWidget(self.tabWidget,"projTab")
        projtabLayout=QVBoxLayout(self.projTab,11,6,"projtabLayout")
  
        ######## Group 1
        self.projBG=QButtonGroup(self.projTab,"projBG")
        self.projBG.setColumnLayout(0,Qt.Vertical)
        self.projBG.layout().setSpacing(6)
        self.projBG.layout().setMargin(6)
        bgLayout=QVBoxLayout(self.projBG.layout())
        bgLayout.setAlignment(Qt.AlignCenter)

        layoutbgV=QVBoxLayout(None,0,6,"layoutbgV")
        self.proj_none=QRadioButton(self.projBG,"proj_none")
        layoutbgV.addWidget(self.proj_none)
        self.proj_sin=QRadioButton(self.projBG,"proj_sin")
        layoutbgV.addWidget(self.proj_sin)
        bgLayout.addLayout(layoutbgV)
        projtabLayout.addWidget(self.projBG)

        ######################
        self.projRABG=QGroupBox(self.projTab,"projRABG")
        self.projRABG.setColumnLayout(0,Qt.Vertical)
        self.projRABG.layout().setSpacing(6)
        self.projRABG.layout().setMargin(11)
        projRALayout= QHBoxLayout(self.projRABG.layout())
        self.textLabel_projRA = QLabel(self.projRABG,"textLabel_projRA")
        projRALayout.setAlignment(Qt.AlignCenter)
        projRALayout.addWidget(self.textLabel_projRA)

        self.lineEdit_projRA = QLineEdit(self.projRABG,"lineEdit_projRA")
        ll=common_utils.radToRA(self.parentWidget().cview.ra0)
        tempstr='%d:%d:%f'%(ll[0],ll[1],ll[2])
        self.lineEdit_projRA.setText(tempstr)
        self.lineEdit_projRA.setReadOnly(0)
        projRALayout.addWidget(self.lineEdit_projRA)
        self.connect( self.lineEdit_projRA, SIGNAL("textChanged(const QString &)"), self.lineEditprojRA)

        projtabLayout.addWidget(self.projRABG)
        ######################
        self.projDecBG=QGroupBox(self.projTab,"projDecBG")
        self.projDecBG.setColumnLayout(0,Qt.Vertical)
        self.projDecBG.layout().setSpacing(6)
        self.projDecBG.layout().setMargin(11)
        projDecLayout= QHBoxLayout(self.projDecBG.layout())
        self.textLabel_projDec = QLabel(self.projDecBG,"textLabel_projDec")
        projDecLayout.setAlignment(Qt.AlignCenter)
        projDecLayout.addWidget(self.textLabel_projDec)

        self.lineEdit_projDec = QLineEdit(self.projDecBG,"lineEdit_projDec")
        ll=common_utils.radToDec(self.parentWidget().cview.dec0)
        tempstr='%d:%d:%f'%(ll[0],ll[1],ll[2])
        self.lineEdit_projDec.setText(tempstr)
        self.lineEdit_projDec.setReadOnly(0)
        projDecLayout.addWidget(self.lineEdit_projDec)
        self.connect( self.lineEdit_projDec, SIGNAL("textChanged(const QString &)"), self.lineEditprojDec)

        projtabLayout.addWidget(self.projDecBG)
 
        if self.parentWidget().cview.proj.isOn():
         self.proj_sin.setChecked(1)
        else:
         self.proj_none.setChecked(1)

        self.connect(self.projBG, SIGNAL("clicked(int)"), self.projBGradioClick)
        ###############################



        self.tabWidget.insertTab(self.projTab,QString.fromLatin1(""))

############## end of Tabs
        LayoutWidget.addWidget(self.tabWidget)

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
######### End of button Bar
        LayoutWidget.addLayout(buttonLayout)

        self.languageChange()

        self.resize(QSize(528,381).expandedTo(self.minimumSizeHint()))
        self.clearWState(Qt.WState_Polished)

        self.connect(self.buttonHelp,SIGNAL("clicked()"),self.help)
        self.connect(self.buttonOk,SIGNAL("clicked()"),self.accept)
        self.connect(self.buttonCancel,SIGNAL("clicked()"),self.reject)


    def languageChange(self):
        self.setCaption(self.__tr("Change Options"))
        self.tabWidget.changeTab(self.axisTab,self.__tr("Axes"))
        self.markerBG.setTitle(self.__tr("Markers"))
        self.marker_rad.setText(self.__tr("Radians"))
        self.marker_deg.setText(self.__tr("Degrees"))
        self.gridBG.setTitle(self.__tr("Grid"))
        self.grid_ON.setText(self.__tr("On"))
        self.grid_OFF.setText(self.__tr("Off"))
        self.xticksBG.setTitle(self.__tr("X Ticks"))
        self.textLabel_xspace.setText(self.__tr("Spacing"))
        self.textLabel_xticks.setText(self.__tr("Count"))
        self.yticksBG.setTitle(self.__tr("Y Ticks"))
        self.textLabel_yspace.setText(self.__tr("Spacing"))
        self.textLabel_yticks.setText(self.__tr("Count"))

        QToolTip.add(self.markerBG,self.__tr("Change Coordinate display","Coords"))
        QWhatsThis.add(self.markerBG,self.__tr("Changes the Display of coordinats"))

        self.fontBG.setTitle(self.__tr("Font"))
        self.fontButton.setText(self.__tr("Choose..."))

        self.tabWidget.changeTab(self.displayTab,self.__tr("Display"))
        self.sourceBG.setTitle(self.__tr("Point Sources"))
        self.source_point.setText(self.__tr("Points"))
        self.source_cross.setText(self.__tr("Crosses"))
        self.source_ext.setText(self.__tr("Proportional Crosses"))
        self.legendBG.setTitle(self.__tr("Legend"))
        self.legend_ON.setText(self.__tr("On"))
        self.legend_OFF.setText(self.__tr("Off"))

        self.tabWidget.changeTab(self.zTab,self.__tr("Z axis"))
        self.plotBG.setTitle(self.__tr("Plot"))
        self.z_A.setText(self.__tr("Apparent Brightness"))
        self.z_I.setText(self.__tr("Stokes I"))
        self.z_Q.setText(self.__tr("Stokes Q"))
        self.z_U.setText(self.__tr("Stokes U"))
        self.z_V.setText(self.__tr("Stokes V"))
        self.scaleBG.setTitle(self.__tr("Scale"))
        self.scale_lin.setText(self.__tr("Linear"))
        self.scale_log.setText(self.__tr("Log"))

        self.tabWidget.changeTab(self.karmaTab,self.__tr("Karma Ann"))
        self.karmaFontBG.setTitle(self.__tr("Fonts"))
        self.textLabel_karmaf0.setText(self.__tr("Names: "))
        self.karmafontButton0.setText(self.__tr("Choose..."))
        self.karmaColBG.setTitle(self.__tr("Colors"))
        self.textLabel_karmaC0.setText(self.__tr("Names: "))


        self.tabWidget.changeTab(self.patchTab,self.__tr("Patches"))
        self.patchBG.setTitle(self.__tr("Phase Center"))
        self.patch_geom.setText(self.__tr("Geometric"))
        self.patch_centroid.setText(self.__tr("Weighted"))
        self.patchMBG.setTitle(self.__tr("Method"))
        self.patch_method1.setText(self.__tr("One"))
        self.patch_method2.setText(self.__tr("Two"))

        self.tabWidget.changeTab(self.projTab,self.__tr("Projection"))
        self.projBG.setTitle(self.__tr("Projection Method"))
        self.proj_sin.setText(self.__tr("SIN"))
        self.proj_none.setText(self.__tr("None"))
        self.projRABG.setTitle(self.__tr("Projection centre RA (hh:mm:ss)"))
        self.projDecBG.setTitle(self.__tr("Projection centre Dec (dd:mm:ss)"))

        self.buttonHelp.setText(self.__tr("&Help"))
        self.buttonHelp.setAccel(self.__tr("F1"))
        self.buttonOk.setText(self.__tr("&OK"))
        self.buttonOk.setAccel(QString.null)
        self.buttonCancel.setText(self.__tr("&Cancel"))
        self.buttonCancel.setAccel(QString.null)


    def __tr(self,s,c = None):
        return qApp.translate("OptionsDialog",s,c)

    ####################### Callbacks
    def markerBGradioClick(self,id):
     #print "marker BG button %d clicked" %id
     self.coord_radians=id


    def gridBGradioClick(self,id):
     #print "grid BG button %d clicked" %id
     if id==0:
      self.turn_gridon=1
     else:
      self.turn_gridon=0

    def lineEditXspace(self,newText):
     if newText.length()>0:
      try:
       fval=float(newText.ascii())
       #print "Line edit xspace text: %f" % fval
      except (TypeError,ValueError):
       print "invalid number "+unicode(newText)

    def lineEditXticks(self,newText):
     #print "Line edit xticks text: " + unicode(newText)
     if newText.length()>0:
      try:
       intval=int(newText.ascii())
       #print "intval =",intval
       if intval>0:
        self.lineEdit_xticks.setText(str(intval))
        tempstr='%9.7f'%((self.parentWidget().cview.x_max-self.parentWidget().cview.x_min)/intval)
        self.lineEdit_xspace.setText(tempstr)
        self.axes_changed=1
       else:
        raise TypeError
      except (TypeError,ValueError):
       #print "invalid number "+unicode(newText)
       self.lineEdit_xticks.setText(str(self.parentWidget().cview.xdivs))

    def lineEditYspace(self,newText):
     #print "Line edit yspace text: " + unicode(newText)
     pass
#     self.lineEdit_yticks.setText(newText)

    def lineEditYticks(self,newText):
     #print "Line edit yticks text: " + unicode(newText)
     if newText.length()>0:
      try:
       intval=int(newText.ascii())
       #print "intval =",intval
       if intval>0:
        self.lineEdit_yticks.setText(str(intval))
        tempstr='%9.7f'%((self.parentWidget().cview.y_max-self.parentWidget().cview.y_min)/intval)
        self.lineEdit_yspace.setText(tempstr)
        self.axes_changed=1
       else:
        raise TypeError
      except (TypeError,ValueError):
       #print "invalid number "+unicode(newText)
       self.lineEdit_yticks.setText(str(self.parentWidget().cview.ydivs))



    def sourceBGradioClick(self,id):
     self.display_source_type=id
     #print "source BG button %d clicked" %id

    def legendBGradioClick(self,id):
     #print "legend BG button %d clicked" %id
     if id==0:
      self.turn_legendon=1
     else:
      self.turn_legendon=0


    def plotBGradioClick(self,id):
     self.plot_z_type=id
     #print "plot BG button %d clicked" %id

    def scaleBGradioClick(self,id):
     #print "scale BG button %d clicked" %id
     pass


    def patchBGradioClick(self,id):
     self.patch_center=id
     #print "plot BG button %d clicked" %id

    def patchMBGradioClick(self,id):
     self.patch_method=id

    def projBGradioClick(self,id):
     self.proj_changed=1
     self.proj_id=id
     #print "plot BG button %d clicked" %id

    #select new font
    def chooseFont_karma0( self ) :
     newfont, ok = QFontDialog.getFont(self)
     if ok:
      namelist=QFontInfo(newfont)
      print namelist.family()
      print dir(namelist)
      #self.karma_font0=namelist[0]+namelist[1];
      print self.karma_font0

    ## parse RA hh:min:sec
    def lineEditprojRA(self,newText):
     if newText.length()>0:
      try:
       sl=string.split(newText.ascii(),':')
       fval=float(sl[0])+(float(sl[2])/60.0+float(sl[1]))/60.0
       fval*=math.pi/12
       self.proj_ra0=fval
       #print "Line edit xspace text: %f" % fval
       self.proj_changed=1
      except (TypeError,ValueError):
       pass
       #print "invalid number "+unicode(newText)

    def lineEditprojDec(self,newText):
     if newText.length()>0:
      try:
       sl=string.split(newText.ascii(),':')
       fval=float(sl[0])+(float(sl[2])/60.0+float(sl[1]))/60.0
       fval*=math.pi/180.0
       self.proj_dec0=fval
       #print "Line edit xspace text: %f" % fval
       self.proj_changed=1
      except (TypeError,ValueError):
       pass
       #print "invalid number "+unicode(newText)



    def accept(self):
     if self.axes_changed==1:
      xticks=int(self.lineEdit_xticks.text().ascii())
      yticks=int(self.lineEdit_yticks.text().ascii())
      self.parentWidget().cview.updateAxes(xticks,yticks)

     if self.turn_gridon!= -1:
      #print "Grid ON %d"%self.turn_gridon
      self.parentWidget().cview.grid_on=self.turn_gridon
      if self.turn_gridon==1:
       self.parentWidget().cview.axes.gridOn()
      else:
       self.parentWidget().cview.axes.gridOff()

     if self.display_source_type != -1:
      #print "Source type change %d"%self.display_source_type
      if self.display_source_type==0:
        self.parentWidget().cview.showPointSources(0)
      elif self.display_source_type==1:
        self.parentWidget().cview.showPointSources(1)
      else:
        self.parentWidget().cview.showPointSources(2)

     if self.coord_radians != -1:
      #print "Coordinate change %d"%self.coord_radians
      if self.coord_radians==0:
       self.parentWidget().cview.default_coords='rad'
       self.parentWidget().cview.axes.switchCoords('rad')
      else:
       self.parentWidget().cview.default_coords='deg'
       self.parentWidget().cview.axes.switchCoords('deg')


     if self.plot_z_type!=-1:
      #0,1,2,3,4: App. brightness, I,Q,U,V
      if self.plot_z_type==0:
       self.parentWidget().cview.updateDisplay('A')
      elif self.plot_z_type==1: 
       self.parentWidget().cview.updateDisplay('I')
      elif self.plot_z_type==2: 
       self.parentWidget().cview.updateDisplay('Q')
      elif self.plot_z_type==3: 
       self.parentWidget().cview.updateDisplay('U')
      elif self.plot_z_type==4: 
       self.parentWidget().cview.updateDisplay('V')

     if  self.turn_legendon!=-1:
      if  self.turn_legendon==1:
       self.parentWidget().cview.showLegend(1)
      else:
       self.parentWidget().cview.showLegend(0)

     if self.patch_center!=-1:
      if self.patch_center==0:
        self.parentWidget().cview.lsm.default_patch_center='G'
      else: 
        self.parentWidget().cview.lsm.default_patch_center='C'

     if self.patch_method!=-1:
        self.parentWidget().cview.lsm.default_patch_method=self.patch_method+1

     
     if self.proj_changed!=0:
       if self.proj_id==0:
        self.parentWidget().cview.proj_on=0
       else:
        self.parentWidget().cview.proj_on=1
       
       self.parentWidget().cview.ra0=self.proj_ra0
       self.parentWidget().cview.dec0=self.proj_dec0
       self.parentWidget().cview.updateCanvas()

     self.parentWidget().canvas.update()

     QDialog.accept(self)

    def reject(self):
     QDialog.reject(self)


    def help(self):
      tmp_str="""<div style="text-align: left;"><span style="font-weight: bold;">Options
               Are:</span><br>
               <ol>
               <li><u>Axes</u><br>
               <ul>
               <li>Change grid markings (Radians or Degrees)</li>
               <li>Change grid spacings: X and Y axes</li>
               <li>Show/Hide grid (when creating patches, the grid must be ON)</li>
               <li>Change the font used for axes</li>
               </li>
               </ul>
               </li>
               <li><u>Display</u><br>
               <ul>
               <li>Select how to plot sources (points, crosses etc.)</li>
               <li>Colour bar or a Legend?</li>
               </ul>
               </li>
               <li><u>Z axis</u><br>
               <ul>
               <li>What to display (Brightness, Stokes IQUV)
               </li>
               <li>Scale to use in plotting and colours</li>
               </ul>
               </li>
               <li><u>Cells</u><br>
               <ul>
               <li>Change the Frequency and Time grid on which to visualise the LSM</li>
               </ul>
               </li>
               </ol>
               </div>
               """
      dialog=SDialog(self)
      dialog.setInfoText(tmp_str)
      dialog.setTitle("Help")
      dialog.show()

def main(args):
  app=QApplication(args)
  win=OptionsDialog(None)
  win.show()
  app.connect(app,SIGNAL("lastWindowClosed()"),
               app,SLOT("quit()"))
  app.exec_loop()

if __name__=="__main__":
   main(sys.argv)
