#!/usr/bin/python

##########################################################
# Everything done with the canvas of the image
#
##########################################################


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
from qt import *
from qtcanvas import *

from LSM import *
import transform
from SDialog import *
from PatchDialog import *
#from PatchGLDialog import *

from common_utils import *
from canvas_inner import *

import Timba.array

##########################################################################



class MyCanvasView(QCanvasView):
  def __init__(self,canvas,parent,name,xlabel,ylabel,zlabel,sliderlabel,lsm_object,parent_window):
    QCanvasView.__init__(self,canvas,parent,name)
    ############ initialize attributes
    
    # number of divisions in axes
    self.xdivs=5
    self.ydivs=5
    self.grid_on=0
    self.legend_on=0
    self.axes_on=1
    self.title_on=0
    self.obsres_on=0 # beam
    self.obswin_on=0
    self.proj_on=0 # projection 0: off, 1:0n
    self.display_point_sources='pcross' #cross,point,pcross

    self.canvas().setDoubleBuffering(True) 
    self.viewport().setMouseTracking(1) 

    self.default_freq_index=0
    self.default_time_index=0
    # display A, I, Q, U, V
    # A : Apparent brightness
    self.default_mode='A'

    # display coordinates in Radians (rad) or degrees (deg)
    self.default_coords='deg' 

    # fonts 
    font=QFont(QApplication.font())
    font.setPointSize(10)
    font1=QFont(QApplication.font())
    font1.setPointSize(11)
    self.fonts={}
    self.fonts['default']=font
    self.fonts['title']=font1
    self.fonts['xlabel']=font1
    self.fonts['ylabel']=font1
    self.fonts['xcoord']=font
    self.fonts['ycoord']=font


    self.font=font

    self.lsm=lsm_object
    self.parent=parent_window

    (tmp_max,tmp_min,tmp_abs_min)=self.lsm.getBrightnessLims()
    self.max_brightness=tmp_max
    self.min_brightness=tmp_min
    self.abs_min_brightness=tmp_abs_min

    # sanity check
    if self.max_brightness==self.abs_min_brightness:
      self.max_brightness=self.abs_min_brightness+1.0
    if self.abs_min_brightness==0.0:
      self.abs_min_brightness=1e-6
    # get boundaries (using ObsWin  ?)
    # 1. create Cell using ObsWin, 2. send request to meqserver, 
    # 3. scan the result set to determine boundaries of
    #    RA, Dec, I, Q, U, V 
    bounds=self.lsm.getBounds()
    # now determine actual geometric bundaries
    # refined boundaries are
    #    ______________________
    #    h1 
    #    ----------------------
    #    
    #
    #    h2
    #  
    #    ----------------------
    #    h3  
    #    ----------------------
    #    h4 
    #    ----------------------
    #    ^ w1 ^ w2 ^    w3   ^ w4 ^ w5(legend)
    # this reduces to the following boundaries for transformation
    # so d1=w1+w2, d2=w4, d3=h1, d4= h3+h4
    # lets take w4 as padding.
    #    +-----------------+
    #    |         d3      |
    #    | d1 |---------|d2|
    #    .    .         .  .
    #    .    .         .  .
    #    |    +---------+  |
    #    +_____ d4_________+

    # setup 
    self.h1=self.h2=self.h3=self.h4=0
    self.w1=self.w2=self.w3=self.w4=self.w5=0
    self.pad=0 # padding : added to width w3
    # now calculate the real values
    self.calculateDefaultDims()

    self.d1=self.w1+self.w2
    self.d2=self.w4
    self.d3=self.h1
    self.d4=self.h3+self.h4

    # limits in real coordinates
    self.ra_min=bounds['min_RA']
    self.ra_max=bounds['max_RA']
    self.dec_min=bounds['min_Dec']
    self.dec_max=bounds['max_Dec']
    # add a margin of 2% of the real size
    xmarg=abs(self.ra_max-self.ra_min)*0.02
    ymarg=abs(self.dec_max-self.dec_min)*0.02
    if xmarg==0:
      xmarg=0.001
    if ymarg==0:
      ymarg=0.001
    self.ra_min=self.ra_min-xmarg
    self.ra_max=self.ra_max+xmarg
    self.dec_min=self.dec_min-ymarg
    self.dec_max=self.dec_max+ymarg

    ############# create phase center and a projector
    self.ra0=(self.ra_min+self.ra_max)*0.5
    self.dec0=(self.dec_min+self.dec_max)*0.5
    self.proj=transform.Projector(self.ra0,self.dec0,0)
    if self.proj_on==0: self.proj.Off() # turn off projection
    ### now recalculate the bounds according to the projection params
    left_down=self.proj.give_limits(self.ra_min,self.ra_max,self.dec_min,self.dec_max)
    self.x_min=left_down[0]
    self.x_max=left_down[1]
    self.y_min=left_down[2]
    self.y_max=left_down[3]


    # sanity check
    if self.x_min==self.x_max:
      self.x_max=self.x_min+0.1
    if self.y_min==self.y_max:
      self.y_max=self.y_min+0.1
    # canvas size
    H=self.canvas().height()
    W=self.canvas().width()
    # scales
    if (self.x_max != self.x_min):
     self.x_scale=(W-self.d1-self.d2)/(self.x_max-self.x_min)
    else: 
     self.x_scale=1.0

    if (self.y_max != self.y_min):
     self.y_scale=(H-self.d3-self.d4)/(self.y_max-self.y_min)
    else:
     self.y_scale=1.0

    # Tool Tip
    self.tooltip=DynamicTip(self)
    # printer
    self.printer=QPrinter()

    ############ create p-unit list
    # table for all p-units plots on canvas
    self.p_tab={}
    ## query only selected punits
    if self.lsm.display_punits>0:
      pulist=self.lsm.queryLSM(count=self.lsm.display_punits)
    else:
      pulist=self.lsm.queryLSM(all=1)

    # plot all p-units/sources
    for punit in pulist:
     sname=punit.name
     mytype=punit.getType()
     if mytype==POINT_TYPE:
      #xys=self.globalToLocal(punit.sp.getRA(),punit.sp.getDec())
      self.p_tab[sname]=PointSourceDisplay(sname,self)
     elif mytype==GAUSS_TYPE:
      self.p_tab[sname]=GaussianSourceDisplay(sname,self)
     else:
      # we have a patch
      retval=punit.getLimits()
      self.p_tab[sname]=PatchDisplay(sname,self,retval[0],retval[1],\
                 retval[2],retval[3])


    ############ create axes/grid
    if self.axes_on==1:
     self.axes=Axes(self,bounds,self.xdivs,self.ydivs)
     self.axes.gridOff()
    else:
     self.axes=None

    ########### create title
    if self.title_on==1:
     self.title=FontHorizImage("Title",self.canvas(),self.fonts['title'],self.pad)
     xys=self.globalToLocal(self.x_max, self.y_max)
     self.title.move(xys[0]+(self.w3-self.title.width())/2,xys[1]-self.h1)
     self.title.setZ(0)
     self.title.show()
    else:
     self.title=None

    ############ create legend
    self.legend=None

    ############ create zoom-window
    self.zwindow=ZWindow(self.canvas())
    #self.zwindow.hide()
    self.zoom_status=GUI_ZOOM_NONE
    # store transformation matrices zoom in/out
    self.tmstack=None



    ############ create slider
    # done by the parent

    ############ create x,y,z labels
    self.xlabel=xlabel
    self.ylabel=ylabel
    self.zlabel=zlabel
    self.resetFTindices()
    #head=self.lsm.getCurrentFreqTime(self.default_freq_index,self.default_time_index)

    #tmpval=stdForm(head[0],'%3.4f')
    #headstr="("+tmpval[0]+tmpval[1]+"Hz,"
    #tmpval=stdForm(head[1],'%3.4f')
    #headstr=headstr+tmpval[0]+tmpval[1]+"s)"
    #self.zlabel.setText("<font color=\"blue\">"+headstr+"</font>")


    self.slabel=sliderlabel

    # anything moving is kept here
    self.__moving=0
 

    
  def contentsMouseMoveEvent(self,e):
    point = self.inverseWorldMatrix().map(e.pos())
    xys=self.localToGlobal(point.x(),point.y())

    if self.default_coords=='rad':
     x_val='%2.6f'%xys[0]
     y_val='%2.6f'%xys[1]
     self.xlabel.setText("<font color=\"blue\">"+str(x_val)+"</font>")
     self.ylabel.setText("<font color=\"blue\">"+str(y_val)+"</font>")
    else: # degrees
     tmpval=radToRA(xys[0])
     self.xlabel.setText("<font color=\"blue\">"+str(tmpval[0])+"<sup>o</sup>"+str(tmpval[1])+"<sup>'</sup>"+str(tmpval[2])+"<sup>''</sup>""</font>")
     tmpval=radToDec(xys[1])
     self.ylabel.setText("<font color=\"blue\">"+str(tmpval[0])+"<sup>o</sup>"+str(tmpval[1])+"<sup>'</sup>"+str(tmpval[2])+"<sup>''</sup>""</font>")

    # move zoom window
    if self.zoom_status==GUI_ZOOM_START or\
       self.zoom_status==GUI_SELECT_START:
     self.zwindow.setLowerRight(point.x(),point.y())
     self.zwindow.show()
     self.canvas().update()


    if self.zoom_status==GUI_MOVE_START and self.__moving:
     point = self.inverseWorldMatrix().map(e.pos());
     self.__moving.cross.moveBy(point.x() - self.__moving_start.x(),point.y() - self.__moving_start.y())
     self.__moving.pcross.moveBy(point.x() - self.__moving_start.x(),point.y() - self.__moving_start.y())
     self.__moving.circle.moveBy(point.x() - self.__moving_start.x(),point.y() - self.__moving_start.y())
     self.__moving_start = point
     self.canvas().update()

    return
   
  def contentsMousePressEvent(self,e):
   point = self.inverseWorldMatrix().map(e.pos())
   button=e.button()
   #button==Qt.LeftButton:
   if self.zoom_status==GUI_ZOOM_NONE:
    ilist = self.canvas().collisions(point) #QCanvasItemList ilist
    head=self.lsm.getCurrentFreqTime(self.default_freq_index,self.default_time_index)
    tmpstr=stdForm(head[0],'%3.4f')
    headstr="<font color=\"blue\">("+tmpstr[0]+tmpstr[1]+"Hz, "
    tmpstr=stdForm(head[1],'%3.4f')
    headstr=headstr+tmpstr[0]+tmpstr[1]+"s)</font><br/>Sources:<br/>"
    tmp_str=headstr+"<ul>"
    found_anything=0
    #print ilist
    for each_item in ilist:
     if each_item.rtti()==POINT_SOURCE_RTTI:
      # set flag
      found_anything=1
      tmp_str=tmp_str+"<li>Label : "+each_item.name+" "
      # print brightness value as well
      punit=self.lsm.p_table[each_item.name]
      br=""
      if punit != None:
       if self.default_coords=='rad':
         br="<br/>At RA,DEC [%5.4f, %5.4f] (rad)<br/>App. Brightness: "%(punit.sp.getRA(),punit.sp.getDec())
         ab=stdForm(punit.getBrightness(self.default_mode,self.default_freq_index, self.default_time_index),"%5.3f")
         br=br+ab[0]+ab[1]+"Jy"
       else:
         tmpval=radToRA(punit.sp.getRA())
         br="<br/>At RA,DEC ["+str(tmpval[0])+"<sup>o</sup>"+str(tmpval[1])+"<sup>'</sup>"+str(tmpval[2])+"<sup>''</sup>,"
         tmpval=radToDec(punit.sp.getDec())
         br=br+str(tmpval[0])+"<sup>o</sup>"+str(tmpval[1])+"<sup>'</sup>"+str(tmpval[2])+"<sup>''</sup>] (deg)<br/>App. Brightness: "
         ab=stdForm(punit.getBrightness(self.default_mode,self.default_freq_index, self.default_time_index),"%5.3f")
         br=br+ab[0]+ab[1]+"Jy"
       # extract MeqParms if any
       [qq,uu,vv]=extract_polarization_parms(punit.getSixpack(),self.lsm.getNodeScope())
       [ra,dec,sI]=extract_parms(punit.getSixpack(),self.lsm.getNodeScope())
       br=br+"<br/>Polarization I="+str(sI)
       if (qq!=0 or uu!=0 or vv !=0):
         br=br+"<br/>Q="+str(qq)+"<br/>U="+str(uu)+"<br/>V="+str(vv)
      tmp_str+=br
      
      if punit._patch_name !=None:
       tmp_str+="<br/>patch <font color=\"blue\">"+punit._patch_name+"</font></li>"
      else:
       tmp_str+="</li>"
 

    if found_anything != 0: # not empty
     tmp_str=tmp_str+"</ul>"
     dialog=SDialog(self)
     dialog.setInfoText(tmp_str)
     dialog.show()
    else: # list is empty, no point sources, may be patches
     # re-scan the list for patches 
     found_anything=0
     tmp_str=""
     for each_item in ilist:
      if each_item.rtti()==PATCH_IMAGE_RTTI:
       #print "Found a Patch name %s"%each_item.parent.name
       #print each_item.image
       found_anything=1
       mimes=QMimeSourceFactory()
       patch_img=each_item.image.scale(100,100,QImage.ScaleMin)
       patch_img.invertPixels()
       mimes.setImage("img"+each_item.parent.name,patch_img)
       tmp_str+="<p><u>"+each_item.parent.name+"</u><br/>"
       tmp_str+="Image: <img source=\"img"+each_item.parent.name+"\" alt=\"image\""
       tmp_str+=" title=\"Image Title\" border=\"1\" style=\"width: 262px; height: 300px;\"/>"
       tmp_str+="<br/><br/></p>"
       #print tmp_str
     if found_anything !=0: 
      dialog=SDialog(self)
      dialog.textEdit.setMimeSourceFactory(mimes)
      dialog.setInfoText(tmp_str)
      dialog.setTitle("Patch Info")
      dialog.show()
     """try:
        # get vellsets if any
        punit=self.lsm.p_table[each_item.parent.name]
        lims=punit.sp.getValueSize(self.default_mode,\
         self.default_freq_index,\
         self.default_time_index)
        # the return type should be a Timba.array
        my_arr=punit.sp.getValue(self.default_mode,\
         self.default_freq_index,\
         self.default_time_index)
        # create GL window
        #wn=PatchGLDialog(self,"Patch",1,0,my_arr,lims)
        #wn.show()
       except:
        pass
     """

   if self.zoom_status==GUI_ZOOM_WINDOW:
    self.zoom_status=GUI_ZOOM_START
    self.zwindow.setUpperLeft(point.x(),point.y())
    self.zwindow.setLowerRight(point.x(),point.y())
    self.zwindow.show()
    self.canvas().update()
   # select window
   if self.zoom_status==GUI_SELECT_WINDOW:
    self.zoom_status=GUI_SELECT_START
    self.zwindow.setUpperLeft(point.x(),point.y())
    self.zwindow.setLowerRight(point.x(),point.y())
    self.zwindow.show()
    self.canvas().update()

   # moving things
   if self.zoom_status==GUI_MOVE_START:
    ilist = self.canvas().collisions(point) #QCanvasItemList ilist
    item=ilist.pop(0)
    if item.rtti()==POINT_SOURCE_RTTI:
     # enable moving
     print item.name
     # add all canvas items by that name to the move list
     # self.p_tab[sname]
     self.__moving=self.p_tab[item.name]
     self.__moving_start=point
   return

  def contentsMouseReleaseEvent(self,e):
   point = self.inverseWorldMatrix().map(e.pos())
   #cpoint=self.contentsToViewport(point)
   if self.zoom_status==GUI_ZOOM_START:
    self.zoom_status=GUI_ZOOM_WINDOW
    self.zwindow.hide()
    # zoom window
    mm=QWMatrix()
    # limit to stop accidental zooming
    win_delta=5
    if abs(point.x()-self.zwindow.left_x)>win_delta:
     xsc=float(self.visibleWidth())/(point.x()-self.zwindow.left_x)
    else:
     xsc=1
    if abs(point.y()-self.zwindow.left_y)>win_delta:
     ysc=float(self.visibleHeight())/(point.y()-self.zwindow.left_y)
    else:
     ysc=1
    if xsc<ysc:
      ysc=xsc
    else:
      xsc=ysc
    # if window is too small, do nothing
    if (xsc==1) and (ysc==1): return;
    # if we have negative scales also do nothing
    if (xsc<0) or (ysc<0): return;
    mm.scale(xsc,ysc)
    mm.translate(-self.zwindow.left_x,-self.zwindow.left_y)
    #mm*=self.worldMatrix()
    # store old matrix 
    #g = self.worldMatrix()
    g=QWMatrix()
    #print "inverse translation ",self.zwindow.left_x,self.zwindow.left_y
    g.translate(self.zwindow.left_x,self.zwindow.left_y)
    #print "inverse scaling ",1/xsc, 1/ysc
    g.scale( 1/xsc, 1/ysc)
    #g*=self.worldMatrix()

    self.tmstack=g
    self.setWorldMatrix(mm)
    self.canvas().update()

   if self.zoom_status==GUI_SELECT_START:
    self.zoom_status=GUI_ZOOM_NONE
    self.zwindow.hide()
    ilist = self.canvas().collisions(self.zwindow.getRect()) #QCanvasItemList ilist
    tmp_str=""
    # create a list of point source names,
    # if a patch is created
    psource_list=[]
    for each_item in ilist:
     if each_item.rtti()==POINT_SOURCE_RTTI:
      tmp_str+=" "+each_item.name+"<br>"
      psource_list.append(each_item.name)
    # if empty, do nothing
    if len(psource_list)==0:
     return
    dialog=PatchDialog(self)
    dialog.setInfoText(tmp_str)
    if dialog.exec_loop() == QDialog.Accepted:
     # create a patch
     retval=self.lsm.createPatch(psource_list)
     if retval != None:
      # successfully created patch
      #print "created patch %s"%retval[0]
      # remove these sources from PUnit table on main window
      self.parent.removePUnitRows(psource_list)
      # update the GUI
      self.p_tab[retval[0]]=PatchDisplay(retval[0],self,retval[1],retval[2],\
                 retval[3],retval[4])
      # update PUnit table on main window
      self.parent.insertPUnitRow(retval[0])
    self.canvas().update()


   if self.zoom_status==GUI_MOVE_START:
    self.zoom_status=GUI_ZOOM_NONE
    #print self.__moving.name
    xys=self.localToGlobal(point.x(),point.y())
    # update the LSM
    if self.__moving!=0:
     self.lsm.move_punit(self.__moving.name,xys[0],xys[1])
    # reset move
    self.__moving=0
    self.__moving_start=0
   return


  # return local coordinates for given global coordinates
  def localToGlobal(self,x,y):
   ll=[-(x-self.d1)/self.x_scale+self.x_max,-(y-self.d3)/self.y_scale+self.y_max]
   # do a projection from spherical to rectangular
   l=self.proj.rt_to_sp(ll[0],ll[1])
   return l

  # return global coordinates for given local coordinates
  def globalToLocal(self,x,y):
   # do a projection from rectangular to spherical
   l=self.proj.sp_to_rt(x,y)
   ll=[(-l[0]+self.x_max)*self.x_scale+self.d1,(self.y_max-l[1])*self.y_scale+self.d3]
   return ll

  # return a colour according to brightness
  def getColor(self,z):
    if self.max_brightness==0:
     return QColor(1,1,1) # no color : when Q,U,V is zero

    cl=(z-self.min_brightness)/(self.max_brightness-self.min_brightness) # normalized in [0,1] 
    if cl < 0.25:
     return QColor(0,int(cl*256*4),255)
    elif cl < 0.5:
     return QColor(0,255,int((2-cl*4)*256))
    elif cl < 0.75:
     return QColor(int((4*cl-2)*256),255,0)
    else:
     return QColor(255,int((4-4*cl)*256),0)

  # if at the point there is something, return it
  def anythingToTip(self,pos):
   point = self.inverseWorldMatrix().map(self.viewportToContents(pos))
   ilist = self.canvas().collisions(point) #QCanvasItemList ilist
   if len(ilist)==0:
    return [None,None]
   tmp_str=""
   for each_item in ilist:
     if each_item.rtti()==POINT_SOURCE_RTTI:
      tmp_str+=" "+each_item.name+" "
      # see if it belongs to a patch
      punit=self.lsm.p_table[each_item.name]
      if punit._patch_name!=None:
       tmp_str+="("+punit._patch_name+")"

   cr=QRect(self.contentsToViewport(self.worldMatrix().map(point)),QSize(len(tmp_str),2))
   return [cr,tmp_str]

  # Explicitly delete the tooltip
  def __del__(self):
    del self.tooltip
    self.tooltip=None


  def display(self):
   #print self.parent.slider1.value()
   # rememeber the number of sources showing
   showcount=0;
   if (self.min_brightness<0):
     abs_min=1e-10
   else: abs_min=self.min_brightness
   # hide all 
   for sname in self.p_tab.keys():
      self.p_tab[sname].hide()

   # check to see if need to limit display of sources
   if self.lsm.display_punits>0:
     pulist=self.lsm.queryLSM(count=self.lsm.display_punits)
   else:
     pulist=self.lsm.queryLSM(all=1)
   for punit in pulist:
       sname=punit.name
       if punit.getType()==POINT_TYPE and\
         math.log(abs(punit.getBrightness(self.default_mode,self.default_freq_index, self.default_time_index))/abs_min)/math.log(abs(self.max_brightness)/abs_min) <self.parent.sliderCut.value()/100.0: 
         pass
       else:
         self.p_tab[sname].show()
         showcount+=1;
 
   tmpstr="%4.3f"%(math.pow(math.e,math.log(abs(self.max_brightness)/abs_min)*self.parent.sliderCut.value()/100.0)*abs_min)
   self.slabel.setText(tmpstr+"/"+str(showcount)) 
   self.canvas().update()  

  # update axes
  def updateAxes(self,x_count,y_count):
    bounds=self.lsm.getBounds()
    if self.axes != None:
     self.axes.gridOff()
     self.axes.hide()
     del self.axes
    self.axes=Axes(self,bounds,x_count,y_count)
    self.xdivs=x_count
    self.ydivs=y_count
 
    if self.grid_on==0:
     self.axes.gridOff()

  # display a new (f,t) with a new I,Q,U,V or app_brightness
  # value
  # type='I','Q','U','V', or 'A' for app_brightness
  def updateDisplay(self,type='A',f_index=0,t_index=0):
   self.default_freq_index=f_index
   self.default_time_index=t_index
   self.default_mode=type
   #print "Update display ",type
   # first, set min,max limits for
   # brightness
   #print "Current limits [%f,%f]"%(self.max_brightness,self.min_brightness)
   self.max_brightness=self.lsm.getMaxBrightness(type,f_index,t_index)
   self.min_brightness=self.lsm.getMinBrightness(type,f_index,t_index)
   
   #print "Updated limits [%f,%f]"%(self.max_brightness,self.min_brightness)
   # next update p-unit table (colours)
   for sname in self.lsm.p_table.keys():
    punit=self.lsm.p_table[sname]
    mytype=punit.getType()
    if mytype==POINT_TYPE or mytype==GAUSS_TYPE :
     # update size and colour both, if pcrosses are displayed 
     self.p_tab[sname].updateDisplayProperties(self.getColor(punit.getBrightness(type,f_index,t_index)), punit.getBrightness(type,f_index,t_index))
    else: #PATCH_TYPE
     self.p_tab[sname].updateDisplayProperties()

   # update indicator
   head=self.lsm.getCurrentFreqTime(self.default_freq_index,self.default_time_index)
   tmpval=stdForm(head[0],'%3.4f')
   headstr="("+tmpval[0]+tmpval[1]+"Hz,"
   tmpval=stdForm(head[1],'%3.4f')
   headstr=headstr+tmpval[0]+tmpval[1]+"s)"
   #headstr="f=%5.4f GHz,t=%5.4f s"%(head[0]/1.0e9,head[1])
   self.zlabel.setText("<font color=\"blue\">"+headstr+"</font>")
   # update legend
   if self.legend!=None:
    self.legend.hide()
    self.showLegend(1)

   self.canvas().update()  

  def showPointSources(self,flag):
    # flag:0==circle(point),1==cross
    if flag==0:
     self.display_point_sources='point'
    elif flag==1:
     self.display_point_sources='cross'
    else:
     self.display_point_sources='pcross'
    for sname in self.lsm.p_table.keys():
      self.p_tab[sname].showType(flag)

  # resets F,T back to 0,0
  # needed when the cells are updated
  def resetFTindices(self):
    self.default_freq_index=0
    self.default_time_index=0
    head=self.lsm.getCurrentFreqTime(self.default_freq_index,self.default_time_index)

    tmpval=stdForm(head[0],'%3.4f')
    headstr="("+tmpval[0]+tmpval[1]+"Hz,"
    tmpval=stdForm(head[1],'%3.4f')
    headstr=headstr+tmpval[0]+tmpval[1]+"s)"
    self.zlabel.setText("<font color=\"blue\">"+headstr+"</font>")


  # create patches from the grid, only if the grid is ON
  def createPatchesFromGrid(self,min_bright=0.0,max_bright=10.0,min_sources=1):
   # create two arrays for x divisions
   # and y divisions and send to the LSM to create 
   # patches
   if self.grid_on==1:
    #print "creating patches from grid"
    stp=(self.ra_max-self.ra_min)/self.xdivs
    x_array=[self.ra_min]
    for ii in range(self.xdivs):
     x_array.append(self.ra_min+(ii+1)*stp)
    stp=(self.dec_max-self.dec_min)/self.ydivs
    y_array=[self.dec_min]
    for ii in range(self.ydivs):
     y_array.append(self.dec_min+(ii+1)*stp)

    retval_arr=self.lsm.createPatchesFromGrid(x_array,y_array,min_bright,max_bright,min_sources)
    #print retval_arr
    if retval_arr != None:
      for retval in retval_arr:
       if retval !=None:
        # successfully created patch
        #print "created patch",retval
        # get the sources of this patch
        punit=self.lsm.getPUnit(retval[0])
        psource_list=punit.getSources()
        # remove these sources from PUnit table on main window
        self.parent.removePUnitRows(psource_list)
        # update the GUI
        self.p_tab[retval[0]]=PatchDisplay(retval[0],self,retval[1],retval[2],\
                 retval[3],retval[4])
        # update PUnit table on main window
        self.parent.insertPUnitRow(retval[0])
      self.canvas().update()

  def showLegend(self,flag):
   """if flag==1, show legend, else hide legend"""
   # get dimensions needed
   [char_w,char_h]=self.getTextDims("10000.000",self.font)
   if flag==1:
    self.canvas().resize(self.canvas().width()+30+char_w,self.canvas().height())
    # get limits from the boundary of main plot
    qp=self.axes.rect.rect()
    print qp.right(),qp.bottom(),qp.top()
    #rr=QCanvasRectangle(qp.right()+self.d2+5,qp.top(),24,qp.bottom()-qp.top(),self.canvas())
    self.legend=Legend(qp.right()+self.d2+5,qp.top(),24,qp.bottom()-qp.top(),\
        self.canvas(),self,"%8.3f")
    self.legend.show()
    self.legend_on=1
   elif flag==0 and self.legend_on==1:
    self.canvas().resize(self.canvas().width()-30-char_w,self.canvas().height())
    self.legend_on=0
    if self.legend !=None:
     self.legend.hide()
     self.legend=None

  # give the text width,height in pixels using the default font
  def getTextDims(self,txt,myfont):
    fm=QFontMetrics(myfont)
    # find width in pixels
    label=QString(txt)
    char_width=fm.width(label)
    char_height=fm.height()
    return (char_width,char_height)

  #select new font
  def chooseFont( self ) :
   ok = 0
   oldfont = QFont( self.font )
   newfont, ok = QFontDialog.getFont(oldfont,self)
   if ok:
    self.font=newfont
    self.axes.updateFont(newfont)
    if self.legend !=None:
     self.legend.updateFont(newfont)


  # save the canvas as a pixmap
  def getPixmap(self):
    pm=QPixmap(self.canvas().width(),self.canvas().height())
    pn=QPainter(pm)
    self.canvas().drawArea(self.canvas().rect(),pn)
    return pm
  # return label as a pixmap
  def createLabel(self,label):
    pmap=QPixmap(len(label)+5,20) 
    painter=QPainter(pmap)
    m=QWMatrix()
    m.rotate(90)
    #painter.begin(pmap)
    painter.setWorldMatrix(m)
    painter.drawText(0,0,label)
    #painter.end() 
    return pmap


  # calculate default dimentions
  def calculateDefaultDims(self):
    # set default padding
    self.pad=5
    # title height and width
    [tw,th]=self.getTextDims("DEFAULT TITLE",self.fonts['title'])
    self.h1=th+self.pad*2
    # ylabel height and width
    [ylw,ylh]=self.getTextDims("Declination",self.fonts['ylabel'])
    # y coords (in radians and degrees)
    [char_w,char_h]=self.getTextDims("10000.000",self.fonts['ycoord']) # radians
    [char_w1,char_h1]=self.getTextDims("66:66:66",self.fonts['ycoord']) # degrees
    ycw=max(char_w,char_w1)
    ych=self.ydivs*char_h # this is calculated to see if we clobber everything
                          # we cannot do anything about it
    self.h2=max(ylw,ych)
    # we will recalculate the above to fit any free space left from 
    # canvas height

    self.w1=ylh+2*self.pad
    self.w2=ycw+2*self.pad

    # x label and coords
    [xlw,xlh]=self.getTextDims("Right Ascension",self.fonts['xlabel'])
    self.h3=xlh+self.pad*2
    [char_w,char_h]=self.getTextDims("10000.000",self.fonts['xcoord']) # radians
    [char_w1,char_h1]=self.getTextDims("66:66:66",self.fonts['xcoord']) # degrees
    self.h4=max(char_h,char_h1)+self.pad*2
     
    xcw=self.xdivs*max(char_w,char_w1)
   
    self.w3=max(xcw,xlw)
    # we will recalculate the above to fit any free space left from 
    # canvas height

    # this is the space required to print trailing part of x coords
    self.w4=max(char_w,char_w1)

    # calculate only if legend is on
    if self.legend_on==0:
     self.w5=self.pad

    if (self.axes_on==0) and (self.obsres_on==0):
     self.w1=self.w2=self.pad
     self.h3=self.h4=self.pad

    if self.title_on==0:
     self.h1=self.pad+5
 
    # get canvas size
    H=self.canvas().height()
    W=self.canvas().width()
    # recalculate drawing area
    self.h2=max(H-self.h1-self.h3-self.h4,self.h2)
    self.w3=max(W-self.w1-self.w2-self.w4-self.w5,self.w3)
 
   
    
  # redraw everything - called when contents have changed
  # or geometry has changed
  def updateCanvas(self):
    # first delete everything on the canvas
    ilist=self.canvas().allItems()
    for itm in ilist:
     if itm:
      itm.setCanvas(None)
      del itm
    bounds=self.lsm.getBounds()

    # setup 
    self.h1=self.h2=self.h3=self.h4=0
    self.w1=self.w2=self.w3=self.w4=self.w5=0
    self.pad=0 # padding
    # now calculate the real values
    self.calculateDefaultDims()

    self.d1=self.w1+self.w2
    self.d2=self.w4
    self.d3=self.h1
    self.d4=self.h3+self.h4

    # limits in real coordinates
    self.ra_min=bounds['min_RA']
    self.ra_max=bounds['max_RA']
    self.dec_min=bounds['min_Dec']
    self.dec_max=bounds['max_Dec']

    # projector
    self.proj=transform.Projector(self.ra0,self.dec0,0)
    if self.proj_on==0: self.proj.Off() # turn off projection
    ### now recalculate the bounds according to the projection params
    left_down=self.proj.give_limits(self.ra_min,self.ra_max,self.dec_min,self.dec_max)
    self.x_min=left_down[0]
    self.x_max=left_down[1]
    self.y_min=left_down[2]
    self.y_max=left_down[3]

    # sanity check
    if self.x_min==self.x_max:
      self.x_max=self.x_min+0.1
    if self.y_min==self.y_max:
      self.y_max=self.y_min+0.1
    # canvas size
    H=self.canvas().height()
    W=self.canvas().width()
    # scales
    if (self.x_max != self.x_min):
     self.x_scale=(W-self.d1-self.d2)/(self.x_max-self.x_min)
    else: 
     self.x_scale=1.0

    if (self.y_max != self.y_min):
     self.y_scale=(H-self.d3-self.d4)/(self.y_max-self.y_min)
    else:
     self.y_scale=1.0

    ############ re-create p-unit list
    # table for all p-units plots on canvas
    self.p_tab={}
    # plot all p-units/sources
    for sname in self.lsm.p_table.keys():
     punit=self.lsm.p_table[sname]
     if punit.getType()==POINT_TYPE:
      self.p_tab[sname]=PointSourceDisplay(sname,self)
     elif punit.getType()==GAUSS_TYPE:
      self.p_tab[sname]=GaussianSourceDisplay(sname,self)
     else:
      # we have a patch
      retval=punit.getLimits()
      self.p_tab[sname]=PatchDisplay(sname,self,retval[0],retval[1],\
                 retval[2],retval[3])


    ############ re-create axes/grid
    if self.axes_on==1:
     self.axes=Axes(self,bounds,self.xdivs,self.ydivs)
    else:
     self.axes=None
    if self.grid_on==1:
     self.axes.gridOn()
    else:
     self.axes.gridOff()

    ########### create title
    if self.title_on==1:
     self.title=FontHorizImage("Title",self.canvas(),self.fonts['title'],self.pad)
     xys=self.globalToLocal(self.x_max, self.y_max)
     self.title.move(xys[0]+(self.w3-self.title.width())/2,xys[1]-self.h1)
     self.title.setZ(0)
     self.title.show()
    else:
     self.title=None

    ############ create zoom-window
    self.zwindow=ZWindow(self.canvas())
    #self.zwindow.hide()
    self.zoom_status=GUI_ZOOM_NONE
    # use old transformation matrices zoom in/out
    #self.tmstack=None
    
    self.__moving=0
    self.__moving_start=0
    self.canvas().update()


   
  # enlarge the sources by a factor of 2
  def enlarge_display(self):
    for sname in self.lsm.p_table.keys():
       punit=self.lsm.p_table[sname]
       if punit.getType()==POINT_TYPE:
          # update size and colour both, if pcrosses are displayed 
          self.p_tab[sname].resizeItems(2)
       else: #PATCH_TYPE
          pass

    self.canvas().update()

  # shrink the sources by a factor of 2
  def shrink_display(self):
    for sname in self.lsm.p_table.keys():
       punit=self.lsm.p_table[sname]
       if punit.getType()==POINT_TYPE:
          # update size and colour both, if pcrosses are displayed 
          self.p_tab[sname].resizeItems(0.5)
       else: #PATCH_TYPE
          pass

    self.canvas().update()

##################################################################
