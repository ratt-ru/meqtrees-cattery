#!/usr/bin/python

##########################################################
# inner classes needed by the canvas
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

import Timba.array

##########################################################################
#################################################################
# class for drawing point (and elliptic - extended) sources
class Circle(QCanvasEllipse):
  def __init__(self,center_x,center_y,radius,color,name,canvas):
    QCanvasEllipse.__init__(self,radius*2,radius*2,canvas)
    self.name=name
    self.canvas=canvas
    self.color=color
    self.radius=radius
    self.setBrush(QBrush(self.color))
    self.setX(center_x)
    self.setY(center_y)
    self.myRTTI=POINT_SOURCE_RTTI

  def rtti(self):
    return self.myRTTI

  def updateColor(self,color):
   self.color=color
   self.setBrush(QBrush(self.color))

  def getDims(self):
   return self.radius


#################################################################
# class for drawing point (and elliptic - extended) sources
class Ellipse1(QCanvasEllipse):
  def __init__(self,center_x,center_y,major,minor,PA,color,name,canvas):
    QCanvasEllipse.__init__(self,major,minor,canvas)
    self.name=name
    self.canvas=canvas
    self.color=color
    self.major=major
    self.minor=minor
    self.setBrush(QBrush(self.color))
    self.setX(center_x)
    self.setY(center_y)
    self.myRTTI=POINT_SOURCE_RTTI

  def rtti(self):
    return self.myRTTI

  def updateColor(self,color):
   self.color=color
   self.setBrush(QBrush(self.color))

  def getDims(self):
   return max(self.major,minor)

class Ellipse(QCanvasPolygonalItem):
  def __init__(self,center_x,center_y,major,minor,PA,color,name,canvas):
    QCanvasPolygonalItem.__init__(self,canvas)
    self.name=name
    self.canvas=canvas
    self.color=color
    self.major=major
    self.minor=minor
    self.PA=PA
    self.myRTTI=POINT_SOURCE_RTTI
    self.setPen(QPen(self.color))
    self.setBrush(QBrush(self.color, Qt.Dense3Pattern))
    self.pen().setWidth(0) # have no effect because polygonal item does not use it
    self.move(center_x,center_y)

  ## destructor must call hide
  def __del__(self):
    self.hide()

  def rtti(self):
    return self.myRTTI
 
  def updateColor(self,color):
   self.color=color

  def getDims(self):
   return (self.major,self.minor)

  ## area of the bounding rectangle
  def areaPoints(self):
    pa=QPointArray(4)
    pa.setPoint(0,QPoint(self.x()+self.major/2,self.y()+self.major/2))
    pa.setPoint(1,QPoint(self.x()-self.major/2,self.y()+self.major/2))
    pa.setPoint(2,QPoint(self.x()-self.major/2,self.y()-self.major/2))
    pa.setPoint(3,QPoint(self.x()+self.major/2,self.y()-self.major/2))
    return pa

  ## draw the ellipse
  def drawShape(self,p):
   plist=self.calculatePoints(20)
   (X0,Y0)=plist.pop() 
   Xs=X0
   Ys=Y0 
   for (X,Y) in plist:
    pa=QPointArray(3)
    pa.setPoint(0,QPoint(self.x(),self.y()))
    pa.setPoint(1,QPoint(X0,Y0))
    pa.setPoint(2,QPoint(X,Y))
    p.drawPolygon(pa)
    #p.drawLine(X0,Y0,X,Y)
    X0=X
    Y0=Y
   pa=QPointArray(3)
   pa.setPoint(0,QPoint(self.x(),self.y()))
   pa.setPoint(1,QPoint(X0,Y0))
   pa.setPoint(2,QPoint(Xs,Ys))
   p.drawPolygon(pa)
   #p.drawLine(X0,Y0,Xs,Ys)


  def updateColor(self,color):
   self.color=color
   self.setBrush(QBrush(self.color, Qt.Dense3Pattern))
   self.setPen(QPen(self.color))

  def updateSize(self,scale):
    pass

  def calculatePoints(self,np):
    dtheta=math.pi*2/np
    dang=0
    cosa=math.cos(self.PA)
    sina=math.sin(self.PA)
    plist=[]
    for ci in range(np):
     theta=dang+ci*dtheta
     cosT=math.cos(theta)
     sinT=math.sin(theta)
     X=self.major/2*cosT*cosa-self.minor/2*sinT*sina
     Y=self.major/2*cosT*sina+self.minor/2*sinT*cosa
     plist.append((X+self.x(),Y+self.y()))

    return plist


##################################################################
class pCross(QCanvasPolygon):
  def __init__(self,center_x,center_y,w,l,color,name,canvas):
    QCanvasPolygon.__init__(self,canvas)
    self.name=name
    self.canvas=canvas
    self.color=color
    self.width=w
    self.length=l
    self.setBrush(QBrush(self.color))
    self.myRTTI=POINT_SOURCE_RTTI
    pa=QPointArray(12)
    pa.setPoint(0, QPoint(w,w))
    pa.setPoint(1, QPoint(w,l))
    pa.setPoint(2, QPoint(-w,l))
    pa.setPoint(3, QPoint(-w,w))
    pa.setPoint(4, QPoint(-l,w))
    pa.setPoint(5, QPoint(-l,-w))
    pa.setPoint(6, QPoint(-w,-w))
    pa.setPoint(7, QPoint(-w,-l))
    pa.setPoint(8, QPoint(w,-l))
    pa.setPoint(9, QPoint(w,-w))
    pa.setPoint(10, QPoint(l,-w))
    pa.setPoint(11, QPoint(l,w))
    self.setPoints(pa)
    self.myRTTI=POINT_SOURCE_RTTI
    self.setBrush(QBrush(color))
    self.move(center_x,center_y)

  def rtti(self):
    return self.myRTTI
 
  def updateColor(self,color):
   self.color=color
   self.setBrush(QBrush(self.color))

  def getDims(self):
   return (self.width,self.length)


##################################################################
class tCross(QCanvasPolygonalItem):
  def __init__(self,center_x,center_y,w,l,color,name,canvas):
    QCanvasPolygonalItem.__init__(self,canvas)
    self.name=name
    self.canvas=canvas
    self.color=color
    self.width=w
    self.length=abs(l)
    self.myRTTI=POINT_SOURCE_RTTI
    self.setPen(QPen(self.color))
    self.setBrush(QBrush(self.color))
    self.pen().setWidth(w)
    self.move(center_x,center_y)

  ## destructor must call hide
  def __del__(self):
    self.hide()

  def rtti(self):
    return self.myRTTI
 
  def updateColor(self,color):
   self.color=color

  def getDims(self):
   return (self.width,self.length)

  ## area of the bounding rectangle
  def areaPoints(self):
    pa=QPointArray(4)
    pa.setPoint(0,QPoint(self.x()+self.length,self.y()+self.length))
    pa.setPoint(1,QPoint(self.x()-self.length,self.y()+self.length))
    pa.setPoint(2,QPoint(self.x()-self.length,self.y()-self.length))
    pa.setPoint(3,QPoint(self.x()+self.length,self.y()-self.length))
    return pa

  ## draw the cross
  def drawShape(self,p):
   # draw the lines
   length=self.length
   if length<1: length=1
   p.drawLine(self.x()-length,self.y(),self.x()+length,self.y())
   p.drawLine(self.x(),self.y()-length,self.x(),self.y()+length)

  def updateColor(self,color):
   self.color=color
   self.setBrush(QBrush(self.color))
   self.setPen(QPen(self.color))

  def updateSize(self,scale):
   self.length=self.length*scale
   if self.length<1: self.length=1



#################################################################
class DynamicTip( QToolTip ):
    def __init__( self, parent ):
        QToolTip.__init__( self, parent )

    def maybeTip( self, pos ):
        rs =QToolTip(self).parentWidget().anythingToTip(pos)
        if rs[0]==None:
          return
        
        QToolTip(self).tip( rs[0], rs[1] )
        self.parentWidget().canvas().update()

####################################################
# class for zoom window
class ZWindow(QCanvasRectangle):
  def __init__(self,canvas):
    QCanvasRectangle.__init__(self,canvas) 
    self.left_x=0
    self.left_y=0
    self.right_x=0
    self.right_y=0
  def setUpperLeft(self,left_x,left_y):
    self.left_x=left_x
    self.left_y=left_y
    # move to this place
    self.move(left_x,left_y)
    # set size to zero
    self.setSize(0,0)
  def setLowerRight(self,right_x,right_y):
    self.right_x=right_x
    self.right_y=right_y
    # resize the rectangle to get correct size
    self.setSize(right_x-self.left_x,right_y-self.left_y)

  def getRect(self):
   return QRect(QPoint(self.left_x,self.left_y),QPoint(self.right_x,self.right_y))



##########################################################
class Axes:
  """ canvas_view: QCanvasview object
      bounds: dict with min_RA,max_RA,min_Dec,max_Dec values
      x_ticks,y_ticks: number of division in x and y direction
  """
  def __init__(self,canvas_view,bounds,x_ticks,y_ticks):
    self.cview=canvas_view
    # draw boundary
    self.rect=QCanvasRectangle(self.cview.w1+self.cview.w2, self.cview.h1,self.cview.w3+self.cview.pad,self.cview.h2,self.cview.canvas())
    self.rect.show()

    # grid
    self.grid=[]

    # font for axes ticks
    self.axfont=QFont(canvas_view.font)
    #self.axfont.setPointSize(10)
    #self.axfont.setWeight( QFont.Bold )
    #self.axfont.setUnderline( True )


    # draw X axis
    self.xax=[]
    self.xax_text=[]
    self.xax_degtext=[]
    ll=self.getTicks(bounds['min_RA'],bounds['max_RA'],x_ticks,'%1.3f')
    xgridvals=[]
    for tck in ll:
     xval=tck[0]
     xstr=tck[1]
     degstr=tck[2]+":"+tck[3]+":"+tck[4]
     # draw line
     ff=QCanvasLine(self.cview.canvas())
     xys=self.cview.globalToLocal(xval,bounds['min_Dec'])
     # change y coordinate
     xys[1]=self.cview.h1+self.cview.h2 
     ff.setPoints(xys[0],xys[1],xys[0],xys[1]-10)
     ff.show()
     self.xax.append(ff)
     # text in radians 
     rt=QCanvasText(self.cview.canvas())
     rt.setText(xstr)
     rt.setFont( self.axfont )
     rt.move(xys[0]-25,xys[1]+9)
     rt.setZ(0)
     rt.hide()
     self.xax_text.append(rt)
     # text in degrees
     dt=QCanvasText(self.cview.canvas())
     dt.setText(degstr)
     dt.setFont( self.axfont )
     dt.move(xys[0]-25,xys[1]+9)
     dt.setZ(0)
     dt.hide()
     self.xax_degtext.append(dt)
     if self.cview.default_coords=='rad':
      rt.show()
     else:
      dt.show()

     xgridvals.append(xval)

    # draw Y axis
    self.yax=[]
    self.yax_text=[]
    self.yax_degtext=[]
    ll=self.getTicks(bounds['min_Dec'],bounds['max_Dec'],y_ticks,'%1.4f','y')
    ygridvals=[]
    for tck in ll:
     yval=tck[0]
     xstr=tck[1]
     degstr=tck[2]+"."+tck[3]+"."+tck[4]
     # draw line
     ff=QCanvasLine(self.cview.canvas())
     xys=self.cview.globalToLocal(bounds['max_RA'],yval)
     # change x coordinate
     xys[0]=self.cview.w1+self.cview.w2
     ff.setPoints(xys[0],xys[1],xys[0]+10,xys[1])
     ff.show()
     self.yax.append(ff)
     # text in radians 
     rt=QCanvasText(self.cview.canvas())
     rt.setText(xstr)
     rt.setFont( self.axfont )
     fm=QFontMetrics(self.axfont)
     char_width=fm.width(rt.text())
     # fixme 10 should be equal to tick width
     rt.move(xys[0]-char_width-10,xys[1]-9)
     rt.setZ(0)
     rt.hide()
     self.yax_text.append(rt)
     # text in degrees
     dt=QCanvasText(self.cview.canvas())
     dt.setText(degstr)
     dt.setFont( self.axfont )
     char_width=fm.width(dt.text())
     dt.move(xys[0]-char_width-10,xys[1]-9)
     dt.setZ(0)
     dt.hide()
     self.yax_degtext.append(dt)
     if self.cview.default_coords=='rad':
      rt.show()
     else:
      dt.show()

     # draw grid lines
     ygridvals.append(yval)


    # create two arrays
    xcarr=Timba.array.zeros((x_ticks+1,y_ticks+1))
    ycarr=Timba.array.zeros((x_ticks+1,y_ticks+1))
    ii=0
    for xc in xgridvals:
     jj=0
     for yc in ygridvals:
       # do the projection
       xys=self.cview.globalToLocal(xc,yc)
       xcarr[ii,jj]=xys[0]
       ycarr[ii,jj]=xys[1]
       jj=jj+1
     ii=ii+1
    # draw grid lines
    for ii in range(0,x_ticks+1):
     for jj in range(0,y_ticks):
       ff=QCanvasLine(self.cview.canvas())
       ff.setPoints(xcarr[ii,jj],ycarr[ii,jj],xcarr[ii,jj+1],ycarr[ii,jj+1])
       ff.setPen(QPen(QColor(0,0,255)))
       self.grid.append(ff)
       ff.show()
    for jj in range(0,y_ticks+1):
     for ii in range(0,x_ticks):
       ff=QCanvasLine(self.cview.canvas())
       ff.setPoints(xcarr[ii,jj],ycarr[ii,jj],xcarr[ii+1,jj],ycarr[ii+1,jj])
       ff.setPen(QPen(QColor(0,0,255)))
       self.grid.append(ff)
       ff.show()


    # draw axes labels
    self.ylabel=FontVertImage("Declination",self.cview.canvas(),self.cview.fonts['ylabel'],self.cview.pad)
    xys=[self.cview.w1+self.cview.w2,self.cview.h1+self.cview.h2]
    self.ylabel.move(xys[0]-self.cview.d1,xys[1]-(self.cview.canvas().height()-self.cview.d3-self.cview.d4)/2-self.ylabel.height()/2)
    self.ylabel.setZ(0)
    self.ylabel.show()

    self.xlabel=FontHorizImage("Right Ascension",self.cview.canvas(),self.cview.fonts['xlabel'],self.cview.pad)
    self.xlabel.move(xys[0]+(self.cview.canvas().width()-self.cview.d1-self.cview.d2)/2-self.xlabel.width()/2,xys[1]+self.cview.h3)
    self.xlabel.setZ(0)
    self.xlabel.show()

  # return tickmarks for the axis
  # returns [coordinate value, rad_value, hr, min, sec]
  def getTicks(self,start_x,end_x,divisions=10,format='%e',axis='x'):
     inc=(end_x-start_x)/float(divisions)
     output=[]
     # start point
     st=format%start_x
     # degrees
     if axis=='x':
      tmpval=radToRA(start_x)
     else:
      tmpval=radToDec(start_x)
     output.append([start_x, st, str(tmpval[0]), str(tmpval[1]), str(tmpval[2])])
     for i in range(divisions):
      val=start_x+(i+1)*inc
      st=format%val
      if axis=='x':
       tmpval=radToRA(val)
      else:
       tmpval=radToDec(val)
      output.append([val, st, str(tmpval[0]), str(tmpval[1]), str(tmpval[2])])
      #output.append([val,st])

     return output

  # hide axes 
  def hide(self):
   for i in range(len(self.xax)):
    self.xax[i].hide()
   for i in range(len(self.xax_text)):
    self.xax_text[i].hide()
    self.xax_degtext[i].hide()
   for i in range(len(self.yax)):
    self.yax[i].hide()
   for i in range(len(self.yax_text)):
    self.yax_text[i].hide()
    self.yax_degtext[i].hide()

  # show axes 
  def show(self):
   for i in range(len(self.xax)):
    self.xax[i].show()
   for i in range(len(self.xax_text)):
    self.xax_text[i].show()
    self.xax_degtext[i].show()
   for i in range(len(self.yax)):
    self.yax[i].show()
   for i in range(len(self.yax_text)):
    self.yax_text[i].show()
    self.yax_degtext[i].show()


  # hide the grid
  def gridOff(self):
    for i in range(len(self.grid)):
      self.grid[i].hide()
  # show the grid
  def gridOn(self):
    for i in range(len(self.grid)):
      self.grid[i].show()

  def updateFont(self,newfont):
    self.axfont=newfont
    for i in range(len(self.xax_text)):
     self.xax_text[i].setFont( self.axfont )
     self.xax_degtext[i].setFont( self.axfont )
    for i in range(len(self.yax_text)):
     self.yax_text[i].setFont( self.axfont )
     self.yax_degtext[i].setFont( self.axfont )

  # display coords either in radians (rad) or in degrees (deg)
  def switchCoords(self,coords='rad'):
   if coords=='rad':
    for i in range(len(self.xax_text)):
     self.xax_text[i].show()
     self.xax_degtext[i].hide()
    for i in range(len(self.yax_text)):
     self.yax_text[i].show()
     self.yax_degtext[i].hide()
   else:
    for i in range(len(self.xax_text)):
     self.xax_text[i].hide()
     self.xax_degtext[i].show()
    for i in range(len(self.yax_text)):
     self.yax_text[i].hide()
     self.yax_degtext[i].show()
   
#############################################################
class PointSourceDisplay:
 def __init__(self,name,parent):
  self.name=name
  self.cview=parent
  # get corresponding PUnit
  punit=self.cview.lsm.p_table[self.name]
  #print self.name
  # get coords 
  xys=self.cview.globalToLocal(punit.sp.getRA(),punit.sp.getDec())
  self.x=xys[0]
  self.y=xys[1]
  #print xys
  self.cross=self.addCross(xys[0],xys[1],1,5,self.cview.getColor(punit.getBrightness()),self.name,punit.getBrightness())
  #print "Brightness: min,max, absmin, current",self.cview.min_brightness,self.cview.max_brightness, self.cview.abs_min_brightness, punit.getBrightness()
  br=abs(punit.getBrightness())
  if self.cview.max_brightness!=self.cview.min_brightness:
   rat=abs(float(self.cview.max_brightness))/self.cview.abs_min_brightness
  else:
   rat=1.1
  if br >0.0:
   length=int(math.log(br/self.cview.abs_min_brightness)/math.log(rat)*10)
   if length==0:
    length=1 # one pixel
  else:
   length=1
  self.pcross=self.addCross(xys[0],xys[1],1,length,self.cview.getColor(punit.getBrightness()),self.name,punit.getBrightness())
  self.circle=self.addCircle(xys[0],xys[1],2,self.cview.getColor(punit.getBrightness()),self.name,punit.getBrightness())

  self.show()

 # circle
 def addCircle(self,center_x,center_y,radius,color,name,z_value):
  i = Circle(center_x,center_y,radius,color,name,self.cview.canvas())
  i.setZ(z_value)
  return i

 # cross - polygon
 def addCross(self,center_x,center_y,w,l,color,name,z_value):
  c=tCross(center_x,center_y,w,l,color,name,self.cview.canvas())
  c.setZ(z_value)
  return c

 def show(self):
  if self.cview.display_point_sources=='cross':
   self.cross.show()
  elif self.cview.display_point_sources=='point':
   self.circle.show()
  else:
   self.pcross.show()


 def hide(self):
  if self.cview.display_point_sources=='cross':
   self.cross.hide()
  elif self.cview.display_point_sources=='point':
   self.circle.hide()
  else:
   self.pcross.hide()

 def hideAll(self):
  self.cross.hide()
  self.circle.hide()
  self.pcross.hide()

 def showType(self,flag):
  if flag==0:
   self.cross.hide()
   self.pcross.hide()
   self.circle.show()
  elif flag==1:
   self.circle.hide()
   self.pcross.hide()
   self.cross.show()
  elif flag==2:
   self.circle.hide()
   self.pcross.show()
   self.cross.hide()


 def updateDisplayProperties(self,newcolor,new_value):
  self.circle.updateColor(newcolor)
  self.cross.updateColor(newcolor)
  # Neet to adjust the size of pcross as well
  # instead of adjusting current pcross, recreate a new one
  # self.pcross.updateColor(newcolor)
  self.pcross.hide()
  del self.pcross
  if new_value==0 or\
    self.cview.max_brightness==0:
   length=0
  else:
   new_value=abs(new_value)
   if self.cview.max_brightness!=self.cview.min_brightness:
     rat=abs(self.cview.max_brightness)/self.cview.abs_min_brightness
   else:
     rat=1.1
   if new_value>0.0:
     length=int(math.log(new_value/self.cview.abs_min_brightness)/math.log(rat)*10)
     if length==0:
       length=1 # one pixel
   else:
     length=1

  self.pcross=self.addCross(self.x,self.y,1,length,newcolor,self.name,new_value)
 
  if self.cview.display_point_sources=='cross':
   self.cross.show()
  elif self.cview.display_point_sources=='point':
   self.circle.show()
  else:
   self.pcross.show()


 # shrink/enlarge each object by the given scale
 def resizeItems(self,new_scale):
  # Neet to adjust the size of pcross as well
  # instead of adjusting current pcross, recreate a new one
  # self.pcross.updateColor(newcolor)
  oldc=self.pcross
  (w,l)=oldc.getDims()
  new_l=l*new_scale
  if (int(new_l)>=1): 
    oldc.hide()
    # if the new with is less than 1, change to a skeleton
    self.pcross=tCross(self.x,self.y,w,new_l,oldc.color,self.name,self.cview.canvas())
    self.pcross.setZ(oldc.z())
    del oldc 
  else:
    # do nothing
    return

  if self.cview.display_point_sources=='cross':
   self.cross.show()
  elif self.cview.display_point_sources=='point':
   self.circle.show()
  else:
   self.pcross.show()


################################################################
#############################################################
class PatchDisplay:
 def __init__(self,name,parent,x_min,y_min,x_max,y_max):
  self.name=name
  self.cview=parent
  self.x_min=x_min
  self.x_max=x_max
  self.y_min=y_min
  self.y_max=y_max
  # create a rectangle
  xys=self.cview.globalToLocal(self.x_max,self.y_max)
  topLeft=QPoint(xys[0],xys[1])
  xys=self.cview.globalToLocal(self.x_min,self.y_min)
  bottomRight=QPoint(xys[0],xys[1])
  rectangle=QRect(topLeft,bottomRight)
  self.rect=QCanvasRectangle(rectangle,self.cview.canvas())
  self.rect.setPen(QPen(QColor(255,0,0)))
  self.rect.setZ(0)

  try:
   # get vellsets if any
   punit=self.cview.lsm.p_table[self.name]
   lims=punit.sp.getValueSize(self.cview.default_mode,\
    self.cview.default_freq_index,\
    self.cview.default_time_index)
 
   # the return type should be a Timba.array
   byarr=punit.sp.getValue(self.cview.default_mode,\
    self.cview.default_freq_index,\
    self.cview.default_time_index)
   # create an image
   self.image=ImageItem(self.createArrayImage(byarr,lims[0],lims[1]),\
           self.cview.canvas(),self)
   self.image.move(self.rect.x(),self.rect.y())
   self.image.setZ(-1)
   # change RTTI to a custom value
   self.image.setRTTI(PATCH_IMAGE_RTTI)
  except:
   self.image=None
   pass
  self.show()

 def hide(self):
  self.rect.hide()

 def show(self):
  self.rect.show()
  if self.image !=None:
   self.image.show()

 def hideAll(self):
  self.rect.hide()

 def showAll(self):
  self.rect.show()
  if self.image !=None:
   self.image.show()

 def showType(self,flag):
  self.show()

 def updateDisplayProperties(self):
  try:
   punit=self.cview.lsm.p_table[self.name]
   lims=punit.sp.getValueSize(self.cview.default_mode,\
    self.cview.default_freq_index,\
    self.cview.default_time_index)
 
   # the return type should be a Timba.array
   byarr=punit.sp.getValue(self.cview.default_mode,\
    self.cview.default_freq_index,\
    self.cview.default_time_index)
   # create an image
   self.image=ImageItem(self.createArrayImage(byarr,lims[0],lims[1]),\
           self.cview.canvas())
   self.image.move(self.rect.x(),self.rect.y())
   self.image.setZ(-1)
  except:
   pass

  # create an Image from the given Timba.array
 # size x_dim by y_dim
 def createArrayImage(self,narray,x_dim,y_dim):
   # first find min,max values
   minval=narray.min()
   maxval=narray.max()
   #print "array size is %d,%d with values %f,%f"%(x_dim,y_dim,minval,maxval)
   # create an Image of the size of this array
   # create image from this array, 32 bits depth, 2^24 colours
   im=QImage(x_dim,y_dim,32)
   # fill image with White
   im.fill(qRgb(255,255,255))
   # set pixel values for only non-zero elements
   [nz_x,nz_y]=Timba.array.nonzero(narray)
   for ci in range(len(nz_x)):
     cl=self.getRGB(narray[nz_x[ci]][nz_y[ci]]-minval,maxval)
     # flip the image in the y direction
     im.setPixel(nz_x[ci],y_dim-1-nz_y[ci],cl)
     #print "value %f at %d,%d"%(narray[nz_x[ci]][nz_y[ci]],nz_x[ci],nz_y[ci])

   # resize image to fit the rectangle
   im2=im.smoothScale(self.rect.width(),self.rect.height())
   #print "created image size %d,%d"%(im2.width(),im2.height())
   return im2




 # return a colour (RGB) according to brightness
 def getRGB(self,z,max_brightness):
  if max_brightness==0:
    return qRgb(1,1,1) # no color : when Q,U,V is zero

  cl=z/max_brightness # normalized in [0,1] 
  if cl < 0.25:
    return qRgb(0,int(cl*256*4),255)
  elif cl < 0.5:
    return qRgb(0,255,int((2-cl*4)*256))
  elif cl < 0.75:
    return qRgb(int((4*cl-2)*256),255,0)
  else:
    return qRgb(255,int((4-4*cl)*256),0)

##########################################################################
##########################################################################
class GaussianSourceDisplay:
 def __init__(self,name,parent):
  self.name=name
  self.cview=parent
  # get corresponding PUnit
  punit=self.cview.lsm.p_table[self.name]
  #print self.name
  # get coords 
  xys=self.cview.globalToLocal(punit.sp.getRA(),punit.sp.getDec())
  self.x=xys[0]
  self.y=xys[1]
  # get major, minor, PA values
  elist=punit.getExtParms()
  eX=elist[0]
  eY=elist[1]
  eP=elist[2]
  # take out scale
  if self.cview.y_scale > self.cview.x_scale:
   eY=eY/(self.cview.y_scale)*self.cview.x_scale
  else:
   eX=eX/(self.cview.x_scale)*self.cview.y_scale

  ## no need to do this anymore, because we store the true value
  # flip sign or PA
  #eP=-eP

  xys1=self.cview.globalToLocal(punit.sp.getRA()+eX/2,punit.sp.getDec()+eY/2)
  major=int(abs(xys[0]-xys1[0]))*2
  minor=int(abs(xys[1]-xys1[1]))*2
  # set a default minimum value
  if major<2: major=2
  if minor<2: minor=2
  
  self.ellipse=self.addCircle(xys[0],xys[1],major,minor,eP,self.cview.getColor(punit.getBrightness()),self.name,punit.getBrightness())

  self.show()

 # circle
 def addCircle(self,center_x,center_y,major,minor,PA,color,name,z_value):
  i = Ellipse(center_x,center_y,major,minor,PA,color,name,self.cview.canvas())
  i.setZ(z_value)
  return i

 def show(self):
   self.ellipse.show()

 def hide(self):
  self.ellipse.hide()

 def hideAll(self):
  self.ellipse.hide()

 def showType(self,flag):
  self.ellipse.show()

 def updateDisplayProperties(self,newcolor,new_value):
  self.ellipse.updateColor(newcolor)


 # shrink/enlarge each object by the given scale
 def resizeItems(self,new_scale):
  pass

################################################################

###############################################################
class ImageItem(QCanvasRectangle):
    def __init__(self,img,canvas,parent=None):
        QCanvasRectangle.__init__(self,canvas)
        self.imageRTTI=984376
        self.image=img
        self.parent=parent
        self.pixmap=QPixmap()
        self.setSize(self.image.width(), self.image.height())
        self.pixmap.convertFromImage(self.image, Qt.OrderedAlphaDither);

    def setRTTI(self,rtti):
     self.imageRTTI=rtti

    def rtti(self):
        return self.imageRTTI

    def hit(self,p):
        ix = p.x()-self.x()
        iy = p.y()-self.y()
        if not self.image.valid( ix , iy ):
            return False
        self.pixel = self.image.pixel( ix, iy )
        return  (qAlpha( self.pixel ) != 0)

    def drawShape(self,p):
        p.drawPixmap( self.x(), self.y(), self.pixmap )

#################################################################
#produce text rotated by 90
class FontVertImage(QCanvasRectangle):
    def __init__(self,label,canvas,font,v_space):
        QCanvasRectangle.__init__(self,canvas)
        self.imageRTTI=984376
        self.label=label
        self.font=font
        fm=QFontMetrics(self.font)
        margin=20
        # find width in pixels
        char_width=fm.width(self.label)
        char_height=fm.height()
        #v_space=2
        self.pixmap=QPixmap(v_space+char_height+v_space,char_width+2*margin)
        self.pixmap.fill(QColor(255,255,255))
        painter=QPainter(self.pixmap)
        m=QWMatrix() 
        m.rotate(-90)
        #m.scale(2,2)
        painter.setWorldMatrix(m)
        painter.setFont(self.font)
        tmp_str=QString(self.label)
        painter.drawText(-char_width-margin,char_height,QString(self.label))
        painter.end()
        self.setSize(self.pixmap.width(), self.pixmap.height())

    def rtti(self):
        return self.imageRTTI

    def drawShape(self,p):
        p.drawPixmap(self.x(),self.y(), self.pixmap )

#################################################################
#produce normal text 
class FontHorizImage(QCanvasRectangle):
    def __init__(self,label,canvas,font,v_space):
        QCanvasRectangle.__init__(self,canvas)
        self.imageRTTI=984376
        self.label=label
        self.font=font
        fm=QFontMetrics(self.font)
        margin=20
        #v_space=2
        # find width in pixels
        char_width=fm.width(self.label)
        char_height=fm.height()
        self.pixmap=QPixmap(char_width+2*margin, v_space+char_height+v_space)
        self.pixmap.fill(QColor(255,255,255))
        painter=QPainter(self.pixmap)
        painter.setFont(self.font)
        tmp_str=QString(self.label)
        painter.drawText(margin,char_height/2+2*v_space,QString(self.label))
        painter.end()
        self.setSize(self.pixmap.width(), self.pixmap.height())

    def rtti(self):
        return self.imageRTTI

    def drawShape(self,p):
        p.drawPixmap(self.x(),self.y(), self.pixmap )


##################################################################
# class to draw legend colourbar
class Legend:
    def __init__(self,left,top,width,height,canvas,cview,format="%8.3f"):
       self.cview=cview
       self.canvas=canvas
       self.rect=QCanvasRectangle(left,top,width,height,canvas)
       self.top=top
       self.left=left
       self.width=width
       self.height=height
       self.format=format
   
       # fill this with coloured rectangles
       self.ncolors=10
       # height of each rectangle
       h=height/self.ncolors
       zheight=(self.cview.max_brightness-self.cview.min_brightness)/float(self.ncolors)
       self.rts=[]
       self.txt=[]
       for ci in range(1,self.ncolors):
         y=top+height-h*ci
         rt=QCanvasRectangle(left,y,width,h,canvas)
         zlevel=zheight*ci+self.cview.min_brightness
         rt.setBrush(QBrush(self.cview.getColor(zlevel)))
         self.rts.append(rt)
         # ad label
         dt=QCanvasText(canvas)
         dt.setText(format%zlevel)
         dt.setFont(self.cview.axes.axfont)
         dt.move(left+width,y-3)
         dt.setZ(0)
         self.txt.append(dt)

 
       # the last rectangle
       h=y-top
       y=top
       rt=QCanvasRectangle(left,y,width,h,canvas)
       zlevel=zheight*(ci+1)+self.cview.min_brightness
       rt.setBrush(QBrush(self.cview.getColor(zlevel)))
       self.rts.append(rt)
       dt=QCanvasText(canvas)
       dt.setText(format%zlevel)
       dt.setFont(self.cview.axes.axfont)
       dt.move(left+width,y-3)
       dt.setZ(0)
       self.txt.append(dt)

       # last, text at zero level
       zlevel=self.cview.min_brightness
       dt=QCanvasText(canvas)
       dt.setText(format%zlevel)
       y=top+height
       dt.setFont(self.cview.axes.axfont)
       dt.move(left+width,y-3)
       dt.setZ(0)
       self.txt.append(dt)




    def show(self):
      for cr in self.rts:
         cr.show()
      self.rect.show()
      for cr in self.txt:
         cr.show()


    def hide(self):
      for cr in self.rts:
         cr.hide()
      self.rect.hide()
      for cr in self.txt:
         cr.hide()
 

    def update(self):
       # height of each rectangle
       h=self.height/self.ncolors
       zheight=(self.cview.max_brightness-self.cview.min_brightness)/float(self.ncolors)
       self.rts=[]
       self.txt=[]
       for ci in range(1,self.ncolors):
         y=self.top+self.height-h*ci
         rt=QCanvasRectangle(self.left,y,self.width,h,self.canvas)
         zlevel=zheight*ci+self.cview.min_brightness
         rt.setBrush(QBrush(self.cview.getColor(zlevel)))
         self.rts.append(rt)
         # ad label
         dt=QCanvasText(self.canvas)
         dt.setText(self.format%zlevel)
         dt.setFont(self.cview.axes.axfont)
         dt.move(self.left+self.width,y-3)
         dt.setZ(0)
         self.txt.append(dt)

 
       # the last rectangle
       h=y-self.top
       y=self.top
       rt=QCanvasRectangle(self.left,y,self.width,h,self.canvas)
       zlevel=zheight*(ci+1)+self.cview.min_brightness
       rt.setBrush(QBrush(self.cview.getColor(zlevel)))
       self.rts.append(rt)
       dt=QCanvasText(self.canvas)
       dt.setText(self.format%zlevel)
       dt.setFont(self.cview.axes.axfont)
       dt.move(self.left+self.width,y-3)
       dt.setZ(0)
       self.txt.append(dt)

       # last, text at zero level
       zlevel=self.cview.min_brightness
       dt=QCanvasText(self.canvas)
       dt.setText(self.format%zlevel)
       y=self.top+self.height
       dt.setFont(self.cview.axes.axfont)
       dt.move(self.left+self.width,y-3)
       dt.setZ(0)
       self.txt.append(dt)

    def updateFont(self,newfont):
     for i in range(len(self.txt)):
      self.txt[i].setFont(newfont)
