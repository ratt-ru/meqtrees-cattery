#!/usr/bin/python

#############################################
#### GL window to display a patch (l,m) plane
#### for given freq,time in 3D
#############################################

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

import math
from qt import *
try:
 from qtgl import *
except:
 pass
#from OpenGL.GL import *
#from OpenGL.GLU import *
import sys
import Timba.array

#############################################################
class PatchGL(QGLWidget):
 def __init__(self,parent=None,name=None,my_arr=None,my_lims=None):
  QGLWidget.__init__(self,parent,name)
  self.arr=my_arr
  self.lims=my_lims
  # normalize such that max value is below ~1
  self.arr=(self.arr-self.arr.min())
  self.arr=self.arr/self.arr.max()

  print self.arr.min(),self.arr.max()

  # create a GLList
  self.surf=None


 # painting routine
 def paintGL(self):
   glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
   # call display list
   glCallList(self.surf) 
   glFlush()
 

 # create the plot, used to generate display list
 def plot_surface(self):
   glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
   # the drawing area is limited to [-1,1] [-1,1] in x and y
   # so find a suitable x,y scale or cell size
   # each cell size
   dx=2.0/self.lims[0]
   dy=2.0/self.lims[1]
   xoff=(self.lims[0]-1)/2.0*dx
   yoff=(self.lims[1]-1)/2.0*dy
   contour_levels=30
   glBegin(GL_LINES)
   for ci in range(self.lims[0]-1):
    for cj in range(self.lims[1]-1):
    # glVertex3f(ci*dx-xoff,cj*dy-yoff,self.arr[ci][cj])
    # glVertex3f((ci+1)*dx-xoff,cj*dy-yoff,self.arr[ci+1][cj])
    # glVertex3f((ci+1)*dx-xoff,(cj+1)*dy-yoff,self.arr[ci+1][cj+1])
    # glVertex3f(ci*dx-xoff,(cj+1)*dy-yoff,self.arr[ci][cj+1])

     hh=_NodeNA( ci*dx-xoff,cj*dy-yoff,self.arr[ci][cj],\
              (ci+1)*dx-xoff,cj*dy-yoff,self.arr[ci+1][cj],\
              (ci+1)*dx-xoff,(cj+1)*dy-yoff,self.arr[ci+1][cj+1],\
              ci*dx-xoff,(cj+1)*dy-yoff,self.arr[ci][cj+1])
     # set up some contours
     
     for ck in range(contour_levels):
      clevel=self.arr.min()+ck*self.arr.max()/contour_levels
      # get color
      color=self.getColor(clevel)
      glColor3f(color[0], color[1], color[2])
      self.plot_contour(hh,clevel)
      #print "Point %f,%f"%(ci*dx-xoff,cj*dy-yoff)
   glEnd()
   



 # initialization
 def initializeGL(self):
  glClearColor(0.0,0.0,0.0,0.0)
  # create a GL Display List
  self.surf=glGenLists(1)
  glNewList(self.surf,GL_COMPILE)
  self.plot_surface()
  glEndList()




 # scale with resize
 def resizeGL(self,width,height):
    glViewport( 0, 0, width, height )
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    gluOrtho2D(-1,1,-1,1)
    glMatrixMode(GL_MODELVIEW)


 def mouseMoveEvent(self,e):
  # refresh
  #self.repaint()
  pass

 # return a colour according to z, in OpenGL format
 def getColor(self,z):
   if self.arr.max()==0:
     return (0,0,0) # no color : when Q,U,V is zero

   cl=z/self.arr.max()# normalized in [0,1] 
   if cl < 0.25:
     return (0,cl*4.0,1.0)
   elif cl < 0.5:
     return (0,1.0,2.0-cl*4.0)
   elif cl < 0.75:
     return (4.0*cl-2.0,1.0,0)
   else:
     return (1.0,4.0-4.0*cl,0)


 # hh: _NodeNA object with the (x,y,z) data of 
 # points 0,1,2,3
 # points are taken counterclockwise
 # phi: current potential of contour 
 def plot_contour(self,hh,phi):
  E=0
  H=0
  L=0
  nl=0  # no. of lines to plot
  for i in range(len(hh.h)):
    if (hh.h[i].z==phi):
     E+=1
    elif (hh.h[i].z<phi):
     L+=1
    else:
     H+=1
  
  # now determine no of lines in each rectangle
  ############# case 1
  if ((H==4) or (L==4) or ((E==1)and (H==3))\
   or ((E==1)and (L==3))):
    # no lines to plot
    return 0
  ############# case 2
  if E==4:
   nl=4
   # plot
   for i in range(len(hh.h)-1):
    glVertex3f(hh.h[i].x,hh.h[i].y,hh.h[i].z)
    glVertex3f(hh.h[i+1].x,hh.h[i+1].y,hh.h[i+1].z)
   glVertex3f(hh.h[3].x,hh.h[3].y,hh.h[3].z)
   glVertex3f(hh.h[0].x,hh.h[0].y,hh.h[0].z)
   return nl
  
  ############# case 3
  if E==3:
   i=0
   while(hh.h[i].z == phi):
    i+=1 # find the one not equal
   nl=2
   glVertex3f(hh.h[i].next.x,hh.h[i].next.y,hh.h[i].next.z)
   glVertex3f(hh.h[i].opp.x,hh.h[i].opp.y,hh.h[i].opp.z)
   glVertex3f(hh.h[i].prev.x,hh.h[i].prev.y,hh.h[i].prev.z)
   glVertex3f(hh.h[i].opp.x,hh.h[i].opp.y,hh.h[i].opp.z)
   return nl

  ############# case 4
  if E==2:
   i=0
   while(hh.h[i].z != phi):
    i+=1
   if(hh.h[i].z==hh.h[i].opp.z):
     if L>0:
      nl=1
      glVertex3f(hh.h[i].x,hh.h[i].y,hh.h[i].z)
      glVertex3f(hh.h[i].opp.x,hh.h[i].opp.y,hh.h[i].opp.z)
     else:
      return 0
   else:
    if((H==1) and(L==1)) :
      high=0
      low=0
      while(hh.h[high].z<=phi):
        high+=1
      while(hh.h[low].z>=phi):
        low+=1
      # now interpolate to find the next point
      tempx1=(phi-hh.h[low].z)*(hh.h[high].x-hh.h[low].x)/(hh.h[high].z-hh.h[low].z)
      tempy1=(phi-hh.h[low].z)*(hh.h[high].y-hh.h[low].y)/(hh.h[high].z-hh.h[low].z)
      nl=2
      if (hh.h[i].z ==hh.h[i].next.z):
       glVertex3f(hh.h[i].x,hh.h[i].y,phi)
       glVertex3f(tempx1,tempy1,phi)
       glVertex3f(hh.h[i].next.x,hh.h[i].next.y,phi)
       glVertex3f(tempx1,tempy1,phi)
      else:
       glVertex3f(hh.h[i].x,hh.h[i].y,phi)
       glVertex3f(tempx1,tempy1,phi)
       glVertex3f(hh.h[i].prev.x,hh.h[i].prev.y,phi)
       glVertex3f(tempx1,tempy1,phi)
  
      return nl
    else:
     # only one line
     nl=1
     if (hh.h[i].z ==hh.h[i].next.z):
       glVertex3f(hh.h[i].next.x,hh.h[i].next.y,phi)
       glVertex3f(hh.h[i].x,hh.h[i].y,phi)
     else:
       glVertex3f(hh.h[i].prev.x,hh.h[i].prev.y,phi)
       glVertex3f(hh.h[i].x,hh.h[i].y,phi)
     return nl
   # will not get here
   return 0
  ############# case 5
  if ((H==2) and (L==2)):
   i=0
   while(hh.h[i].z <phi):
    i+=1
   if (hh.h[i].opp.z > phi):
    nl=2
    tempx1=(hh.h[i].x-hh.h[i].next.x)*(phi-hh.h[i].next.z)/(hh.h[i].z-hh.h[i].next.z)+hh.h[i].next.x
    tempy1=(hh.h[i].y-hh.h[i].next.y)*(phi-hh.h[i].next.z)/(hh.h[i].z-hh.h[i].next.z)+hh.h[i].next.y
    tempx2=(hh.h[i].x-hh.h[i].prev.x)*(phi-hh.h[i].prev.z)/(hh.h[i].z-hh.h[i].prev.z)+hh.h[i].prev.x
    tempy2=(hh.h[i].y-hh.h[i].prev.y)*(phi-hh.h[i].prev.z)/(hh.h[i].z-hh.h[i].prev.z)+hh.h[i].prev.y
    glVertex3f(tempx1,tempy1,phi)
    glVertex3f(tempx2,tempy2,phi)
    # to this again for the other line
    tempx2=(hh.h[i].opp.x-hh.h[i].opp.next.x)*(phi-hh.h[i].opp.next.z)/(hh.h[i].opp.z-hh.h[i].opp.next.z)+hh.h[i].opp.next.x
    tempy2=(hh.h[i].opp.y-hh.h[i].opp.next.y)*(phi-hh.h[i].opp.next.z)/(hh.h[i].opp.z-hh.h[i].opp.next.z)+hh.h[i].opp.next.y
    tempx1=(hh.h[i].opp.x-hh.h[i].opp.prev.x)*(phi-hh.h[i].opp.prev.z)/(hh.h[i].opp.z-hh.h[i].opp.prev.z)+hh.h[i].opp.prev.x
    tempy1=(hh.h[i].opp.y-hh.h[i].opp.prev.y)*(phi-hh.h[i].opp.prev.z)/(hh.h[i].opp.z-hh.h[i].opp.prev.z)+hh.h[i].opp.prev.y
    glVertex3f(tempx1,tempy1,phi)
    glVertex3f(tempx2,tempy2,phi)
    return nl
   else:
    nl=1
    if (hh.h[i].next.z >phi):
     tempx2=(hh.h[i].x-hh.h[i].prev.x)*(phi-hh.h[i].prev.z)/(hh.h[i].z-hh.h[i].prev.z)+hh.h[i].prev.x
     tempy2=(hh.h[i].y-hh.h[i].prev.y)*(phi-hh.h[i].prev.z)/(hh.h[i].z-hh.h[i].prev.z)+hh.h[i].prev.y
     tempx1=(hh.h[i].opp.x-hh.h[i].next.x)*(phi-hh.h[i].next.z)/(hh.h[i].opp.z-hh.h[i].next.z)+hh.h[i].next.x
     tempy1=(hh.h[i].opp.y-hh.h[i].next.y)*(phi-hh.h[i].next.z)/(hh.h[i].opp.z-hh.h[i].next.z)+hh.h[i].next.y
    else:
     tempx1=(hh.h[i].x-hh.h[i].next.x)*(phi-hh.h[i].next.z)/(hh.h[i].z-hh.h[i].next.z)+hh.h[i].next.x
     tempy1=(hh.h[i].y-hh.h[i].next.y)*(phi-hh.h[i].next.z)/(hh.h[i].z-hh.h[i].next.z)+hh.h[i].next.y
     tempx2=(hh.h[i].opp.x-hh.h[i].prev.x)*(phi-hh.h[i].prev.z)/(hh.h[i].opp.z-hh.h[i].prev.z)+hh.h[i].prev.x
     tempy2=(hh.h[i].opp.y-hh.h[i].prev.y)*(phi-hh.h[i].prev.z)/(hh.h[i].opp.z-hh.h[i].prev.z)+hh.h[i].prev.y

    glVertex3f(tempx1,tempy1,phi)
    glVertex3f(tempx2,tempy2,phi)
    return nl
  ############# case 6
  if E==1:
   i=0
   while (hh.h[i].z != phi):
    i+=1
   if (H==1) :
    high=0
    while(hh.h[high].z <= phi) :
     high+=1
    if (hh.h[high].z ==phi) :
     # 2 lines
     nl=2
     tempx1=(hh.h[high].x-hh.h[high].next.x)*(phi-hh.h[high].next.z)/(hh.h[high].z-hh.h[high].next.z)+hh.h[high].next.x
     tempy1=(hh.h[high].y-hh.h[high].next.y)*(phi-hh.h[high].next.z)/(hh.h[high].z-hh.h[high].next.z)+hh.h[high].next.y
     tempx2=(hh.h[high].x-hh.h[high].prev.x)*(phi-hh.h[high].prev.z)/(hh.h[high].z-hh.h[high].prev.z)+hh.h[high].prev.x
     tempy2=(hh.h[high].y-hh.h[high].prev.y)*(phi-hh.h[high].prev.z)/(hh.h[high].z-hh.h[high].prev.z)+hh.h[high].prev.y
     glVertex3f(tempx1,tempy1,phi)
     glVertex3f(hh.h[i].x,hh.h[i].y,phi)
     glVertex3f(tempx2,tempy2,phi)
     glVertex3f(hh.h[i].x,hh.h[i].y,phi)
     return nl
    else:
     # 1 line only
     nl=1
     if (hh.h[high].next.z==phi):
      tempx1=(hh.h[high].x-hh.h[high].prev.x)*(phi-hh.h[high].prev.z)/(hh.h[high].z-hh.h[high].prev.z)+hh.h[high].prev.x
      tempy1=(hh.h[high].y-hh.h[high].prev.y)*(phi-hh.h[high].prev.z)/(hh.h[high].z-hh.h[high].prev.z)+hh.h[high].prev.y
     else:
      tempx1=(hh.h[high].x-hh.h[high].next.x)*(phi-hh.h[high].next.z)/(hh.h[high].z-hh.h[high].next.z)+hh.h[high].next.x
      tempy1=(hh.h[high].y-hh.h[high].next.y)*(phi-hh.h[high].next.z)/(hh.h[high].z-hh.h[high].next.z)+hh.h[high].next.y
     glVertex3f(tempx1,tempy1,phi)
     glVertex3f(hh.h[i].x,hh.h[i].y,phi)
     return nl
   else: # L==1
    low=0
    while (hh.h[low].z >=phi):
     low+=1
    if(hh.h[low].opp.z==phi): 
     nl=2
     tempx1=(hh.h[low].x-hh.h[low].next.x)*(phi-hh.h[low].next.z)/(hh.h[low].z-hh.h[low].next.z)+hh.h[low].next.x
     tempy1=(hh.h[low].y-hh.h[low].next.y)*(phi-hh.h[low].next.z)/(hh.h[low].z-hh.h[low].next.z)+hh.h[low].next.y
     tempx2=(hh.h[low].x-hh.h[low].prev.x)*(phi-hh.h[low].prev.z)/(hh.h[low].z-hh.h[low].prev.z)+hh.h[low].prev.x
     tempy2=(hh.h[low].y-hh.h[low].prev.y)*(phi-hh.h[low].prev.z)/(hh.h[low].z-hh.h[low].prev.z)+hh.h[low].prev.y
     glVertex3f(tempx1,tempy1,phi)
     glVertex3f(hh.h[i].x,hh.h[i].y,phi)
     glVertex3f(tempx2,tempy2,phi)
     glVertex3f(hh.h[i].x,hh.h[i].y,phi)
     return nl
    else:
     nl=1
     if(hh.h[low].next.z==phi):
      tempx1=(hh.h[low].x-hh.h[low].prev.x)*(phi-hh.h[low].prev.z)/(hh.h[low].z-hh.h[low].prev.z)+hh.h[low].prev.x
      tempy1=(hh.h[low].y-hh.h[low].prev.y)*(phi-hh.h[low].prev.z)/(hh.h[low].z-hh.h[low].prev.z)+hh.h[low].prev.y
     else:
      tempx1=(hh.h[low].x-hh.h[low].next.x)*(phi-hh.h[low].next.z)/(hh.h[low].z-hh.h[low].next.z)+hh.h[low].next.x
      tempy1=(hh.h[low].y-hh.h[low].next.y)*(phi-hh.h[low].next.z)/(hh.h[low].z-hh.h[low].next.z)+hh.h[low].next.y

     glVertex3f(tempx1,tempy1,phi)
     glVertex3f(hh.h[i].x,hh.h[i].y,phi)
     return nl
  ############# case 7
  if ((H==1)and(L==3)):
   nl=1
   high=0
   while(hh.h[high].z<=phi):
    high+=1
   tempx1=(hh.h[high].x-hh.h[high].next.x)*(phi-hh.h[high].next.z)/(hh.h[high].z-hh.h[high].next.z)+hh.h[high].next.x
   tempy1=(hh.h[high].y-hh.h[high].next.y)*(phi-hh.h[high].next.z)/(hh.h[high].z-hh.h[high].next.z)+hh.h[high].next.y
   tempx2=(hh.h[high].x-hh.h[high].prev.x)*(phi-hh.h[high].prev.z)/(hh.h[high].z-hh.h[high].prev.z)+hh.h[high].prev.x
   tempy2=(hh.h[high].y-hh.h[high].prev.y)*(phi-hh.h[high].prev.z)/(hh.h[high].z-hh.h[high].prev.z)+hh.h[high].prev.y
   glVertex3f(tempx1,tempy1,phi)
   glVertex3f(tempx2,tempy2,phi)
   return nl
  ############# case 8
  if ((L==1)and(H==3)):
   nl=1
   low=0
   while(hh.h[low].z>=phi):
    low+=1
   tempx1=(hh.h[low].x-hh.h[low].prev.x)*(phi-hh.h[low].prev.z)/(hh.h[low].z-hh.h[low].prev.z)+hh.h[low].prev.x
   tempy1=(hh.h[low].y-hh.h[low].prev.y)*(phi-hh.h[low].prev.z)/(hh.h[low].z-hh.h[low].prev.z)+hh.h[low].prev.y
   tempx2=(hh.h[low].x-hh.h[low].next.x)*(phi-hh.h[low].next.z)/(hh.h[low].z-hh.h[low].next.z)+hh.h[low].next.x
   tempy2=(hh.h[low].y-hh.h[low].next.y)*(phi-hh.h[low].next.z)/(hh.h[low].z-hh.h[low].next.z)+hh.h[low].next.y
   glVertex3f(tempx1,tempy1,phi)
   glVertex3f(tempx2,tempy2,phi)
   return nl
  # will not get here
  return nl
####################################################################
# helper class for contour plotting
class _NodeN:
  def __init__(self,x,y,z):
   self.x=x
   self.y=y
   self.z=z
   self.next=None
   self.prev=None
   self.opp=None
class _NodeNA:
   def __init__(self,x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3):
    self.h=[]
    self.h.append(_NodeN(x0,y0,z0))
    self.h.append(_NodeN(x1,y1,z1))
    self.h.append(_NodeN(x2,y2,z2))
    self.h.append(_NodeN(x3,y3,z3))
    self.h[0].next=self.h[1]
    self.h[0].prev=self.h[3]
    self.h[0].opp=self.h[2]
    self.h[1].next=self.h[0]
    self.h[1].prev=self.h[2]
    self.h[1].opp=self.h[3]
    self.h[2].next=self.h[3]
    self.h[2].prev=self.h[1]
    self.h[2].opp=self.h[0]
    self.h[3].next=self.h[0]
    self.h[3].prev=self.h[2]
    self.h[3].opp=self.h[1]



####################################################################
class PatchGLDialog(QDialog):
    def __init__(self,parent = None,name = None,modal = 0,fl = 0,\
          my_arr=None,my_lims=None):
        QDialog.__init__(self,parent,name,modal,fl)
        if not name:
            self.setName("PatchGLDialog")


        FormLayout = QHBoxLayout(self,11,6,"FormLayout")

        # create GL window
        self.glwidget=PatchGL(self,"Patch",my_arr,my_lims)
        FormLayout.addWidget(self.glwidget)
        self.languageChange()

        self.resize(QSize(400,400).expandedTo(self.minimumSizeHint()))
        self.clearWState(Qt.WState_Polished)


    def languageChange(self):
        self.setCaption(self.__tr("PatchGLDialog"))


    def __tr(self,s,c = None):
        return qApp.translate("PatchGLDialog",s,c)

########################################################################
def main(args):
  app=QApplication(args)
  my_arr=Timba.array.fromfunction(lambda x,y: (x**2+y**2)/100.0,(12,14)) 
  my_lims=my_arr.shape
  win=PatchGLDialog(None,"Patch",1,0,my_arr,my_lims)
  win.show()
  app.connect(app,SIGNAL("lastWindowClosed()"),
               app,SLOT("quit()"))
  app.exec_loop()

if __name__=="__main__":
   main(sys.argv)
