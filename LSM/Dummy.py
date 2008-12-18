
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
from LSM_GUI import *

######################################################
# Helper classes and methods for the GUI
#Singleton Class to make sure there is only one 
# instance of QApplication class
class Singleton:
  __single=None
  app=None # both these are class attributes, not instance attributes
  def __init__(self,args):
   if Singleton.__single:
     raise Singleton.__single
   Singleton.__single=self
   Singleton.app=QApplication(args)

def Handle(args,x=Singleton):
 try:
  single=x(args)
 except Singleton,s:
  single=s
 return single.app

class Dummy:
    def __init__(self,lsm_object,args):
     self.lsm=lsm_object
     self.myargs=args
     self.app=None
     self.win=None
    def display(self,**kw):
     if kw.has_key('app') and (kw['app']=='create'):
       self.app=Handle(self.myargs)
     else:
       self.app=qApp
       # check validity of qApp
       try:
         self.app.name()
       except:
         print "LSM: Unable to  create display, call with display(app='create')"
         return

     self.win=LSMWindow(self.lsm) 
     self.win.show()

     # stand alone mode
     if kw.has_key('app') and (kw['app']=='create'):
        self.app.connect(self.app,SIGNAL("lastWindowClosed()"),
          self.app, SLOT("quit()"))
        self.app.exec_loop()
     else: # running under browser
      self.app.connect(self.app,SIGNAL("lastWindowClosed()"),
        self.win, SLOT("close()"))


#######################################################
