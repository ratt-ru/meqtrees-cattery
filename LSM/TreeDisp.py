
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
import common_utils

from SDialog import *
from Timba.TDL import *
image0_data = \
    "\x89\x50\x4e\x47\x0d\x0a\x1a\x0a\x00\x00\x00\x0d" \
    "\x49\x48\x44\x52\x00\x00\x00\x10\x00\x00\x00\x10" \
    "\x08\x06\x00\x00\x00\x1f\xf3\xff\x61\x00\x00\x00" \
    "\x6b\x49\x44\x41\x54\x38\x8d\x63\x60\x18\x4c\xc0" \
    "\x99\x81\x81\x61\x3e\x03\x03\xc3\x7f\x02\xb8\x83" \
    "\x81\x81\xc1\x10\x9b\x01\x2b\x18\x18\x18\xf4\x88" \
    "\xb0\x48\x8f\x81\x81\xa1\x1f\x9b\xc4\x7f\x22\x34" \
    "\x33\xe1\x53\x4b\x8c\x01\x18\x6a\x99\xf0\xa9\x22" \
    "\x06\x8c\x1a\x40\x05\x03\x98\x91\xd8\x9a\x0c\x0c" \
    "\x0c\xb7\x18\x18\x18\x5e\x12\xd0\xa3\xc7\xc0\xc0" \
    "\xc0\xc8\xc0\xc0\xb0\x15\x5d\xc2\x81\x81\xb8\xa4" \
    "\x3c\x9f\x81\x81\xc1\x9d\x32\x77\x53\x13\x00\x00" \
    "\xab\xab\x1b\xc1\x78\x98\xf0\xf1\x00\x00\x00\x00" \
    "\x49\x45\x4e\x44\xae\x42\x60\x82"
image1_data = \
    "\x89\x50\x4e\x47\x0d\x0a\x1a\x0a\x00\x00\x00\x0d" \
    "\x49\x48\x44\x52\x00\x00\x00\x10\x00\x00\x00\x10" \
    "\x08\x06\x00\x00\x00\x1f\xf3\xff\x61\x00\x00\x02" \
    "\x32\x49\x44\x41\x54\x38\x8d\x95\x93\xcf\x4b\x94" \
    "\x71\x10\x87\x9f\x77\x7b\x93\x79\x65\xc5\x77\x29" \
    "\x62\x57\x48\x5c\x11\x64\xd3\x02\x03\x09\xeb\xa4" \
    "\x97\x40\xa8\x83\xe0\xc5\x6e\x11\x94\x46\xd7\x4e" \
    "\x41\x81\xf4\x1f\xa8\x20\xe2\xde\xdc\x08\x69\x2f" \
    "\x42\x87\xba\x54\x17\x0f\x6e\xee\x41\x4a\xa1\xec" \
    "\x3d\xb8\xec\x2e\xa9\xfb\xbe\xca\xae\x3b\xc8\xb7" \
    "\xde\x2e\xfd\x32\x12\xf5\xb9\xcc\x61\xe6\xf9\xcc" \
    "\x5c\xc6\xe2\x3f\xf4\x8c\x7e\x10\xb1\xe9\x55\xa3" \
    "\x1d\x82\x58\x8a\xae\x8b\x09\x72\x29\xe6\x6a\xe9" \
    "\x99\xf4\x81\x59\xeb\x1f\x31\x0a\x3c\xc2\x76\xee" \
    "\x61\x4b\x73\xf2\x2c\x88\xad\x78\x5b\x2e\x5a\x2d" \
    "\xed\x80\x3e\x17\x13\x3c\x4d\x56\x67\x37\x32\xcf" \
    "\x32\xe1\x81\x80\xbe\x07\x5f\xe2\x6a\xea\x6f\x12" \
    "\xf1\x44\xe7\xed\x7e\x61\xa0\xdb\x41\x04\x1a\x1b" \
    "\x2d\xab\xb2\x1b\x86\xab\x9f\xeb\x4c\xbe\x0e\xf0" \
    "\x0a\xfe\x16\xea\x0d\x27\x83\xcc\xdb\xcc\x7c\x26" \
    "\xb4\xfe\xda\x9c\xeb\xeb\x4e\x76\xde\xbf\x01\x5d" \
    "\xad\x8d\x07\x2e\xfb\x45\xf1\x6b\x18\x4e\x2e\x94" \
    "\x78\x99\xf3\xb7\x45\xf3\x97\x12\xd5\x6c\x31\xf2" \
    "\xb3\xf7\x24\x11\x4f\x74\xde\xea\x3f\x5c\x06\x68" \
    "\x39\x67\x59\x23\xfd\x09\x52\x6d\xb1\x33\x48\xea" \
    "\xb1\xeb\xba\x0d\x91\x9e\xd1\xf7\x82\xed\xdc\x19" \
    "\xec\x81\x6b\x17\x0e\x97\x7f\xd1\xd5\x6e\x59\x03" \
    "\x17\x05\x45\x46\xca\xe6\x72\x8b\x2d\x48\x2f\x22" \
    "\xcd\xf1\xa8\x1c\xe5\xfe\x66\x76\x29\x86\x44\x89" \
    "\xfa\xa6\xed\x6a\x44\xd1\x8e\x58\x14\xbc\x20\x38" \
    "\x76\x00\x80\x02\x62\x3b\x1d\xb6\x20\x56\xdd\x08" \
    "\x65\xff\x44\x3e\x62\x0b\x6a\x88\x44\x14\x5d\x57" \
    "\x85\xe5\x72\xec\x04\x32\xa8\x06\x08\xf5\x42\x44" \
    "\x4c\x90\x53\x2d\xed\x38\xf6\xf1\xb7\x3b\x36\x88" \
    "\xd1\x5a\xcc\xe4\xf3\x91\x14\x73\x35\x8c\xa6\xd1" \
    "\x12\x63\xd3\x7b\xe1\x51\xf2\xd8\x74\x25\x44\x7d" \
    "\x50\x6f\x21\x6e\x97\x8a\x91\xf4\x4c\x1a\x31\xc1" \
    "\xb8\xbf\xe5\xaf\x63\x02\xc6\xe7\x2a\x87\x86\x8c" \
    "\xcf\x55\x42\xc7\x28\xa5\xb2\x57\x48\xb0\x3c\xa1" \
    "\xe8\xb6\x0d\xd0\xc3\x0b\x3f\xaf\x83\xd7\xf3\x6b" \
    "\xbc\xd2\x36\xa7\x7d\x22\x5b\x09\xeb\x80\x03\xfc" \
    "\xa9\x8a\xb7\x11\xb0\xb8\xe6\x17\x92\x66\xf1\xae" \
    "\xcb\xea\x8a\x20\xfb\xa7\x00\x96\x96\x97\xb8\x79" \
    "\xc5\xf5\x83\xbd\xa6\xf9\xf2\xae\x34\x7d\x2a\x6a" \
    "\xa7\xd9\xff\xd6\xa0\xfb\x4a\x55\xc1\x2b\x6e\xb2" \
    "\xf4\x51\x6b\x85\xc2\x56\xf6\xfc\xf7\x77\x0f\x5d" \
    "\x56\x72\x62\xcb\xee\xd4\xf4\xd4\xc1\x6f\x1c\x1a" \
    "\x1e\xc2\x8d\xba\x0d\x65\x93\x6a\xf5\x49\x0c\x88" \
    "\xb8\x6d\x0a\xfb\x02\x85\x98\xc9\xe7\x5d\xbc\xcd" \
    "\xc0\x0f\xb6\x39\x4d\x2d\x3b\x9f\x05\xe0\x07\xf4" \
    "\x5f\xfe\x71\x73\xe0\xf9\x0d\x00\x00\x00\x00\x49" \
    "\x45\x4e\x44\xae\x42\x60\x82"

class TreeDisp(QDialog):
    def __init__(self,parent = None,name = None,modal = 0,fl = 0,\
            root=None):
        QDialog.__init__(self,parent,name,modal,fl)

        self.image0 = QPixmap()
        self.image0.loadFromData(image0_data,"PNG")
        self.image1 = QPixmap()
        self.image1.loadFromData(image1_data,"PNG")
        if not name:
            self.setName("TreeDisp")

        self.root=root

        FormLayout = QVBoxLayout(self,11,6,"FormLayout")

        self.lv= QListView(self,"listView")
        self.lv.addColumn(self.__tr("Tree"))
        self.lv.addColumn(self.__tr("Node Type"))
        FormLayout.addWidget(self.lv)

        self.pushButton = QPushButton(self,"pushButton")
        FormLayout.addWidget(self.pushButton)
        self.connect(self.pushButton,SIGNAL("clicked()"),self.accept)

        self.languageChange()

        #self.resize(QSize(600,480).expandedTo(self.minimumSizeHint()))
        self.clearWState(Qt.WState_Polished)
        # handle clicks
        self.connect(self.lv,SIGNAL("clicked(QListViewItem *)"),
             self.displayInitrec)


    def languageChange(self):
        self.setCaption(self.__tr("TreeDisp"))
        self.lv.header().setLabel(0,self.__tr("node"))
        self.lv.header().setLabel(1,self.__tr("class"))
        self.lv.clear()
        #item=QListViewItem(self.lv,None)
        item=NodeItem(self.lv,None,self.root)
        item.setText(0,self.__tr(self.root.name))
        item.setText(1,self.__tr(self.root.classname))
        item.setPixmap(0,self.image0)
        item.setOpen(1)
        self.plotTree(self.root,item)
        self.pushButton.setText(self.__tr("OK"))


    def setTitle(self,text):
        self.setCaption(self.__tr(text))

    def plotTree(self,node,parent):
     # 'node' is a child of 'parent'
     # plot node and its children, and call this recursively
     # first plot root
     # get children
     clist=node.children
     #print clist
     if clist!=None:
      for ci,nstub in clist:
       #citem=QListViewItem(parent,None)
       citem=NodeItem(parent,None,nstub)
       citem.setText(0,self.__tr(str(ci)+":"+nstub.name))
       citem.setText(1,self.__tr(nstub.classname))
       #citem.setPixmap(0,self.image1)
       self.plotTree(nstub,citem)
      
    def __tr(self,s,c = None):
        return qApp.translate("TreeDisp",s,c)

    
    def displayInitrec(self,item):
       if item !=None:
         dialog=SDialog(self,item.node.name,0)
         print common_utils.get_stokes_I_parms(item.node)
         tmp_str="<u>Init Record:</u><br/><ul>"
         irec=item.node.initrec()
         for kk in irec.keys():
          # do another level of indirection if it is a dict
          if isinstance(irec[kk],dict):
            kstr="<li><font color=\"blue\">"+str(kk)+"::</font><ul>"
            subrec=irec[kk]
            for bb in subrec.keys():
             kstr=kstr+"<li><font color=\"blue\">"+str(bb)+"::</font>"+str(subrec[bb])+"</li>"
            kstr=kstr+"</ul></li>"
          else:
             kstr="<li><font color=\"blue\">"+str(kk)+"::</font>"+str(irec[kk])+"</li>"
          tmp_str=tmp_str+kstr
         tmp_str=tmp_str+"</ul>"
         dialog.setInfoText(tmp_str)
         dialog.setTitle(item.node.name)
         dialog.resize(QSize(450,300))
         dialog.show()


class NodeItem(QListViewItem):
   def __init__(self,parent=None,after=None,node=None):
      QListViewItem.__init__(self,parent,after)
      self.node=node

def main(args):
  ns=NodeScope()
  x=ns.x<<Meq.Parm(0)
  y=ns.y<<Meq.Parm(1)
  rt=ns<<Meq.Sin(x)+Meq.Cos(y+5*x)
  ns.Resolve()
  print ns.AllNodes()
  print ns.RootNodes()
  app=QApplication(args)
  win=TreeDisp(None,None,0,0,rt)
  win.setTitle("Testing...")
  win.show()
  app.connect(app,SIGNAL("lastWindowClosed()"),
               app,SLOT("quit()"))
  app.exec_loop()

if __name__=="__main__":
   main(sys.argv)

