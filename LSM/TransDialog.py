#!/usr/bin/python


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

image0_data = [
"224 36 2 1",
"# c #000000",
". c #ffffff",
"....####.........................................................####................................................................####.....................................................####..............................",
"....#...............................................................#................................................................#...........................................................#..............................",
"....#...............................................................#................................................................#...........................................................#..............................",
"....#...........................#...................................#................................................................#...........................#...............................#..............................",
"....#...........#######.........#...................................#................................................................#............#######........#...............................#..............................",
"....#.............#....#.......###..................................#................................................................#..............#....#......###..............................#..............................",
"....#.............#....#.......###..................................#................................................................#..............#....#......###..............................#..............................",
"....#.............#....#......#..#..................................#................................................................#..............#....#.....#..#..............................#..............................",
"....#.............#####.......#..##.................................#................................................................#..............#####......#..##.........#....##.............#..............................",
"....#.............#...#.......#####.................................#.........................##..........................###........#..............#...#......#####.........#.....#.............#..............................",
"....#.............#....#......#...#.................................#.........................#.............................#........#..............#....#.....#...#.........#.....#.............#..............................",
"....#.............#....#.....#....##..####..###.##.##.##............#.........................#.............................#........#..............#....#....#....##..###...#..####.............#..............................",
"....#...........#####..###..###..####.#..#.#####.#.##.#.............#.........................#..............#..............#........#............#####..###.###..#####...#..#.#...#.............#.....................##.......",
"....#.................................#..#.#......###.#.............#.........................#.............###.............#........#................................#...#..#.#...#.............#...........#.........##.......",
"....#.................................#####.####..#..#..............#.........................#.............###.............#........#.................................###...##.####.............#...........#.........##.......",
"....#...............................................................#.........................#............#.###............#........#...........................................................#...........#.........######...",
"....#...............................................................#........#########........#............#.###............#........#...........................................................#...........#.........##...##..",
"....#...............................................................#.........................#...........#...##............#........#...........................................................#.......#########.....##...##..",
"....#...............................................................#........#########........#...........#######...........#........#...........................................................#...........#.........##...##..",
"....#...............................................................#.........................#...........#....##...........#........#...........................................................#...........#.........##...##..",
"....#..........########.............................................#.........................#.........###...#####.........#........#...........########........................................#...........#.........#.####...",
"....#............#....##............................................#.........................#.............................#........#.............#....##.......................................#...........#..................",
"....#............#......#...........................................#.........................#.............................#........#.............#......#......................................#..............................",
"....#............#......#..####...###...............................#.........................#.............................#........#.............#......#.####...###...........................#..............................",
"....#............#......#.#....#.##.#...............................#.........................#.............................#........#.............#......##....#.##.#........#....##............#..............................",
"....#............#......#.######.#..................................#.........................##..........................###........#.............#......#######.#...........#.....#............#..............................",
"....#............#......#.#......#..................................#................................................................#.............#......##......#...........#.....#............#..............................",
"....#............#....##..##...#.##..#.####..###.##.##.##...........#................................................................#.............#....##.##...#.##..#.###...#..####............#..............................",
"....#..........########.....###...###..#..#.#####.#.##.#............#................................................................#...........########....###...###.#...#..#.#...#............#..............................",
"....#..................................#..#.#......###.#............#................................................................#.................................#...#..#.#...#............#..............................",
"....#..................................#####.####..#..#.............#................................................................#..................................###...##.####............#..............................",
"....#...............................................................#................................................................#...........................................................#..............................",
"....#...............................................................#................................................................#...........................................................#..............................",
"....####.........................................................####................................................................####.....................................................####..............................",
"................................................................................................................................................................................................................................",
"................................................................................................................................................................................................................................"
]

class TransDialog(QDialog):
    def __init__(self,parent = None,name = None,modal = 0,fl = 0, cview=None, lsm=None):
        QDialog.__init__(self,parent,name,modal,fl)

        self.cview=cview
        self.lsm=lsm
        self.A=[[1,0],[0,1]]
        self.b=[0,0]
        self.image0 = QPixmap(image0_data)

        if not name:
            self.setName("TransDialog")

        self.setSizeGripEnabled(1)

        MyDialogLayout = QVBoxLayout(self,0,-1,"MyDialogLayout")

        layout0 = QHBoxLayout(None,0,6,"layout0")
        spacer0 = QSpacerItem(71,41,QSizePolicy.Expanding,QSizePolicy.Minimum)
        spacer1 = QSpacerItem(71,41,QSizePolicy.Expanding,QSizePolicy.Minimum)
        self.pixmapLabel1 = QLabel(self,"pixmapLabel1")
        self.pixmapLabel1.setPixmap(self.image0)
        self.pixmapLabel1.setScaledContents(0)
        layout0.addItem(spacer0)
        layout0.addWidget(self.pixmapLabel1)
        layout0.addItem(spacer1)
        MyDialogLayout.addLayout(layout0)

        layout11 = QHBoxLayout(None,0,6,"layout11")

        self.groupBox1 = QGroupBox(self,"groupBox1")
        self.groupBox1.setColumnLayout(0,Qt.Vertical)
        self.groupBox1.layout().setSpacing(6)
        self.groupBox1.layout().setMargin(11)
        groupBox1Layout = QVBoxLayout(self.groupBox1.layout())
        groupBox1Layout.setAlignment(Qt.AlignCenter)

        layout7 = QHBoxLayout(None,0,6,"layout7")

        self.textEdit3 = QLineEdit(self.groupBox1,"textEdit3")
        layout7.addWidget(self.textEdit3)

        self.textEdit4 = QLineEdit(self.groupBox1,"textEdit4")
        layout7.addWidget(self.textEdit4)
        groupBox1Layout.addLayout(layout7)

        layout8 = QHBoxLayout(None,0,6,"layout8")

        self.textEdit5 = QLineEdit(self.groupBox1,"textEdit5")
        layout8.addWidget(self.textEdit5)

        self.textEdit6 = QLineEdit(self.groupBox1,"textEdit6")
        layout8.addWidget(self.textEdit6)
        groupBox1Layout.addLayout(layout8)
        layout11.addWidget(self.groupBox1)

        self.groupBox2 = QGroupBox(self,"groupBox2")
        self.groupBox2.setColumnLayout(0,Qt.Vertical)
        self.groupBox2.layout().setSpacing(6)
        self.groupBox2.layout().setMargin(11)
        groupBox2Layout = QVBoxLayout(self.groupBox2.layout())
        groupBox2Layout.setAlignment(Qt.AlignCenter)

        layout10 = QVBoxLayout(None,0,6,"layout10")

        self.textEdit7 = QLineEdit(self.groupBox2,"textEdit7")
        layout10.addWidget(self.textEdit7)

        self.textEdit8 = QLineEdit(self.groupBox2,"textEdit8")
        layout10.addWidget(self.textEdit8)
        groupBox2Layout.addLayout(layout10)
        layout11.addWidget(self.groupBox2)
        MyDialogLayout.addLayout(layout11)

        Layout1 = QHBoxLayout(None,0,6,"Layout1")

        self.buttonHelp = QPushButton(self,"buttonHelp")
        self.buttonHelp.setAutoDefault(1)
        Layout1.addWidget(self.buttonHelp)
        Horizontal_Spacing2 = QSpacerItem(20,20,QSizePolicy.Expanding,QSizePolicy.Minimum)
        Layout1.addItem(Horizontal_Spacing2)

        self.buttonOk = QPushButton(self,"buttonOk")
        self.buttonOk.setAutoDefault(1)
        self.buttonOk.setDefault(1)
        Layout1.addWidget(self.buttonOk)

        self.buttonCancel = QPushButton(self,"buttonCancel")
        self.buttonCancel.setAutoDefault(1)
        Layout1.addWidget(self.buttonCancel)
        MyDialogLayout.addLayout(Layout1)

        self.languageChange()

        #self.resize(QSize(709,514).expandedTo(self.minimumSizeHint()))
        self.clearWState(Qt.WState_Polished)


        # set default values
        self.textEdit3.setText(str(1.0))
        self.connect(self.textEdit3,SIGNAL("textChanged(const QString &)"),self.lineEdita0)

        self.textEdit4.setText(str(0.0))
        self.connect(self.textEdit4,SIGNAL("textChanged(const QString &)"),self.lineEdita1)

        self.textEdit5.setText(str(0.0))
        self.connect(self.textEdit5,SIGNAL("textChanged(const QString &)"),self.lineEdita2)

        self.textEdit6.setText(str(1.0))
        self.connect(self.textEdit6,SIGNAL("textChanged(const QString &)"),self.lineEdita3)

        self.textEdit7.setText(str(0.0))
        self.connect(self.textEdit7,SIGNAL("textChanged(const QString &)"),self.lineEditb1)
        self.textEdit8.setText(str(0.0))
        self.connect(self.textEdit8,SIGNAL("textChanged(const QString &)"),self.lineEditb2)



        self.connect(self.buttonOk,SIGNAL("clicked()"),self.accept)
        self.connect(self.buttonCancel,SIGNAL("clicked()"),self.reject)


    def languageChange(self):
        self.setCaption(self.__tr("Linear Transformation"))
        self.groupBox1.setTitle(self.__tr("A"))
        self.groupBox2.setTitle(self.__tr("b"))
        self.buttonHelp.setText(self.__tr("&Help"))
        self.buttonHelp.setAccel(self.__tr("F1"))
        self.buttonOk.setText(self.__tr("&OK"))
        self.buttonOk.setAccel(QString.null)
        self.buttonCancel.setText(self.__tr("&Cancel"))
        self.buttonCancel.setAccel(QString.null)


    def __tr(self,s,c = None):
        return qApp.translate("TransDialog",s,c)

    def lineEdita0(self,newText):
     if newText.length()>0:
      try:
       self.A[0][0]=float(newText.ascii())
      except (TypeError,ValueError):
       self.textEdit3.setText(str(1.0))

    def lineEdita1(self,newText):
     if newText.length()>0:
      try:
       self.A[0][1]=float(newText.ascii())
      except (TypeError,ValueError):
       self.textEdit4.setText(str(0.0))

    def lineEdita2(self,newText):
     if newText.length()>0:
      try:
       self.A[1][0]=float(newText.ascii())
      except (TypeError,ValueError):
       self.textEdit5.setText(str(0.0))

    def lineEdita3(self,newText):
     if newText.length()>0:
      try:
       self.A[1][1]=float(newText.ascii())
      except (TypeError,ValueError):
       self.textEdit6.setText(str(1.0))

    def lineEditb1(self,newText):
     if newText.length()>0:
      try:
       self.b[0]=float(newText.ascii())
      except (TypeError,ValueError):
       self.textEdit7.setText(str(0.0))

    def lineEditb2(self,newText):
     if newText.length()>0:
      try:
       self.b[1]=float(newText.ascii())
      except (TypeError,ValueError):
       self.textEdit8.setText(str(0.0))

    def accept(self):
     if self.lsm!=None:
       self.lsm.linear_transform(self.A,self.b)
     if self.cview!=None:
       self.cview.updateCanvas()
     QDialog.accept(self)

def main(args):
  app=QApplication(args)
  win=TransDialog(None)
  #win.setTitle("Testing...")
  win.show()
  app.connect(app,SIGNAL("lastWindowClosed()"),
               app,SLOT("quit()"))
  app.exec_loop()

if __name__=="__main__":
   main(sys.argv)
