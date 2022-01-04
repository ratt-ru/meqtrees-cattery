#!/usr/bin/env python3
#
# The Local Sky Model (LSM)
#
# This code includes the extra classes used by the LSM internally. So a user
# need not worry about what is included here unless he/she wants trouble.
# 


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
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

import sys,time
import pickle # for serialization and file io
# from Dummy import *

from .common_utils import *
from Timba.Meq import meq
from Timba.TDL import *
from . import LSM_Sixpack
from Timba.Meq import meq

from Timba.Apps import app_nogui
from Timba.Apps import assayer

#############################################
class Source:
 """Source object in source table
 Attributes are
 name=Source name - String
 treeType=type of Template tree - string
 tableName=MeqParm table Name - string
 
 Additional parameters for extended sources
 major, minor, pangle: for Elliptic Gaussians
 filename: for images ??"""

 # Constructor
 def __init__(self,name,treeType=POINT_TYPE,tableName='Table1',major=0, minor=0, pangle=0):
   # Check source name to be a string
   if type(name) != type(""):
    raise TypeError("Name must be a string, not %s"  % type(name).__name__)
   self.name=name
   self.tableName=tableName
   self.eX=major
   self.eY=minor
   self.eP=pangle
   if major==0 and minor==0 and pangle==0:
     # point source
     self.treeType=POINT_TYPE
   else:
     self.treeType=GAUSS_TYPE
     

 # set template tree type
 def setTemplateTree(self,treeType):
   self.treeType=treeType
 
 # return tree type
 def treeType(self):
   return self.treeType
 # return tree type
 def getType(self):
   return self.treeType



 # return extended params
 def extParms(self):
  return (self.eX,self.eY,self.eP)

 # Print
 def __str__(self):
   temp_str="Source name: "+self.name
   temp_str+=" Template tree: "+str(self.treeType)
   temp_str+=" MeqParm Table: "+self.tableName
   temp_str+=" Extended ?: "+str(self.extParms())
   return temp_str
###############################################
# class sixpack is not needed to be defined here (JEN code does that)
# instead we define a class to store a cell and six vellsets
# to store the request/response to the meqtree
# sixpack helper
class SpH:
 """the sixpack helper contains:
    root: root or the composer node of the six subtrees
    lsm: LSM using this sixpack - needed to obtain ObsWin and ObsRes info
    sI:
    sQ:
    sU:
    sV:
    RA:
    Dec: the six vellsets of the current response
 """
 
 # Constructor
 def __init__(self,lsm,Root=None,mqs=None,sI=None,sQ=None,sU=None,sV=None,RA=None,Dec=None):
  self.root=Root
  self.lsm=lsm
  self.sI=sI
  self.sQ=sQ
  self.sU=sU
  self.sV=sV
  self.RA=RA
  self.Dec=Dec
  # static value for debugging
  self.static_RA=0
  self.static_Dec=0

 # traditional getter and setter methods
 # setter methods need the vellset as input
 # getter methods return the value array
 def setI(self,I):
  self.sI=I
 def getI(self):
  if self.sI==None:
   return 0
  return self.sI.value
 def setQ(self,Q):
  self.sQ=Q
 def getQ(self):
  if self.sQ==None:
   return 0
  return self.sQ.value
 def setU(self,U):
  self.sU=U
 def getU(self):
  if self.sU==None:
   return 0
  return self.sU.value
 def setV(self,V):
  self.sV=V
 def getV(self):
  if self.sU==None:
   return 0
  return self.sV.value
 def setRA(self,RA):
  self.RA=RA
 def getRA(self):
  if self.RA==None:
   return self.static_RA
  return self.RA.value[0]
#  return self.static_RA
 def setDec(self,Dec):
  self.Dec=Dec 
 def getDec(self):
  if self.Dec==None:
   return self.static_Dec 
  return self.Dec.value[0]
#  return self.static_Dec

 ### Debugging
 def set_staticRA(self,RA):
  self.static_RA=RA
# def get_staticRA(self):
#  return self.static_RA
 def set_staticDec(self,Dec):
  self.static_Dec=Dec 
# def get_staticDec(self):
#  return self.static_Dec


 # setting the meqtree root
 def setRoot(self,root):
  self.root=root
 # link back to the LSM
 def setLSM(self,lsm):
  self.lsm=lsm

 # update values according to the new cell, 
 # the name of the p-unit is given by 'pname'
 def updateValues(self,pname):
  if (self.lsm!=None) and (self.lsm.cells!=None) and\
   (self.lsm.mqs!=None):
   # create request object
   my_request = meq.request(cells=self.lsm.cells, rqtype='ev')
   # get the correct name from the sixpack object
   punit=self.lsm.getPUnit(pname)
   psixpack=punit.getSP()
   my_name=''
   if psixpack!=None and psixpack.ispoint():
    my_name=psixpack.sixpack().name
   elif psixpack!=None: # patch
    my_name=psixpack.root().name
   elif punit.sp!=None:
    my_name='sixpack:q='+pname
   my_args=meq.record(name=my_name, request=my_request)
   my_result=self.lsm.mqs.meq('Node.execute', my_args,wait=True)
   # update values only if no 'Fail has happened'
   # if you try to create patches with 1 source, it will fail
   self.setRA(my_result.result.vellsets[0])
   self.setDec(my_result.result.vellsets[1])
   self.setI(my_result.result.vellsets[2])
   self.setQ(my_result.result.vellsets[3])
   self.setU(my_result.result.vellsets[4])
   self.setV(my_result.result.vellsets[5])
  else:
   print("Cannot update values: (cell,mqs)",self.lsm.cells,self.lsm.mqs)


 # returns the size of the values corresponding to the given (t,f) pair
 # and the quantity 'A','I','Q','U','V','RA','Dec'
 # returns [l,m] for lxm array size
 def getValueSize(self,type,freq_index,time_index):
  # range error check
  if (self.lsm.cells.segments.freq.start_index > freq_index) or\
    (self.lsm.cells.segments.freq.end_index < freq_index):
    print("Index error, Frequency %d" %freq_index)
    freq_index=self.lsm.cells.segments.freq.start_index
  if (self.lsm.cells.segments.time.start_index > time_index) or\
    (self.lsm.cells.segments.time.end_index < time_index):
    print("Index error, Time %d" %time_index)
    time_index=self.lsm.cells.segments.time.start_index
  if type=='A' or type=='I':
   # expect a vellset
   try:
    shape=self.sI.shape
    if len(shape)==4:
     l=shape[2]
     m=shape[3]
    else:
     l=1
     m=1
    return [l,m]
   except:
    # no set, just a scalar
    # print "Scalar Return "
    return [1,1] 
  elif type=='Q':
   try:
    shape=self.sQ.shape
    if len(shape)==4:
     l=shape[2]
     m=shape[3]
    else:
     l=1
     m=1
    return [l,m]
   except:
    # no set, just a scalar
    # print "Scalar Return "
    return [1,1] 
  elif type=='U':
   try:
    shape=self.sU.shape
    if len(shape)==4:
     l=shape[2]
     m=shape[3]
    else:
     l=1
     m=1
    return [l,m]
   except:
    # no set, just a scalar
    # print "Scalar Return "
    return [1,1] 
  elif type=='V':
   try:
    shape=self.sV.shape
    if len(shape)==4:
     l=shape[2]
     m=shape[3]
    else:
     l=1
     m=1
    return [l,m]
   except:
    # no set, just a scalar
    # print "Scalar Return "
    return [1,1] 

  # will not get here
  return [1,1]

 # returns the value(s) corresponding to the given (t,f) pair
 # and the quantity 'A','I','Q','U','V','RA','Dec'
 # note that in case of a patch, this would return a 2D
 # array in l,m coordinates
 def getValue(self,type,freq_index,time_index):
  # range error check
  if (self.lsm.cells.segments.freq.start_index > freq_index) or\
    (self.lsm.cells.segments.freq.end_index < freq_index):
    print("Index error, Frequency %d" %freq_index)
    freq_index=self.lsm.cells.segments.freq.start_index
  if (self.lsm.cells.segments.time.start_index > time_index) or\
    (self.lsm.cells.segments.time.end_index < time_index):
    print("Index error, Time %d" %time_index)
    time_index=self.lsm.cells.segments.time.start_index

  if type=='A' or type=='I':
   # expect a vellset
   try:
    shape=self.sI.shape
    if shape[0]==1:
     # no time dependence
     time_index=0
    if shape[1]==1:
     # no freq. dependence
     freq_index=0
    #print "Return ",self.sI.value[time_index][freq_index]
    return self.sI.value[time_index][freq_index]
   except:
    # no set, just a scalar
    #print "Scalar Return ",self.sI.value[0]
    return self.sI.value[0]
  elif type=='Q':
   # expect a vellset
   try:
    shape=self.sQ.shape
    if shape[0]==1:
     # no time dependence
     time_index=0
    if shape[1]==1:
     # no freq. dependence
     freq_index=0
    #print "Return ",self.sI.value[time_index][freq_index]
    return self.sQ.value[time_index][freq_index]
   except:
    # no set, just a scalar
    #print "Scalar Return ",self.sI.value[0]
    return self.sQ.value[0]
  elif type=='U':
   # expect a vellset
   try:
    shape=self.sU.shape
    if shape[0]==1:
     # no time dependence
     time_index=0
    if shape[1]==1:
     # no freq. dependence
     freq_index=0
    #print "Return ",self.sI.value[time_index][freq_index]
    return self.sU.value[time_index][freq_index]
   except:
    # no set, just a scalar
    #print "Scalar Return ",self.sI.value[0]
    return self.sU.value[0]
  elif type=='V':
   # expect a vellset
   try:
    shape=self.sV.shape
    if shape[0]==1:
     # no time dependence
     time_index=0
    if shape[1]==1:
     # no freq. dependence
     freq_index=0
    #print "Return ",self.sI.value[time_index][freq_index]
    return self.sV.value[time_index][freq_index]
   except:
    # no set, just a scalar
    #print "Scalar Return ",self.sI.value[0]
    return self.sV.value[0]
  else:
    print("Error request",type)
    return 0




 
 # Print
 def __str__(self):
   temp_str="SPHelper: {I="+str(self.getI())
   temp_str+=",Q="+str(self.getQ())
   temp_str+=",U="+str(self.getU())
   temp_str+=",V="+str(self.getV())
   temp_str+=",RA="+str(self.getRA())
   temp_str+=",Dec="+str(self.getDec())
   temp_str+="}"
   return temp_str

 # clone this SpH without circular reference to the LSM
 # and references to the MeqTree system for saving
 def clone(self):
  newsph=SpH(None)
  newsph.root=None
  #newsph.sI=self.sI
  #newsph.sQ=self.sQ
  #newsph.sU=self.sU
  #newsph.sV=self.sV
  #newsph.RA=self.RA
  #newsph.Dec=self.Dec
  # static value for debugging
  newsph.static_RA=self.static_RA
  newsph.static_Dec=self.static_Dec
  return newsph


 
###############################################
class TemTree:
 """Template tree object"""
 
 # Constructor
 def __init__(self,id):
  self.id=id 

#########################################################################
