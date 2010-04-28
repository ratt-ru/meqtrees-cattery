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

import math
from Timba.Meq import meq
from Timba.TDL import *
import Timba.array
##############################################
### common definitions for the GUI
### and utility functions
##############################################

# GUI Zoom Modes
GUI_ZOOM_NONE=0
GUI_ZOOM_WINDOW=1 # turn on by menu selection
GUI_ZOOM_START=2 # turn on by first mouse click, after ZOOM_WINDOW
GUI_SELECT_WINDOW=3 # turn on by menu, window for selecting sources
GUI_SELECT_START=4 # turn on by first mouse click, after SELECT_WINDOW
GUI_MOVE_START=5 # turn on by menu selection

# after the first mouse  click, in ZOOM_WINDOW mode, draw zoom window (show())
# after the mouse released in this mode, set zoom mode to GUI_ZOOM_WINDOW
# and zoom, hide zoom window. by menu selection, can set to GUI_ZOOM_NONE


# define constants for RTTI values of special canvasview objects
POINT_SOURCE_RTTI=1001
PATCH_RECTANGLE_RTTI=1002
PATCH_IMAGE_RTTI=1003

# PUnit/Source types
POINT_TYPE= 0
PATCH_TYPE= 1
GAUSS_TYPE= 2
IMAGE_TYPE= 3

# default depths  for various items on the canvas
IMAGE_DEPTH=-1


# column indices for PUnit table
PCOL_NAME=0
PCOL_TYPE=1
PCOL_SLIST=2
PCOL_CAT=3
PCOL_BRIGHT=4
PCOL_FOV=5
PCOL_I=6
PCOL_Q=7
PCOL_U=8
PCOL_V=9
PCOL_RA=10
PCOL_DEC=11

# options for file export
EXPORT_NONE=0
EXPORT_IMG_EPS=101
EXPORT_IMG_PNG=102
EXPORT_IMG_BMP=103
EXPORT_IMG_KARMA=104
EXPORT_PT_EPS=201
EXPORT_PT_TEX=202
EXPORT_PT_TXT=203
EXPORT_ST_EPS=301
EXPORT_ST_TEX=302
EXPORT_ST_TXT=303


######## binary search
def bin_search(xarr,x,i_start,i_end):
  # xarr: [-inf,...,+inf] array so x must be within, also xarr 
  #  must be sorted
  # x: the value to search for
  # i_start,i_end: indices of the array to search for x,
  #  at start, these will be [0,len(xarr)-1]
  # return k, such that x_{k} <= x < x_{k+1}
  # in case an error, return -1

  # trivial case
  if i_start==i_end:
   if  xarr[i_start]==x:
     return i_start
   else:
     print "bin search error 1"
     return -1

  # trivial case with length 2 array
  if i_end==i_start+1:
   if x>=xarr[i_start] and\
     x<xarr[i_end]:
    return i_start
   else:
    print "bin search error 2"
    return -2

  # compare with the mid point
  i=(i_start+i_end)/2
  if x >=xarr[i] and\
     x<xarr[i+1]:
   return i
  else:
   # go to lower half or the upper half of array
   if x < xarr[i]:
    return bin_search(xarr,x,i_start,i)
   else:
    return bin_search(xarr,x,i,i_end)


  # will not get here
  print "bin search error 3"
  return -3



####### Radians to RA=[hr,min,sec]
## Rad=(hr+min/60+sec/60*60)*pi/12
def radToRA(rad):
  # convert negative values
  if rad <0:
   rad=rad+2*math.pi

  tmpval=rad*12.0/math.pi
  hr=int(tmpval)
  tmpval=tmpval-hr
  tmpval=tmpval*60
  mins=int(tmpval)
  tmpval=tmpval-mins
  tmpval=tmpval*60
  sec=tmpval
  return [hr%24,mins%60,sec]

####### Radians to Dec=[hr,min,sec]
## Rad=(hr+min/60+sec/60*60)*pi/180
def radToDec(rad):
  if rad<0:
   mult=-1;
   rad=abs(rad)
  else:
   mult=1
  tmpval=rad*180.0/math.pi
  hr=int(tmpval)
  tmpval=tmpval-hr
  tmpval=tmpval*60
  mins=int(tmpval)
  tmpval=tmpval-mins
  tmpval=tmpval*60
  sec=tmpval
  return [mult*(hr%180),mins%60,sec]


########## metric form
### [0,1000)  -
### [1000,1000000) /10^3 k
### [10^6,10^9) /10^6 M
### [10^9,10^12) /10^9 G
### [10^12,..) /10^12 T
# input value=val, format string=format
# output, value+suffix (k,M,G,T)
def stdForm(val,format):
 vl=abs(val)
 if vl <1e3:
  tmpstr=" "
  vlformat=format%vl
 elif vl<1e6:
  tmpstr="k"
  vlformat=format%(vl/1e3)
 elif vl<1e9:
  tmpstr="M"
  vlformat=format%(vl/1e6)
 elif vl<1e12:
  tmpstr="G"
  vlformat=format%(vl/1e9)
 else:
  tmpstr="T"
  vlformat=format%(vl/1e12)

 return [vlformat,tmpstr]


###############################################################
# the following methods are used to 
# serialize trees by brute force. In fact,
# no serialization is done, but only the essence required
# to recreate the whole tree is stored as (recursive) dictionaries.
def rec_parse(myrec):
  """ recursively parse init record 
      and construct a dictionary
  """
  #print myrec
  my_keys=myrec.keys()
  new_dict={}
  for kk in my_keys:
   if kk=='funklet' or kk=='init_funklet':
     ff=myrec[kk]
     if is_meqpolclog(ff):
       #print "create polclog"
       gg={'type':'meqpolclog', 'funk':serialize_funklet(ff)}
     else:
       gg={'type':'meqpolc', 'funk':serialize_funklet(ff)}
    
     #print gg
     new_dict[kk]=gg
   elif isinstance(myrec[kk],meq.record):
     new_dict[kk]=rec_parse(myrec[kk])
   elif isinstance(myrec[kk],type(Timba.array.array(0))): 
     # just serialize the value
     #new_dict[kk+'_isarray']=pickle.dumps(myrec[kk])
     #print myrec[kk].__class__
     #print pickle.dumps(myrec[kk])
     new_dict[kk]=myrec[kk]
   else:
     new_dict[kk]=myrec[kk]
  return new_dict

def traverse(root,node_dict,subscope=None):
  # if subscope  is not None, that means the nodes
  # have a leadins 'subscope':: part in their names. 
  # so we need to strip it before saving
  #print "Subscope=",subscope
  chlist=root.children
  name=root.name
  if subscope:
   nlist=string.split(name,"::")
   name=nlist.pop()# get last item
 
  #print "My name",name
  classname=root.classname
  if not node_dict.has_key(name):
   node_dict[name]={'name':name, 'classname':classname, 'initrec':{},\
        'children':[]}  
   ir=root.initrec()
   #print ir
   myrec=node_dict[name]
   myrec['initrec']=rec_parse(ir)
   # if any children, traverse
   for idx,ch in chlist:
    traverse(ch,node_dict,subscope)
    if subscope:
      nlist=string.split(ch.name,"::")
      chname=nlist.pop()# get last item
      myrec['children'].append(chname)
    else:
      myrec['children'].append(ch.name)
   #print "Rec ",node_dict[name]

def create_node_stub(mydict,stubs,ns,myname):
  myrec=mydict[myname]
  # first, if this node has any children
  # and if they have not being created,
  # create them
  chlist=myrec['children']
  #stublist=[]
  for ch in chlist:
   if not stubs.has_key(ch):
    stubs[ch]=create_node_stub(mydict,stubs,ns,ch)
  # now we have created the child list
  
  # first check the nodescope if the node given by the name already exists
  if ns[myname].initialized():
   print "WARNING: node %s already created"%myname
   got_stub=cname_node_stub(ns,myname)
   # add children
   for child in chlist:
    # handle differently if this is a sub scope
    if ns._name:
      got_stub.add_children(ns._name+'::'+child)
    else:
      got_stub.add_children(child)
   return got_stub

  # now deal with initrec()
  irec=myrec['initrec']
  #print 'My Rec==',myrec
  #print 'Init Rec==',irec
  myclass=myrec.pop('classname')
  mygentype=myclass.replace("Meq","")
  if len(chlist)>0:
   # handle differently if it is a subscope
   if ns._name:
     full_chlist=[]
     for child in chlist:
       full_chlist.append(ns._name+'::'+child)
   else:
       full_chlist=chlist
   fstr="ns['"+myname+"']<<Meq."+mygentype+'(children='+str(full_chlist)+','
  else:
   fstr="ns['"+myname+"']<<Meq."+mygentype+'('
  irec_str=""
  # Remove JUNK! from initrecord()
  # remove class field
  if irec.has_key('class'):
   irec.pop('class')
  # remove node_desctiption
  if irec.has_key('node_description'):
   irec.pop('node_description')
  # remove name
  if irec.has_key('name'):
   irec.pop('name')
  # remove node index
  if irec.has_key('nodeindex'):
   irec.pop('nodeindex')
  # remove children
  if irec.has_key('children'):
   irec.pop('children')



  for kname in irec.keys():
   krec=irec[kname]
   if not isinstance(krec,dict):
    irec_str=irec_str+" "+kname+"="+str(krec)+','
   else:
    if (kname=='init_funklet'):
     if krec['type']=='meqpolclog':
      irec_str=irec_str+" "+"init_funklet="+"default_funklet_value,"
      # deserialize the value
      default_funklet_value=krec['funk']
      #print dir(default_funklet_value)
     else:# assume to be a meqpolc()
      irec_str=irec_str+" "+"init_funklet="+"default_funklet_value,"
      # deserialize the value
      default_funklet_value=krec['funk']
    elif (kname=='funklet'):
     if krec['type']=='meqpolclog':
      irec_str=irec_str+" "+"funklet="+"default_funklet_value,"
      # deserialize the value
      default_funklet_value=krec['funk']
      #print dir(default_funklet_value)
     else:# assume to be a meqpolc()
      irec_str=irec_str+" "+"funklet="+"default_funklet_value,"
      # deserialize the value
      default_funklet_value=krec['funk']
 
     #print default_funklet_value


  total_str=fstr+irec_str+')'
  # MeqParm is special
  #if myclass.lstrip('Meq')=='Parm':
  # total_str="ns['"+myname+"']<<Meq.Parm(default_funklet_value)"
  #print "Total=",total_str
  #exec total_str in globals(),locals()
  exec total_str
  #print ns[myname].initrec()
  #return ns[myname]
  return cname_node_stub(ns,myname)
     
 
# the basic assumption with the following method is 
# the forest has no circular references
def reconstruct(my_dict,ns):
  # temp dictionary to store created node stubs
  nodestub_dict={}
  for sname in my_dict.keys():
   if not nodestub_dict.has_key(sname):
     nodestub_dict[sname]=create_node_stub(my_dict,nodestub_dict,ns,sname)
  return nodestub_dict

def is_meqpolclog(obj):
 if str(obj.__class__)=="<class 'Timba.dmi.MeqPolcLog'>":
  return True
 else:
  return False
def is_meqpolc(obj):
 if str(obj.__class__)=="<class 'Timba.dmi.MeqPolc'>":
  return True
 else:
  return False

def serialize_funklet(fnklt):
 #print fnklt
 return fnklt
 coeff=fnklt['coeff'] # this is a Timba.array
 if coeff.size()>1:
  return coeff.tolist()
 else: # return scalar
  #print "coeff=",coeff
  #print "sum=",coeff.real
  return coeff.real
 # is we have a compiled funklet, we can return a dict
 # to be done
  
#########################################################################
#########################################################################
# utility function to extract MeqParm information from a TDL_Sixpack
# sixpack=sixpack object (in composed form)
# ns=nodescope used for the sixpack
def extract_parms(sixpack,ns):
 myra=sixpack.ra()
 mydec=sixpack.dec()
 mybr=sixpack.stokesI()
 #print "RA,Dec,App.Bri=",myra,mydec,mybr
 return [myra,mydec,mybr]

# utility function to extract polarization percentage [Q,U,V]
# MeqParm information from a TDL_Sixpack
# sixpack=sixpack object (in composed form)
# ns=nodescope used for the sixpack
def extract_polarization_parms(sixpack,ns,absolute=0):
 myii=sixpack.stokesI()
 myqq=sixpack.stokesQ()
 myuu=sixpack.stokesU()
 myvv=sixpack.stokesV()

 if absolute==0:
  ##print "Q,U,V=",myqq,myuu,myvv
  return [myqq,myuu,myvv]

 
 RM=sixpack.rm()
 SI=sixpack.SI()
 f0=sixpack.f0()
 
 return [myii,myqq,myuu,myvv,SI,f0,RM]
 


#########################################################################
# utility function to extract MeqParm information from a TDL_Sixpack
# sixpack=sixpack object (in composed form)
# ns=nodescope used for the sixpack
def node_extract_parms(sixpack,ns):
 # get name
 myname=sixpack.label()
 if ns._name:
  ra_name=ns._name+'::ra:q='+myname
  dec_name=ns._name+'::dec:q='+myname
  I0_name=ns._name+'::I0:q='+myname
  SIFI_name=ns._name+'::SIF_stokesI:q='+myname
 else:
  ra_name='ra:q='+myname
  dec_name='dec:q='+myname
  I0_name='I0:q='+myname
  SIFI_name='SIF_stokesI:q='+myname
 
 #print "looking for params of ",myname
 # RA
 allnodes=ns.Repository()
 ra=allnodes[ra_name]
 # look for the meqparms in this node
 myra=get_default_parms(ra)
 dec=allnodes[dec_name]
 mydec=get_default_parms(dec)
 if allnodes.has_key(I0_name):
  br=allnodes[I0_name]
  mybr=get_default_parms(br)
 elif allnodes.has_key(SIFI_name):
  br=allnodes[SIFI_name]
  mybr=get_default_parms(br)
 else:
  mybr=0.0

 #print "RA,Dec,App.Bri=",myra,mydec,mybr
 return [myra,mydec,mybr]

# utility function to extract polarization percentage [Q,U,V]
# MeqParm information from a TDL_Sixpack
# sixpack=sixpack object (in composed form)
# ns=nodescope used for the sixpack
def node_extract_polarization_parms(sixpack,ns,absolute=0):
 # get name
 print dir(sixpack)
 myname=sixpack.label()
 #print "looking for QUV of ",myname
 # first get I0 for further processing
 if ns._name:
  Q_name=ns._name+'::Qpct:q='+myname
  sQ_name=ns._name+'::stokesQ:q='+myname
  U_name=ns._name+'::Upct:q='+myname
  sU_name=ns._name+'::stokesU:q='+myname
  V_name=ns._name+'::Vpct:q='+myname
  sV_name=ns._name+'::stokesV:q='+myname
 else:
  Q_name='Qpct:q='+myname
  sQ_name='stokesQ:q='+myname
  U_name='Upct:q='+myname
  sU_name='stokesU:q='+myname
  V_name='Vpct:q='+myname
  sV_name='stokesV:q='+myname
 

 allnodes=ns.Repository()

 if ns._name:
  SIFI_name=ns._name+'::SIF_stokesI:q='+myname
  I0_name=ns._name+'::I0:q='+myname
  RM_name=ns._name+'::RM:q='+myname
 else:
  SIFI_name='SIF_stokesI:q='+myname
  I0_name='I0:q='+myname
  RM_name='RM:q='+myname
 if allnodes.has_key(I0_name):
   sI=allnodes[I0_name]
   myii=get_default_parms(sI)
   SI=0
   f0=-1
 elif allnodes.has_key(SIFI_name): 
   sI=allnodes[SIFI_name]
   [myii,SI,f0]=get_stokes_I_parms(sI)
 else:
   myii=SI=0
   f0=-1
 
 if allnodes.has_key(Q_name):
   qq=allnodes[Q_name]
   myqq=get_default_parms(qq)*myii/100
 elif allnodes.has_key(sQ_name):
   qq=allnodes[sQ_name]
   myqq=get_default_parms(qq)
 else:
   myqq=0
 if allnodes.has_key(U_name):
   uu=allnodes[U_name]
   myuu=get_default_parms(uu)*myii/100
 elif allnodes.has_key(sU_name):
   uu=allnodes[sU_name]
   myuu=get_default_parms(uu)
 else:
   myuu=0
 if allnodes.has_key(V_name):
   vv=allnodes[V_name]
   myvv=get_default_parms(vv)*myii/100
 elif allnodes.has_key(sV_name):
   vv=allnodes[sV_name]
   myvv=get_default_parms(vv)
 else:
   myvv=0


 if absolute==0:
  ##print "Q,U,V=",myqq,myuu,myvv
  return [myqq,myuu,myvv]


 if allnodes.has_key(RM_name):
   rmnode=allnodes[RM_name]
   RM=get_default_parms(rmnode)
 else:
   RM=0
 
 return [myii,myqq,myuu,myvv,SI,f0,RM]
 


#
# utility to extract the default parms (if any) of a given MeqParm
# node. It will return the following in the given order, if they
# exist.
# 1. default_value
# 2. First coefficient of init_funklet
# 3. First coefficient of funklet
# 4. 0.0
def get_default_parms(nd):
 # get initrec
 irec=nd.initrec()
 #print irec
 # try to get default_value
 if irec.has_key('default_value'):
  cf=irec['default_value']
  # this need to be a scalar
  if cf.size>=1:
     #my_val=cf.tolist().pop(0)
     my_val=Timba.array.ravel(cf)[0]
  else: #scalar
     my_val=float(cf)
 elif irec.has_key('init_funklet'):
  # get coeff array (or scalar)
  fn=irec['init_funklet']
  if(is_meqpolc(fn)):
   cf=fn['coeff'] 
   if cf.size>=1:
     my_val=Timba.array.ravel(cf)[0]
   else: #scalar
     my_val=float(cf)
  elif(is_meqpolclog(fn)):
   cf=fn['coeff'] 
   if cf.size>=1:
     my_val=Timba.array.ravel(cf)[0]
   else: #scalar
     my_val=float(cf)
   # return exponent
   my_val=math.pow(10,my_val)
  else: # error
   print "WARNING: unable to find a value for funklet",fn
   my_val=-1
 else: # error
   my_val=-1
 return my_val


#
# utility to get flux, spectral index, rest frequency
# from a Stokes I node.
def get_stokes_I_parms(nd):
    irec=nd.initrec()
    if irec.has_key('init_funklet'):
       fn=irec['init_funklet']
       if (is_meqpolc(fn)):
         cf=fn['coeff']
         if cf.size>=1:
            my_I=Timba.array.ravel(cf)[0]
         else: #scalar
            my_I=float(cf)
         my_SI=0
       elif (is_meqpolclog(fn)):
         cf=fn['coeff']
         if cf.size>=1:
           my_I=Timba.array.ravel(cf)[0]
         else: #scalar
           my_I=float(cf)
         # return exponent
         my_I=math.pow(10,my_I)
         if cf.size>=2:
           ## NOTE: the SI can have a polynomial, but for the moment we ignore 
           ## higher order terms
           my_SI=Timba.array.ravel(cf)[1]
         else: #scalar
           my_SI=0
       fr=fn['axis_list']
       if fr and fr.has_key('freq'):
         my_F=fr['freq']
       else:
         my_F=0
       return  [my_I,my_SI,my_F]
    else:
      print "WARNING: unable to find a value for Stokes I in ",nd.name
    return [0,0,0]


######################################################
## extract a node stub from given nodescope using the 
## full name in the format 'a':q='b'
def cname_node_stub(ns,nodename):
  # first strip any subscopes
  sblist=string.split(nodename,"::")
  nodename=sblist.pop()# get last item
  alist=string.split(nodename,":")
  nodestub=None
  if len(alist)==1:
   nodestub=ns[alist[0]]
  else:
    # we have qualifiers
    wstr="nodestub=ns['"+alist.pop(0)+"']"
    for qstr in alist:
      # try to split on the '=' sign
      #blist=string.split(qstr,'=')
      #if len(blist)==1:
      #  wstr=wstr+"'"+blist[0]+"',"
      #else:
      #  wstr=wstr+blist[0]+"='"+blist[1]+"',"
      wstr=wstr+"('"+qstr+"')"

    #print wstr
    exec wstr
  return nodestub


######################################################
## strip subscope name from a given node name
def strip_subscope(nodename):
  nlist=string.split(nodename,"::")
  return nlist.pop()

#######################################################
### change value of MeqParm
def change_parm(ns,node,new_value):
 # first delete old node 
 #try:
 #if ns.__hasattr__(nodename):
 #  #ns.__delattr__(nodename)
 #  pass
 #gg=ns._repository.pop(nodename)
 #del gg
 #ns.Resolve()
 #except AttributeError:
 # print "WARNING: no node exists with name ",nodename
 # pass
 # now recreate a new node
 #work_str="ns['"+nodename+"']<<Meq.Parm(meq.polc("+str(new_value)+"))"
 #print work_str
 #exec work_str in globals(),locals()
 #print "trying to change params of ",node.name
 g=node.initrec()
 #print g
 g['init_funklet']=meq.polc(new_value)
 #ns.Resolve()

# change MeqParm of TDL_Sixpack_Point RA,Dec 
def change_radec(sixpack,new_ra,new_dec,ns):
 myname=sixpack.label()
 ra=sixpack.ra()
 change_parm(ns,ra,new_ra)
 dec=sixpack.dec()
 change_parm(ns,dec,new_dec)

# change MeqParm of TDL_Sixpack_Patch RA,Dec 
def change_radec_patch(patch_name,new_ra,new_dec,ns):
 ra=ns.ra0(q=patch_name)
 change_parm(ns,ra,new_ra)
 dec=ns.dec0(q=patch_name)
 change_parm(ns,dec,new_dec)




################################################################
## convert l,m coordinates to RA,Dec coordinates
## see wng/wnmccv.for WNMCLM for more detail
def lm_to_radec(ra0,dec0,l,m):
    sind0=math.sin(dec0)
    cosd0=math.cos(dec0)
    dl=l
    dm=m
    d0=dm*dm*sind0*sind0+dl*dl-2*dm*cosd0*sind0
    sind=math.sqrt(abs(sind0*sind0-d0))
    cosd=math.sqrt(abs(cosd0*cosd0+d0))
    if (sind0>0):
     sind=abs(sind)
    else:
     sind=-abs(sind)

    dec=math.atan2(sind,cosd)

    if l!=0:
     ra=math.atan2(-dl,(cosd0-dm*sind0))+ra0
    else:
     ra=math.atan2((1e-10),(cosd0-dm*sind0))+ra0


    return (ra,dec)


## convert ra,dec to lm (NCP)
def radec_to_lm(ra0,dec0,ra,dec):
    l=-math.sin(ra-ra0)*math.cos(dec)
    sind0=math.sin(dec0)
    if sind0 != 0:
     m=-(math.cos(ra-ra0)*math.cos(dec)-math.cos(dec0))/math.sin(dec0)
    else:
     m=0
    return (l,m)

## convert ra,dec to lm (SIN)
def radec_to_lm_SIN(ra0,dec0,ra,dec):
    l=-math.sin(ra-ra0)*math.cos(dec)
    m=-(math.cos(ra-ra0)*math.cos(dec)*math.sin(dec0)-math.cos(dec0)*math.sin(dec))
    return (l,m)

#################################################################
if __name__ == '__main__':
  ns=NodeScope()
  my_sixpack=LSM_Sixpack.newstar_source(ns,punit="foo",I0=1.0, f0=1e6,RA=2.0, Dec=2.1,trace=0, Qpct=0.1, Upct=1,Vpct=-0.1)
  my_sixpack.sixpack(ns)
  ns.Resolve()
  print ns.AllNodes()
  print ns["ra:q="+my_sixpack.label()]
  print my_sixpack.ra().name
  my_sixpack.display()
  rr=my_sixpack.ra()
  print rr.name
  print extract_parms(my_sixpack,ns)
  print extract_polarization_parms(my_sixpack,ns)
  change_radec(my_sixpack,3.0,4.5,ns)
  my_sixpack.display()
  print dir(ns)
  print extract_parms(my_sixpack,ns)


  ra0=1.111
  dec0=0.220
  l=0.003
  m=0.4330
  (rr,dd)=lm_to_radec(ra0,dec0,l,m)
  (L,M)=radec_to_lm(ra0,dec0,rr,dd)
  print rr,dd,L,M,l,m
