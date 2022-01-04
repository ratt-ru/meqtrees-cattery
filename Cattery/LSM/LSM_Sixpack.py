#!/usr/bin/env python3


#
#% $Id$ 
#
#
# Copyright (C) 2002-2007
# The MeqTree Foundation & 
# ASTRON (Netherlands Foundation for Research in Astronomy)
# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
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
# along with this program; if not, see <http://www.gnu.org/licenses/>,
# or write to the Free Software Foundation, Inc., 
# 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from Timba.Meq import meq
from Timba import utils
from Timba.TDL import *
import math

_dbg = utils.verbosity(0, name='LSM_Sixpack')
_dprint = _dbg.dprint                    # use: _dprint(2, "abc")
_dprintf = _dbg.dprintf   








##########################################################################################
##########################################################################################
##########################################################################################


class Sixpack:
    """
    Constructors:
    Sixpack(stokesI=sI,stokesQ=sQ,stokesU=sU,stokesV=sV,ra=RA,dec=Dec,label=label): 
    by default, stokesI=1.0 and dec=pi/2 and all other node stubs are zero
    
    Sixpack(stokesI=sI,stokesQ=sQ,stokesU=sU,stokesV=sV,ra=RA,dec=Dec,ns=ns,label=label): 
    by default, stokesI=1.0 and dec=pi/2 and all other node stubs are zero,
    composed into one subtree as well.
    
    Other methods:
    decompose() : decomposes the root into six subtrees
    in composed state, sixpack !=None,
    in decomposed state, sixpack ==None
    sixpack(ns=ns): if already composed, return the sixpack subtree,
    else, first compose it using given nodescope and return it
    iquv(ns=ns): compose the fourpack using the given nodescope 
    and return it or return an already composed subtree
    radec(ns=ns): compose the twopack using the given nodescope and
    return it or return an already composed subtree
    
    stokesI(new_stokesI):
    if called without any input, returns the StokesI,
    else, set StokesI node stub to the new value
    stokesQ(new_stokesQ):
    stokesU(new_stokesQ):
    stokesV(new_stokesQ): same as above stokesI()
    
    ra(new_RA): 
    dec(new_Dec): same as above stokesI()

    Sixpack_Point contains:
    __label: label of root node, if any
    
    node stubs
    __sI:
    __sQ:
    __sU:
    __sV:
    __RA:
    __Dec: six stubs for the six subtrees
    __sixpack: Root of the Sixpack subtree
    __iquv: root of fourpack subtree
    __radec: root of radec subtree
    """
    """
    Constructors:
    Sixpack(root=root,label=label): 
    the above two are mandetory, only in a patch the 'root' input argument
    should be given. Note that by default both are set to None. 
    
    Other methods:
    root() : return the root 
    
    Sixpack_Patch contains:
    __label: label of root node, if any
    __root: root of subtree
    """

    def __init__(self,**pp):
        """Possible input (and defalut values) for the constructor are:
         Sixpack(stokesI=sI,stokesQ=sQ,stokesU=sU,stokesV=sV,ra=RA,dec=Dec,label=label): roots of 
         the six subtrees but not composed
         Sixpack(stokesI=sI,stokesQ=sQ,stokesU=sU,stokesV=sV,ra=RA,dec=Dec,ns=ns): roots of
         the six subtrees, composed into one subtree
        """
        pp.setdefault('label',None)

        pp.setdefault('ns',None)
        pp.setdefault('RA',0)
        pp.setdefault('Dec',math.pi/2)
        pp.setdefault('I0',1.0)
        pp.setdefault('stokesQ',0)
        pp.setdefault('stokesU',0)
        pp.setdefault('stokesV',0)
        pp.setdefault('SI',0)
        pp.setdefault('RM',0)
        pp.setdefault('f0',1e6)
        pp.setdefault('type','point')
        pp.setdefault('lm',None);
        

        self.__label=pp['label']
        self.__sI=pp['I0']
        self.__sQ=pp['stokesQ']
        self.__sU=pp['stokesU']
        self.__sV=pp['stokesV']
        self.__RA=pp['RA']
        self.__Dec=pp['Dec']
        self.__SI=pp['SI']
        self.__f0=pp['f0']
        self.__RM=pp['RM']

        # remember the given nodescope, if any
        self.__ns=pp['ns']
        # do not accept root of sixpack as constructor input
        self.__sixpack=None
        # root of 2pack and 4pack
        self.__radec=None
        self.__iquv=None
        self.__type=pp['type']

        # print my_ns
        # print my_name
        return None




    #-------------------------------------------------------------------------
    # common methods to get/set an item from the Sixpack, if input is None,
    # item is returned, else item is set to new value given as input.

    def ra(self, val=None):
        if val != None:
            self.__RA = val
        return self.__RA

    def dec(self, val=None):
        if val != None:
            self.__Dec = val
        return self.__Dec

    def stokesI(self, val=None):
        if val != None:
            self.__sI = val
        return self.__sI

    def stokesQ(self, val=None):
        if val != None:
            self.__sQ = val
        return self.__sQ

    def stokesU(self, val=None):
        if val != None:
            self.__sU = val
        return self.__sU

    def stokesV(self, val=None):
        if val != None:
            self.__sV = val
        return self.__sV
 
    def nodescope(self, val=None):
        if val != None:
            self.__ns = val
        return self.__ns
 

    def sixpack(self, ns=None):
        return self


    def iquv(self, ns=None):
        """return the 4pack from the six node stubs"""
        return self.__iquv
    
    def lm(self, ns=None):
        return self.__lm


    def radec(self, ns=None):
        """return the 2pack from the six node stubs"""
        return self.__radec

    def rm(self,val=None):
        if val != None:
         self.__RM=val
        return self.__RM

    def SI(self,val=None):
        if val != None:
          self.__SI=val
        return self.__SI

    def f0(self,val=None):
        if val != None:
          self.__f0=val
        return  self.__f0

    def label(self):
       return self.__label

    def ispoint(self):
       return self.__type=='point';

    def display(self, txt=None,full=False):
        print(self.__label)
        return self.__label

    # Generic string
    def __str__(self):
        return self.display()




##################################################################################################
##################################################################################################
##################################################################################################


#############################################################################################
#############################################################################################
def newstar_source (ns=0, predefine=False, flux_att=1.0, slave=False, simul=False, **pp):
   """Make a Sixpack (I,Q,U,V,Ra,Dec) for a source with NEWSTAR parametrisation"""

   # Make the Sixpack and get its ParmSet object:
   punit = pp['punit']
   sixpack = Sixpack(label=punit, **pp)
   if 'parmtable' in pp:
       sixpack.parmtable(pp['parmtable'])
 
   return sixpack




#############################################################################################
#############################################################################################

if __name__=='__main__':
    ns=NodeScope()
    my_name='my_sixpack'
    # create some node stubs for the sixpack
    # first some parameters
    ns.f<<Meq.Parm(meq.polclog([1,0.1,0.01]))
    ns.t<<Meq.Parm(meq.polclog([0.01,0.1,1]))
    # next the node stubs
    stubI=ns['Istub']<<1.1*Meq.Sin(ns.f+ns.t)
    stubQ=ns['Qstub']<<2.0*Meq.Cos(ns.f)
    stubU=ns['Ustub']<<2.1*Meq.Sin(ns.f-2)
    stubV=ns['Vstub']<<2.1*Meq.Cos(ns.f-2)
    stubRA=ns['RAstub']<<2.1*Meq.Cos(ns.f-2)*Meq.Sin(ns.t)
    stubDec=ns['Decstub']<<2.1*Meq.Cos(ns.f-2)*Meq.Sin(ns.t)

    # now create the sixpack
    my_sp=Sixpack(label=my_name, ns=ns,
                              ra=stubRA, dec=stubDec,
                              stokesI=stubI,
                              stokesQ=stubQ,
                              stokesU=stubU,
                              stokesV=stubV)
    my_sp.display()

    # Conversion to 2x2  cohaerency matrix
    print(my_sp.coh22(ns))
    print(my_sp.coh22(ns))
    print(my_sp.coh22(ns, 'circular'))
    print()
  
    my_name='patch0'
    stubR=ns[my_name]<<1.1*Meq.Sin(ns.f+ns.t)+2.0*Meq.Cos(ns.f)-2.1*Meq.Sin(ns.f-2)
 
    my_sp_patch=Sixpack(label=my_name,root=stubR)

    # resolve the forest
    ns.Resolve()
    print("========================= Point ")
    print(dir(my_sp))
    my_sp.display()
    my_sp.decompose()
    print(my_sp.stokesI())
    print(my_sp.stokesQ())
    print(my_sp.stokesU())
    print(my_sp.stokesV())
    print(my_sp.root())
    print(my_sp.type())
    print("========================= Patch")
    print(dir(my_sp_patch))
    my_sp_patch.display()
    my_sp_patch.decompose()
    print(my_sp_patch.stokesI())
    print(my_sp_patch.stokesQ())
    print(my_sp_patch.stokesU())
    print(my_sp_patch.stokesV())
    print(my_sp_patch.root())

#############################################################################################
