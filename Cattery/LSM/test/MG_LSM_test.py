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


from Timba.TDL import *
from Timba.Meq import meq
from Timba.LSM.LSM import LSM

# Create Empty LSM - global
lsm=LSM()
def _define_forest (ns):
 global lsm
 home_dir = os.environ['HOME']
 infile_name = home_dir + '/Timba/LSM/test/3C343_nvss.txt'
 lsm.build_from_catalog(infile_name,ns)
 plist=lsm.queryLSM(count=3)
 ## build your own tree
 source_list=[]
 for pu in plist:
    # essential source params
    (ra,dec,sI,sQ,sU,sV,SIn,f0,RM)=pu.getEssentialParms(ns)
    # extended sources
    (eX,eY,eP)=pu.getExtParms()
    ns.I(pu.name)<<Meq.Constant(sI)
    ns.Q(pu.name)<<Meq.Constant(sQ)
    ns.U(pu.name)<<Meq.Constant(sU)
    ns.V(pu.name)<<Meq.Constant(sV)
    ns.IQUV(pu.name)<<Meq.Composer(ns.I(pu.name),ns.Q(pu.name),ns.U(pu.name),ns.V(pu.name))
    source_list.append(ns.IQUV(pu.name))
    

 lsm.display()

#********************************************************************************
# Initialisation and testing routines
# NB: this section should always be at the end of the script
#********************************************************************************


#-------------------------------------------------------------------------
# Meqforest execution routine (may be called from the browser):
# The 'mqs' argument is a meqserver proxy object.

def _test_forest (mqs, parent):
 pass

#### PUnits
def _tdl_job_query_punits(mqs, parent):
 global lsm
 # obtain the punit list of the 3 brightest ones
 plist=lsm.queryLSM(count=3)
 for pu in plist:
  my_sp=pu.getSP()
  my_sp.display()

#####################################################################
#-------------------------------------------------------------------------
# Test routine to check the tree for consistency in the absence of a server

if __name__ == '__main__':
  print('\n*******************\n** Local test of:',script_name,':\n')
  ns=NodeScope()
  _define_forest(ns)
  ns.Resolve()
  print("Added %d nodes" % len(ns.AllNodes()))
  #display LSM without MeqBrowser
  # create cell
  lsm.display(app='create')
#********************************************************************************
#********************************************************************************




