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

# this contains parallelization-related options
from Timba.TDL import *

mpi_enable = False;
parallelize_by_source = False;
mpi_nproc = False;

_options = [ 
    TDLMenu('Enable MPI',
              toggle='mpi_enable',
      *[ 
          TDLOption('mpi_nproc',"Number of processors to distribute to",[2,4,8],more=int),
          TDLOption('parallelize_by_source',"Enable parallelization by source",False),
        ]
    )
];

def compile_options ():
  global _options;
  return _options;


# This is a function to add a large number of visibilities in a clever way.
# The problem is that having too many children on a node leads to huge cache 
# usage. So instead we make a hierarchical tree to only add N things at a time.
#
# 'nodes' are output (sum) nodes
# 'visibilities' is a list of nodes containing visibilities (per component)
# 'ifrs' is a list of IFRs (so that for each i,p,q, visibilities[i](p,q) is a valid node)
# 'step' is the number of items to add at a time
# 'kw' is passed as-is to the Meq.Add() node
#
def smart_adder (nodes,visibilities,ifrs,step=8,**kw):
  sums = visibilities;
  nsum = 1;
  while sums:
    # if down to 'step' terms or less, generate final nodes
    if len(sums) <= step:
      for ifr in ifrs:
        nodes(*ifr) << Meq.Add(*[x(*ifr) for x in sums],**kw)
      break;
    # else generate intermediate sums
    else:
      newsums = [];
      for i in range(0,len(sums),step):
        # if we're dealing with one odd term out, propagate it to newsums list as-is
        if i == len(sums)-1:
          newsums.append(sums[i]);
        else:
          newnode = nodes('(%d:%d)'%(i*nsum,min(i+step-1,len(sums)-1)*nsum));  # create unique name for intermediate sum node
          for ifr in ifrs:
            newnode(*ifr) << Meq.Add(*[x(*ifr) for x in sums[i:i+step]],**kw);
          newsums.append(newnode);
      sums = newsums;
      nsum *= step;
  return nodes;

def add_visibilities (nodes,vislist,ifrs):
  """Smart method to add a list of visibility nodes in various clever ways
  (depending on our parallelization settings).
  'nodes' is an unqualified output node.
  'vislist' is a list of unqualified visibilities (presumably, per source)
  'ifrs' is a list of p,q pairs, such that for every vis in vislist, vis(p,q) yields a valid node.
  
  Upon return, for each p,q, nodes(p,q) will contain the sum of vis(p,q) for each vis in vislist
  """
  # If Parallelization is enabled, divide visibilities into batches and
  # place on each processor
  if mpi_enable and parallelize_by_source:
    # how many sources per processor?
    nsrc = len(vislist);
    src_per_proc = [nsrc/mpi_nproc]*mpi_nproc;
    # distriebute remainder over a few processors
    remainder = nsrc%mpi_nproc;
    if remainder != 0:
      for i in range(1,remainder+1):
        src_per_proc[i] += 1;
    # now loop over processors
    per_proc_nodes = [];
    isrc0 = 0;
    for proc in range(mpi_nproc):
      isrc1 = isrc0 + src_per_proc[proc];
      # this node will contain the per-processor sum
      procnode = nodes('P%d'%proc);
      per_proc_nodes.append(procnode);
      # now, make nodes to add contributions of every source on that processor
      smart_adder(procnode,vislist[isrc0:isrc1],ifrs,proc=proc);
      isrc0 = isrc1;
    # now, make one final sum of per-processor contributions
    smart_adder(nodes,per_proc_nodes,ifrs,proc=0,mt_polling=True);
  else:
    # No parallelization, all sourced added up on one machine.
    smart_adder(nodes,vislist,ifrs);
  return nodes;