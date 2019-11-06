# -*- coding: utf-8 -*-
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

from Timba.dmi import *
from Timba.TDL import *

_pages = [];

class Page (object):
  def __init__ (self,name,hplots=2,vplots=3,folder=None):
    self.name = name;
    self.hsize = hplots;
    self.vsize = vplots;
    # find bookmarks list in Settings, or create new one
    if folder is None:
      self._bklist = Settings.forest_state.get('bookmarks',None);
      if self._bklist is None:
        self._bklist = Settings.forest_state.bookmarks = [];
    else:
      self._bklist = folder;
    # page number, used when making multiple pages
    self._pagenum = 0;
    # make a new page record
    self._new_pagerec();
      
  def _new_pagerec (self):
    """creates a new, empty page record. Also used to advance a page.""";
    # make name 
    name = self.name;
    self._pagenum += 1;
    if self._pagenum > 1:
      name += ", Page "+str(self._pagenum);
    self._ix = self._iy = 0;
    # create basic record, but defer adding it to bookmarks
    # until we get at least one entry on the page
    self._pagelist = [];
    self._pagerec = record(name=name,page=self._pagelist);
    self._pagerec_added = False;
    
  def _add_pagerec (self):
    """adds page record to bookmarks, if not already added""";
    if not self._pagerec_added:
      self._bklist.append(self._pagerec);
      self._pagerec_added = True;
      
  def add (self,node,viewer="Result Plotter"):
    """Adds a panel to this page. 'node' is a node object or a node name.
    'viewer' may be used to override the default viewer type""";
    # resolve argument to node name
    if is_node(node):
      node = node.name;
    elif not isinstance(node,str):
      raise TypeError("node or node name expected");
    # add page to bookmarks
    self._add_pagerec();
    # add bookmark
    self._pagelist.append(
      record(viewer=viewer,udi="/node/"+node,pos=(self._iy,self._ix))
    );
    self._ix += 1;
    if self._ix >= self.hsize:
      self._ix = 0;
      self._iy += 1;
      # have we filled up a page? start a new one
      if self._iy >= self.vsize:
        self._new_pagerec();
    return self;
        
  def __lshift__ (self,node):
    return self.add(node);

def PlotPage (name,*rows):
  bklist = [];
  
  if not isinstance(name,str):
    # must be just another row...
    rows = [name] + list(rows);
  
  for irow,onerow in enumerate(rows):
    for icol,node in enumerate(onerow):
      bklist.append(record(
        viewer="Result Plotter",
        udi="/node/"+node,
        pos=(irow,icol)));
        
  if not isinstance(name,str):
    return bklist;
  else:
    return record(name=name,page=bklist);

class Folder (object):
  def __init__ (self,name,folder=None):
    self.name = name;
    # find bookmarks list in Settings, or create new one
    self._bklist = folder if folder is not None else Settings.forest_state.get('bookmarks',None);
    if self._bklist is None:
      self._bklist = Settings.forest_state.bookmarks = [];
    self._folder_list = None;
    
  def _flist (self):
    if self._folder_list is None:
      self._folder_list = [];
      self._bklist.append(record(name=self.name,folder=self._folder_list));
    return self._folder_list;
  
  def page (self,*args,**kw):
    return Page(folder=self._flist(),*args,**kw);
  
  def subfolder (self,*args,**kw):
    return Folder(folder=self._flist(),*args,**kw);
  
def _int_or_str(x):
  """helper function, converts string to integer if possible""";
  try: return int(x);
  except: return x;

import math

def _make_folder (folder,nodes,npp,maxmenu):
  """Does the real work of _make_node_folder. If there's over maxmenu
  entries in the nodes list, recursively calls itself to create subfolders.""";
#  print "_make_folder",folder.name,len(nodes);
  # this is how many menu entries we're going to need
  nent = int(math.ceil(len(nodes)/float(npp)));
  # if this is few enough, insert them into folder directly
  if nent <= maxmenu:
    for i in range(len(nodes))[0::npp]:
      i1 = min(i+npp,len(nodes));
      pg = folder.page("%s - %s"%(nodes[i].name,nodes[i1-1].name));
      for node in nodes[i:i1]:
        pg.add(node);
  else:
    sublevels = int(math.ceil(math.log(nent)/math.log(maxmenu)))-1;
    nodes_per_sublevel = (maxmenu**sublevels)*npp;
#    print nent,sublevels,nodes_per_sublevel;
    for i in range(len(nodes))[0::nodes_per_sublevel]:
      i1 = min(i+nodes_per_sublevel,len(nodes));
      subfolder = folder.subfolder("%s - %s"%(nodes[i].name,nodes[i1-1].name));
      _make_folder(subfolder,nodes[i:i1],npp,maxmenu);

def make_node_folder (name,nodes,sorted=False,ncol=2,nrow=2,folder=None,maxmenu=25):
  """Creates a sub-folder with bookmarks for a group of nodes.
  """;
  # sort nodes by name
  if not sorted:
    nodes = list(nodes);
    global _int_or_str;
    from past.builtins import cmp
    from functools import cmp_to_key
    nodes.sort(key=cmp_to_key(lambda a,b:
      cmp(list(map(_int_or_str,a.name.split(':'))),list(map(_int_or_str,b.name.split(':'))))));
  # place in top-level folder or in subfolder
  if folder is None:
    folder = Folder(name);
  else:
    folder = folder.subfolder(name);
  # figure out # of plots per page
  npp = ncol*nrow;
  # call _make_folder() to do the real work
  _make_folder(folder,nodes,npp,maxmenu);
