# -*- coding: utf-8 -*-
# Meow.LSM
# Meow interface to LSM files

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
from Cattery.LSM.LSM import LSM as LSMClass
from Timba.utils import curry
import traceback
import Meow
import Meow.OptionTools
import Meow.Context
import math
from math import *

# constants for available LSM formats
# these are also used as labels in the GUI
NATIVE = "native [.lsm]";
NEWSTAR = "NEWSTAR [.mdl]";
NVSS = "NVSS file";
CLEAN = "clean components";
TEXT_RAD  = "text file (rad)";
TEXT_DMS  = "text file (hms/dms)";
VIZIER = "VizieR file";
OR_GSM = "OR_GSM file";
SKA = "SKA model catalog file";

class MeowLSM (object):
  def __init__ (self,filename=None,format=NATIVE,include_options=True,option_namespace='lsm'):
    """Initializes a MeowLSM object.
    A filename and a format may be specified, although the actual file will\
    only be loaded on demand.
    If include_options=True, immediately instantiates the options. If False, it is up to
    the caller to include the options in his menus.
    """;
    self.tdloption_namespace = option_namespace;
    self._compile_opts = [];
    self._runtime_opts = [];
    self.filename = filename;
    self.format = format;
    self.lsm = None;
    self.show_gui = False;
    # immediately include options, if needed
    if include_options:
      TDLCompileOptions(*self.compile_options());
      TDLRuntimeOptions(*self.runtime_options());

  def compile_options (self):
    """Returns list of compile-time options""";
    if not self._compile_opts:
      self._compile_opts.append(
        TDLOption("filename","LSM file",
                   TDLFileSelect("*.lsm *.txt *.*",default=self.filename,exist=True),
                   namespace=self)
      );
      format_opt = TDLOption("format","File format",
                 [ NATIVE,NEWSTAR,NVSS,CLEAN,TEXT_RAD,TEXT_DMS,VIZIER,OR_GSM,SKA ],
                 namespace=self);
      self._compile_opts.append(format_opt);
      beam_opt = TDLOption("beam_expr","Primary beam expression",
        [ None,
          "cos(min(65*fq*r,1.0881))**6   # WSRT",
          "max(cos(65*fq*r)**6,.01)      # NEWSTAR",
        ],more=str,namespace=self,
        doc="""<P>This is the primary power beam expression used to convert intrinsic <i>I</i> flux 
        as given by the LSM to apparent flux <b>for purposes of sorting only</b>.
        If not specified, LSM sources will be sorted by intrinsic <i>I</i> flux. The order of the sources is 
        relevant to the "LSM subset" option below.</P>
        
        <P>Note again that this beam expression only affects the sort order of the sources. The actual 
        source fluxes used in your tree will always be as per the LSM itself -- it is up to the tree to
        incorporate a proper primary beam model as required.</P>
        
        <P>The beam expression can be any valid Python expression (all functions from the 
        math module are available). The variables 'r' and 'fq' refer to the source distance 
        from field centre (in radians) and reference frequency (in GHz), respectively.</P>
        
        <P>If your LSM uses apparent rather than intrinsic fluxes, or you don't need to sort
        the sources by apparent flux, set this option to 'None'.</P>""");
      self._compile_opts.append(beam_opt);

      subset_doc = """<P>Selects a subset of sources from the
          LSM. You may specify sources:</P>

          <P>* by name</P>
          
          <P>* by ordinal number (from 0, in order of decreasing brightness)</P>

          <P>* as ranges, e.g. "M:N" (M to N inclusive), or ":M" (0 to M), or "N:" (N to last). 

          <P>Prefixing a name or a number or a range by "-" excludes that source from the selection.</P>

          <P>Examples: "-foo" (all sources except 'foo'), "0 3:5 bar" (sources #0,3,4,5 and 'bar'),
          "-3:5 -foo" (all sources except 3,4,5 and 'foo'.)</P>""";

      subset_opt = TDLOption('lsm_subset',"Use subset of LSM sources",
          ["all"],more=str,namespace=self,doc=subset_doc);
      self._subset_parser = Meow.OptionTools.ListOptionParser(minval=0,name="source");
      subset_opt.set_validator(self._subset_parser.validator);
      self._compile_opts.append(subset_opt);
      solve_subset_opt = TDLOption("solve_subset","For which sources",["all"],
            more=str,namespace=self,doc=subset_doc);
      self._solve_subset_parser = Meow.OptionTools.ListOptionParser(minval=0,name="source");
      self._compile_opts.append(
        TDLMenu("Make solvable source parameters",
          solve_subset_opt,
          TDLOption("solve_I","I",False,namespace=self),
          TDLOption("solve_Q","Q",False,namespace=self),
          TDLOption("solve_U","U",False,namespace=self),
          TDLOption("solve_V","V",False,namespace=self),
          TDLOption("solve_spi","spectral index",False,namespace=self),
          TDLOption("solve_pos","position",False,namespace=self),
          TDLOption("solve_RM","rotation measure",False,namespace=self),
          TDLOption("solve_shape","shape (for extended sources)",False,namespace=self),
          toggle='solvable_sources',namespace=self,
        ));
      save_opt = TDLMenu("Save LSM in native format",toggle="save_native",default=False,namespace=self,
                          *( TDLOption("save_native_filename","Filename",
                                     TDLFileSelect("*.lsm",exist=False),
                             namespace=self),));
      self._compile_opts.append(save_opt);
      save_txt_opt = TDLMenu("Save LSM in text (hms/dms) format",toggle="save_text",
                            default=False,namespace=self,
                          *( TDLOption("save_text","Filename",False,namespace=self),));
      self._compile_opts.append(save_txt_opt);
      
      def _select_format (format):
        if format == NATIVE:
          save_opt.set_value(False,save=False);
      format_opt.when_changed(_select_format);

 #     self._compile_opts += [
 #       TDLOption("show_gui","Show LSM GUI",False,namespace=self)
 #     ];
      
    return self._compile_opts;

  def runtime_options (self):
    """Makes and returns list of compile-time options""";
    # no runtime options, for now
    return self._runtime_opts;

  def load (self,ns,filename=None,format=None):
    """Loads LSM file. Filename/format arguments may be used to override those
    specified in the constructor or via options""";
    filename = filename or self.filename;
    format = format or self.format;

    self.lsm = LSMClass();

    # set up table of format readers
    # all are expected to take an lsm object (e.g. self) as arg 1,
    # a filename as arg2, and a node scope as arg3
    FORMAT_READERS = {};
    FORMAT_READERS[NATIVE]   = LSM.load;
    FORMAT_READERS[NEWSTAR]  = LSM.build_from_newstar;
    FORMAT_READERS[NVSS]     = LSM.build_from_catalog;
    FORMAT_READERS[CLEAN]    = LSM.build_from_complist;
    FORMAT_READERS[TEXT_RAD] = LSM.build_from_extlist_rad;
    FORMAT_READERS[TEXT_DMS] = LSM.build_from_extlist;
    FORMAT_READERS[VIZIER]   = LSM.build_from_vizier;
    FORMAT_READERS[OR_GSM]      = LSM.build_from_orgsm;
    FORMAT_READERS[SKA]      = LSM.build_from_ska;

    # read LSM using the selected format reader
    reader = FORMAT_READERS.get(format,None);
    if reader is None:
      raise TypeError,"Unknown LSM format '%s'"%format;

    try:
      reader(self.lsm,filename,ns);
    except:
      traceback.print_exc();
      raise RuntimeError,"""There was an error reading the LSM file %s. Either the file format is set incorrectly, or the LSM file is corrupt."""%filename;

    # save if needed
    if self.save_native and self.save_native_filename:
      self.lsm.save(self.save_native_filename);
    # save if needed
    if self.save_text and self.save_text_filename:
      self.lsm.save_as_extlist(self.save_text_filename,ns,prefix='');
      
#    if self.show_gui:
#      self.lsm.display()

  def source_list (self,ns,max_sources=None,**kw):
    """Reads LSM and returns a list of Meow objects.
    ns is node scope in which they will be created.
    Keyword arguments may be used to indicate which of the source attributes are to be
    created as Parms, use e.g. I=Meow.Parm(tags="flux") for this.
    The use_parms option may override this.
    """;
    if self.filename is None:
      return [];
    if self.lsm is None:
      self.load(ns);
    # all=1 returns unsorted list, so use a large count instead, to get a sorted list
    plist = self.lsm.queryLSM(count=9999999);
    
    # parse the beam expression
    if self.beam_expr is not None:
      try:
        beam_func = eval("lambda r,fq:"+self.beam_expr);
      except:
        raise RuntimeError,"invalid beam expression";
    else:
      beam_func = None;
    
    # make list of direction,punit,I,I_apparent tuples
    parm = Meow.Parm(tags="source solvable");
    srclist = [];
    for pu in plist:
      ra,dec,I,Q,U,V,spi,freq0,RM = pu.getEssentialParms(ns);
      if self.solve_pos:
        ra = parm.new(ra);
        dec = parm.new(dec);
      direction = Meow.Direction(ns,pu.name,ra,dec,static=not self.solve_pos);
      Iapp = I;
      if beam_func is not None:
      # if phase centre is already set (i.e. static), then lmn will be computed here, and we
      # can apply a beam expression
        lmn = direction.lmn_static();
        if lmn is not None:
          r = sqrt(lmn[0]**2+lmn[1]**2);
          Iapp = I*beam_func(r,freq0*1e-9 or 1.4);  # use 1.4 GHz if ref frequency not specified
      # append to list
      srclist.append((pu.name,direction,pu,I,Iapp));
    # sort list by decreasing apparent flux
    srclist.sort(lambda a,b:cmp(b[4],a[4]));
    
    srclist_full = srclist;
    # extract active subset
    srclist = self._subset_parser.apply(self.lsm_subset,srclist_full,names=[src[0] for src in srclist_full]);
    # extract solvable subset
    solve_subset = self._subset_parser.apply(self.solve_subset,srclist_full,names=[src[0] for src in srclist_full]);
    solve_subset = set([src[0] for src in solve_subset]);

    # make copy of kw dict to be used for sources not in solvable set
    kw_nonsolve = dict(kw);
    # and update kw dict to be used for sources in solvable set
    if self.solvable_sources:
      if self.solve_I:
        kw.setdefault("I",parm);
      if self.solve_Q:
        kw.setdefault("Q",parm);
      if self.solve_U:
        kw.setdefault("U",parm);
      if self.solve_V:
        kw.setdefault("V",parm);
      if self.solve_spi:
        kw.setdefault("spi",parm);
      if self.solve_RM:
        kw.setdefault("RM",parm);
      if self.solve_pos:
        kw.setdefault("ra",parm);
        kw.setdefault("dec",parm);
      if self.solve_shape:
        kw.setdefault("sx",parm);
        kw.setdefault("sy",parm);
        kw.setdefault("phi",parm);

    # make Meow list
    source_model = []

  ## Note: conversion from AIPS++ componentlist Gaussians to Gaussian Nodes
  ### eX, eY : multiply by 2
  ### eP: change sign
    for name,direction,pu,I,Iapp in srclist:
#      print "%-20s %12f %12f"%(pu.name,I,Iapp);
      src = {};
      ( src['ra'],src['dec'],
        src['I'],src['Q'],src['U'],src['V'],
        src['spi'],src['freq0'],src['RM']    ) = pu.getEssentialParms(ns)
      (eX,eY,eP) = pu.getExtParms()
      # scale 2 difference
      src['sx'] = eX*2
      src['sy'] = eY*2
      src['phi'] = -eP
      # override zero values with None so that Meow can make smaller trees
      if not src['RM']:
        src['RM'] = None;
      if not src['spi']:
        src['spi'] = None;
        if src['RM'] is None:
          src['freq0'] = None;
      ## construct parms or constants for source attributes
      ## if source is in solvable set (solvable_source_set of None means all are solvable),
      ## use the kw dict, else use the nonsolve dict for source parameters
      if name in solve_subset:
        solvable = True;
        kwdict = kw;
      else:
        solvable = False;
        kwdict = kw_nonsolve;
      for key,value in src.iteritems():
        meowparm = kwdict.get(key);
        if isinstance(meowparm,Meow.Parm):
          src[key] = meowparm.new(value);
        elif meowparm is not None:
          src[key] = value;

      if eX or eY or eP:
        # Gaussians
        if eY:
          size,phi = [src['sx'],src['sy']],src['phi'];
        else:
          size,phi = src['sx'],None;
        src = Meow.GaussianSource(ns,name=pu.name,
                I=src['I'],Q=src['Q'],U=src['U'],V=src['V'],
                direction=direction,
                spi=src['spi'],freq0=src['freq0'],RM=src['RM'],
                size=size,phi=phi);
      else:
        src = Meow.PointSource(ns,name=pu.name,
                I=src['I'],Q=src['Q'],U=src['U'],V=src['V'],
                direction=direction,
                spi=src['spi'],freq0=src['freq0'],RM=src['RM']);
                
      # check for beam LM
      if pu._lm is not None:
        src.set_attr('beam_lm',pu._lm);
              
      src.solvable = solvable;
      src.set_attr('Iapp',Iapp);
      source_model.append(src);
      
    return source_model;

