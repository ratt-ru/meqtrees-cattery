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

from Timba.TDL import *
from Timba.array import array
import Context

import re
import math

_wsrt_list = [ str(i) for i in range(10)] + ['A','B','C','D','E','F'];
_vla_list = [ str(i) for i in range(1,28) ];

_uvw_from_ms = "from MS";
_uvw_compute_mirror = "compute (VLA convention)";
_uvw_compute = "compute (WSRT convention)";

uvw_source_opt = TDLOption('uvw_source',"UVW coordinates",
      [_uvw_from_ms,_uvw_compute_mirror,_uvw_compute],
      doc="""UVW coordinates can be read from the MS, or recomputed on the fly.
      In the latter case, you have a choice of two opposite sign conventions.""");
uvw_refant_opt = TDLOption('uvw_refant',"Reference antenna #",0,more=int,
      doc="""This is the reference antenna used to compute antenna-based UVWs.
      Specify a 0-based antenna index. Note that this antenna must be present in
      all timeslots where data is available.""");

_options = [ uvw_source_opt,uvw_refant_opt ];
uvw_source_opt.when_changed(lambda x:uvw_refant_opt.show(x==_uvw_from_ms));


class IfrArray (object):
  def compile_options ():
    return _options;
  compile_options = staticmethod(compile_options);

  # syntax of IfrSet.subset() and IfrArray.subset() strings
  ifr_spec_syntax = """<P>Specify interferometers as a list of tokens, separated by commas or spaces. 
    The format of each token is as follows:</P>

    <P>"*" or "all": selects all ifrs.</P>

    <P>"&lt;=N", "&lt;N", "&gt;=N", "&gt;N": selects by baseline length.</P>

    <P>"P-Q" or "P:Q" or "PQ": selects interferometer P-Q. (The third form is useful when 
    station names are a single character, as for WSRT.) Q may be "*", in which case all 
    interferometers formed with antenna P are selected.</P>

    <P>"-xxx": excludes baselines specified by token xxx. E.g. "all -4* -&lt;144 0-1 AB" selects all 
    baselines, then throws out those with antenna 4 and those shorter than 144m, then adds 0-1 
    and A-B. Note also that as a shortcut, if the first token in the list is negated, then a "*" 
    is implicitly prepended. That is, a "-AB" at the start is equivalent to "all -AB".</P>
    """;

  class IfrSet (object):
    """IfrSet is an object representing a set of interferometers.
    """
    def __init__ (self,station_list=None,ifr_index=None,positions=None,observatory=None):
      """Creates a set of interferometers from either 
      * 'station_list': a list of station IDs
      * 'station_list': a list of (ip,p): station numbers and ID tuples.
      * 'ifr_index': a list of (ip,p),(iq,q) pairs,
      These parameters are optional. They're used in subset() selection strings:
        'positions' is a [nant,3] array of antenna positions, if known.
        'observatory' is an optional observatory name, if known.
      """
      self.observatory = observatory;
      if station_list:
        # make _station_index: list of (ip,p) pairs.
        if isinstance(station_list[0],(list,tuple)):
          self._station_index = list(station_list);
        else:
          self._station_index = list(enumerate(station_list));
        # _ifr_index is a list of all IFRs, as (ip,p),(iq,q) pairs.
        self._ifr_index = [ ((ip,p),(iq,q)) for ip,p in self._station_index
                                            for iq,q in self._station_index if ip<iq ];
      elif ifr_index:
        self._ifr_index = list(sorted(ifr_index));
        # rebuild station index from that
        stations = set([ ifrs[0] for ifrs in self._ifr_index ]);
        stations.update([ ifrs[1] for ifrs in self._ifr_index ]);
        self._station_index = list(sorted(stations));
      # _number_to_name: maps ip to p (since ip's are not necessarily contiguous)
      self._number_to_name = dict(self._station_index);
      # make reverse dictionary: from name (also in uppercase) to number.
      self._name_to_number = dict([(p.upper(),ip) for ip,p in self._station_index]);
      self._name_to_number.update([(p,ip) for ip,p in self._station_index]);
      # make a dict of positions and baseline lengths
      # _positions[ip] gives the position of antenna #ip
      # _baselines[ip,iq] gives the baseline length for ip,iq
      self.positions = positions;
      if positions is not None:
        self._baselines = dict([
          ((ip,iq),math.sqrt(((positions[ip,:]-positions[iq,:])**2).sum()))
          for ip,p in self._station_index for iq,q in self._station_index
        ]);
      else:
        self._baselines = None;
      # make a few more useful lists
      # _stations is just a list of all station anmes
      self._stations = [ p for ip,p in self._station_index ];
      # _ifrs is a list of all IFR names, as (p,q) pairs
      self._ifrs = [ (px[1],qx[1]) for px,qx in self._ifr_index ];
      # fill in dictionary of aliases for baseline specifications
      self._ifr_spec_aliases = dict();
      # '*' or 'all' selects all baselines
      self._ifr_spec_aliases['*'] = self._ifr_spec_aliases['ALL'] = set(self._ifr_index);
      # form up doc string
      self.subset_doc = IfrArray.ifr_spec_syntax;
      # special aliases for WSRT
      if observatory and observatory.upper() == "WSRT":
        fixed = set([ (ip,p) for ip,p in self._station_index if ip < 10 ]);
        movable = set([ (ip,p) for ip,p in self._station_index if ip >= 10 ]);
        self._ifr_spec_aliases["FF"] = self._ifr_spec_aliases["F-F"] = \
          set([(px,qx) for px,qx in self._ifr_index if px in fixed and qx in fixed]);
        self._ifr_spec_aliases["FM"] = self._ifr_spec_aliases["F-M"] = \
          set([(px,qx) for px,qx in self._ifr_index if px in fixed and qx in movable]);
        self._ifr_spec_aliases["MM"] = self._ifr_spec_aliases["M-M"] = \
          set([(px,qx) for px,qx in self._ifr_index if px in movable and qx in movable]);
        self.subset_doc += """
         <P>"FF", "FM", "MM" (or "F-F", "F-M", "M-M"): selects the stadard WSRT sets
         of fixed-fixed, fixed-movable and movable-movable baselines. Can also be negated,
         so "* -MM" selects all baselines except the movable-movable ones.
         </P>
        """;

    def stations (self):
      return self._stations;

    def ifrs (self):
      return self._ifrs;

    def station_index (self):
      return self._station_index;

    def ifr_index (self):
      return self._ifr_index;

    def number_of_station (self,name):
      """Returns the ordinal number (i.e. ANTENNA_ID) of the named station.""";
      return self._name_to_number[name];

    def station_position (self,ip):
      return self.positions[ip,:] if self.positions is not None else array((0,0,0));

    _comparison_predicates = {
      '<':  (lambda a,b: a<b),
      '<=': (lambda a,b: a<=b),
      '>':  (lambda a,b: a>b),
      '>=': (lambda a,b: a>=b),
    };

    def subset (self,specstr,strict=False):
      """Parses an interferometer subset specification string, and returns an index of 
      interferometers found in the string. If strict is False, ignores any errors in the
      string, if true then throws an exception on error.
      
      %s
      
      Returns an IFR set object
      """%IfrArray.ifr_spec_syntax;
      result = set();
      for ispec,spec in enumerate(re.split("[\s,]+",specstr.strip())):
        # default action is to add to set, but a "-" prefix negates this
        add_or_remove = result.update;
        if spec.startswith("-"):
          # special case: "-" at the start means we begin with all baselines selected
          if not ispec:
            result = set(self._ifr_index);
          add_or_remove = result.difference_update;
          spec = spec[1:];
        elif spec.startswith("+"):
          spec = spec[1:];
        # check for named aliases
        named_set = self._ifr_spec_aliases.get(spec.upper(),None);
        if named_set:  
          add_or_remove(named_set);
          continue;
        # now check for baseline length specification
        match = re.match("^(<|<=|>|>=)([^=]+)$",spec);
        if match:
          if self._baselines:
            try:
              length = float(match.group(2));
            except:
              if strict:
                raise ValueError,"invalid ifr specification '%s'"%spec;
              print "Ignoring invalid ifr specification '%s'"%spec;
              continue;
            pred = self._comparison_predicates[match.group(1)];
            add_or_remove([(px,qx) for px,qx in self._ifr_index if pred(self._baselines[px[0],qx[0]],length)]);
          elif strict:
            raise ValueError,"can't use ifr specification %s: baseline information not available"%spec;
          else:
            print "Ignoring ifr specification %s: baseline information not available"%spec;
          continue;
        # else spit into antenna pair
        if len(spec) == 2:
          p,q = spec;
        else:
          p,q = re.split("[-:]",spec,1);
        ip = self._name_to_number.get(p.upper(),None);
        iq = self._name_to_number.get(q.upper(),None);
        # first token must be a valid antenna
        if ip is None:
          if strict:
            raise ValueError,"invalid ifr specification '%s' (station '%s' not known)"%(spec,p);
          print "Ignoring invalid ifr specification '%s' (station '%s' not known)"%(spec,p);
          continue;
        # second token may be a wildcard
        if q == "*":
          add_or_remove([(px,qx) for px,qx in self._ifr_index if px[0]==ip or qx[0]==ip ]);
        elif iq is None:
          if strict:
            raise ValueError,"invalid ifr specification '%s' (station '%s' not known)"%(spec,q);
          print "Ignoring invalid ifr specification '%s' (station '%s' not known)"%(spec,q);
        elif ip<iq:
          add_or_remove([(px,qx) for px,qx in self._ifr_index if (px[0],qx[0])==(ip,iq)]);
        elif ip>iq:
          add_or_remove([(px,qx) for px,qx in self._ifr_index if (px[0],qx[0])==(iq,ip)]);

      return IfrArray.IfrSet(ifr_index=result,observatory=self.observatory,positions=self.positions);

  def __init__(self,ns,station_list,station_index=None,uvw_table=None,
               observatory=None,
               ms_uvw=None,mirror_uvw=None,resamplers=False,positions=None):
    """Creates an IfrArray object, representing an interferometer array.
    'station_list' is a list of station IDs, not necessarily numeric.
        It can also be a list of (index,ID) tuples, when a subset of antennas
        is specified.
    'station_index' is an optional list of numeric station indices. If not given,
      [0,...,N-1] will be used. This is an alternative to passing in tuples
      to station_list.
    'uvw_table' is a path to a MEP table containing station UVWs.
    'ms_uvw' if True, causes UVWs to be read from the MS. If False, causes UVWs to
      be computed with a Meq.UVW node. If None, uses the global uvw_source TDLOption.
    'mirror_uvw' only applicable if UVWs are being computed. If True, uses the VLA
      UVW sign definition, if False, uses the WSRT one.
    'resamplers': if True (and ms_uvw=True), Resampler nodes will be put on the UVWs.
      This is useful when the tree is being executed under a ModRes.
    'positions' is a [nant,3] array of antenna positions, if known.
    'observatory': observatory identifier, if known
    """;
    self.ns = ns;
    self.observatory = observatory;
    # select UVW options
    if ms_uvw is None:
      ms_uvw = (uvw_source == _uvw_from_ms);
      if not ms_uvw:
        mirror_uvw = (uvw_source == _uvw_compute_mirror);
    # make list of station pairs: (0,p0),(1,p1),... etc.
    if isinstance(station_list,IfrArray.IfrSet):
      self.ifrset = station_list;
    else:
      if station_index:
        stataion_list = zip(station_index,station_list);
      self.ifrset = IfrArray.IfrSet(station_list,positions=positions,observatory=observatory);
    # expose some methods of IfrSet directly via our object
    for method in ('stations','ifrs','station_index','ifr_index','subset',
                   'station_position','number_of_station'):
      setattr(self,method,getattr(self.ifrset,method));
    # other init
    self._uvw_table = uvw_table;
    self._ms_uvw = ms_uvw;
    self._mirror_uvw = mirror_uvw;
    self._resamplers = resamplers;
    self._jones = [];

  @staticmethod
  def WSRT (ns,stations=14,uvw_table=None,mirror_uvw=False):
    """Creates and returns an IfrArray for WSRT, i.e., with proper labels.
    The 'stations' argument can be either a number of stations (then only
    the first N antennas will be used), or a list of station indices.
    If not given, then the full WSRT array is used.
    """;
    if isinstance(stations,int):
      stations = enumerate(_wsrt_list[:stations]);
    elif isinstance(stations,(list,tuple)):
      stations = [ (i,_wsrt_list[i]) for i in stations ];
    else:
      raise TypeError,"WSRT 'stations' argument must be a list of stations, or a number";
    return IfrArray(ns,stations,uvw_table=uvw_table,mirror_uvw=mirror_uvw,
                    observatory="WSRT");

  @staticmethod
  def VLA (ns,stations=27,uvw_table=None,mirror_uvw=False):
    """Creates and returns an IfrArray for VLA, i.e., with proper labels.
    The 'stations' argument can be either a number of stations (then only
    the first N antennas will be used), or a list of station indices.
    If not given, then the full WSRT array is used.
    """;
    if isinstance(stations,int):
      stations = enumerate(_vla_list[:stations]);
    elif isinstance(stations,(list,tuple)):
      stations = [ (i,_vla_list[i]) for i in stations ];
    else:
      raise TypeError,"VLA 'stations' argument must be a list of stations, or a number";
    return IfrArray(ns,stations,station_index=index,uvw_table=uvw_table,mirror_uvw=mirror_uvw,
                    observatory="VLA");

  def num_stations (self):
    return len(self.stations());

  def num_ifrs (self):
    return len(self.ifrs());

  ifr_spec_syntax = """<P>Specify interferometers as a list of tokens, separated by commas or spaces. 
  The format of each token is as follows:</P>

  <P>"*" or "all": selects all ifrs.</P>

  <P>"&lt;=N", "&lt;N", "&gt;=N", "&gt;N": selects by baseline length.</P>

  <P>"P-Q" or "P:Q" or "PQ": selects interferometer P-Q. (The third form is useful when 
  station names are a single character, as for WSRT.) Q may be "*", in which case all 
  interferometers formed with antenna P are selected.</P>

  <P>"-xxx": excludes baselines specified by token xxx. E.g. "all -4* -&lt;144 0-1 AB" selects all 
  baselines, then throws out those with antenna 4 and those shorter than 144m, then adds 0-1 
  and A-B. Note also that as a shortcut, if the first token in the list is negated, then a "*" 
  is implicitly prepended. That is, a "-AB" at the start is equivalent to "all -AB".</P>
  """;


  def spigots (self,node=None,corr=[0,1,2,3],column="DATA",**kw):
    """Creates (if necessary) and returns data spigots, as an unqualified node.
    Extra keyword arguments will be passed to Spigot node."""
    if node is None:
      node = self.ns.spigot if column == "DATA" else self.ns.spigot(column.lower());
    (ip0,p0),(iq0,q0) = self.ifr_index()[0];
    if len(corr) == 4:
      dims = [2,2];
    else:
      dims = [ len(corr) ];
    if not node(p0,q0).initialized():
      for (ip,p),(iq,q) in self.ifr_index():
        node(p,q) << Meq.Spigot(station_1_index=ip,station_2_index=iq,input_col=column,
                                corr_index=corr,dims=dims,**kw);
    return node;

  def sinks (self,children,node=None,**kw):
    """Creates (if necessary) and returns sinks, as an unqualified node.
    The 'children' argument should be a list of child nodes, which
    will be qualified with an interferometer pair.
    Extra keyword arguments will be passed to Sink node."""
    if node is None:
      node = self.ns.sink;
    (ip0,p0),(iq0,q0) = self.ifr_index()[0];
    if not node(p0,q0).initialized():
      for (ip,p),(iq,q) in self.ifr_index():
        node(p,q) << Meq.Sink(children=children(p,q),
                               station_1_index=ip,station_2_index=iq,
                               **kw);
    return node;

  def xyz0 (self):
    """Returns array reference position node""";
    self.xyz();
    return self.ns.xyz0;

  def xyz (self,*quals):
    """Returns unqualified station position nodes,
    If a station is supplied, returns XYZ node for that station""";
    xyz0 = self.ns.xyz0;
    #print "xyz",  xyz0.initialized(),self.station_index();
    if not xyz0.initialized():
      for (ip,p) in self.station_index():
        # since the Meow.ReadVisHeader script knows nothing about
        # our station labels, the x/y/z nodes themselves are
        # indexed with station _numbers_ instead.
        # to avoid confusion, we call them "x:num0", etc.
        num = 'num'+str(ip);
        # create XYZ nodes
        x,y,z = self.station_position(ip);
        xyz = self.ns.xyz(p) << Meq.Composer(
          self.ns.x(num) << x,
          self.ns.y(num) << y,
          self.ns.z(num) << z
        );
        if not xyz0.initialized():
          xyz0 << Meq.Selector(xyz); # xyz0 == xyz first station essentially
    return self.ns.xyz(*quals);

  def x_station(self,*quals):
    self.xyz();
    if not self.ns.x_station(self.stations()[0]).initialized():
      for station in  self.stations():
        self.ns.x_station(station) << Meq.Selector(self.ns.xyz(station),index=0);
    return self.ns.x_station(*quals);

  def y_station(self,*quals):
    self.xyz();
    if not self.ns.y_station(self.stations()[0]).initialized():
      for station in  self.stations():
        self.ns.y_station(station) << Meq.Selector(self.ns.xyz(station),index=1);
    return self.ns.y_station(*quals);

  def z_station(self,*quals):
    self.xyz();
    if not self.ns.z_station(self.stations()[0]).initialized():
      for station in  self.stations():
        self.ns.z_station(station) << Meq.Selector(self.ns.xyz(station),index=2);
    return self.ns.z_station(*quals);

  def uvw (self,dir0=None,*quals):
    """returns station UVW node(s) for a given phase centre direction,
    or using the global phase center if None is given.
    For the global phase center (dir0=None), can also use UVWs from the MS or from
    MEP tables, according to whatever was specified in the constructor.
    For other directions, UVWs are always computed.
    If a station is supplied, returns UVW node for that station""";
    if dir0 is not None:
      radec0 = Context.get_dir0(dir0).radec();
      uvw = self.ns.uvw.qadd(radec0);
    else:
      uvw = self.ns.uvw;
    if not uvw(self.stations()[0]).initialized():
      if self._ms_uvw:
        # read UVWs from MS
        # first station gets (0,0,0), the rest is via subtraction
        ip0,p0 = self.station_index()[uvw_refant];
        uvw(p0) << Meq.Composer(0,0,0);
        uvw(p0)._multiproc = True; # hint to parallelizer to clone this node on all processors
        for iq,q in self.station_index():
          if iq < ip0:
            m_uvw(q) << Meq.Spigot(station_1_index=iq,station_2_index=ip0,input_col='UVW');
            spigdef = -m_uvw(q);
          elif iq > ip0:
            spigdef = Meq.Spigot(station_1_index=ip0,station_2_index=iq,input_col='UVW');
          else:
            continue;
          if self._resamplers:
            uvw(q) << Meq.Resampler(uvw(q,"res") << spigdef,mode=1);
          else:
            uvw(q) << spigdef;
      elif self._uvw_table:
        # read UVWs from MEP table
        for station in self.stations():
          uvw(station) << Meq.Composer(
            self.ns.u.qadd(radec0)(station) << Meq.Parm(table_name=self._uvw_table),
            self.ns.v.qadd(radec0)(station) << Meq.Parm(table_name=self._uvw_table),
            self.ns.w.qadd(radec0)(station) << Meq.Parm(table_name=self._uvw_table)
          );
      else:
        # compute UVWs via a UVW node
        radec0 = Context.get_dir0(dir0).radec();
        xyz0 = self.xyz0();
        xyz = self.xyz();
        for station in self.stations():
          uvw_def = Meq.UVW(radec = radec0,
                            xyz_0 = xyz0,
                            xyz   = xyz(station));
          if self._mirror_uvw:
            uvw(station) << Meq.Negate(self.ns.m_uvw(station) << uvw_def );
          else:
            uvw(station) << uvw_def;
    return uvw(*quals);

  def uvw_ifr (self,dir0,*quals):
    """returns interferometer UVW node(s) for a given phase centre direction,
    or using the global phase center if None is given.
    If an IFR is supplied, returns UVW node for that IFR""";
    dir0 = Context.get_dir0(dir0);
    radec0 = dir0.radec();
    uvw_ifr = self.ns.uvw_ifr.qadd(radec0);
    if not uvw_ifr(*(self.ifrs()[0])).initialized():
      uvw = self.uvw(dir0);
      for sta1,sta2 in self.ifrs():
        uvw_ifr(sta1,sta2) << uvw(sta2) - uvw(sta1);
    return uvw_ifr(*quals);

  def uv (self,dir0,*quals):
    """returns station UV node(s) for a given phase centre direction,
    or using the global phase center if None is given.
    If a station is supplied, returns UV node for that station""";
    dir0 = Context.get_dir0(dir0);
    radec0 = dir0.radec();
    uv = self.ns.uv.qadd(radec0);
    if not uv(self.stations()[0]).initialized():
      uvw = self.uvw(dir0);
      for station in self.stations():
        uv(station) << Meq.Selector(uvw(station),index=(0,1),multi=True);
    return uv(*quals);

  def uv_ifr (self,dir0,*quals):
    """returns interferometer UV node(s) for a given phase centre direction.
    or using the global phase center if None is given.
    If an IFR is supplied, returns UVW node for that IFR""";
    dir0 = Context.get_dir0(dir0);
    radec0 = dir0.radec();
    uv_ifr = self.ns.uv_ifr.qadd(radec0);
    if not uv_ifr(*(self.ifrs()[0])).initialized():
      uvw_ifr = self.uvw_ifr(dir0);
      for ifr in self.ifrs():
        uv_ifr(*ifr) << Meq.Selector(uvw_ifr(*ifr),index=(0,1),multi=True);
    return uv_ifr(*quals);
