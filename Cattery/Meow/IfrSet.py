# -*- coding: utf-8 -*-
#
#% $Id: IfrArray.py 7467 2010-01-15 22:57:48Z oms $
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
import re
import math
import fnmatch
import traceback

from Timba.array import array

ifr_spec_syntax = """<P>Specify interferometers as a list of tokens, separated by commas or spaces.
  The format of each token is as follows:</P>

  <P>"*" or "all" or "ALL": selects all ifrs.</P>

  <P>"&lt;=N" or ".le.N", "&lt;N" or ".lt.N", "&gt;=N" or ".ge.N", "&gt;N" or ".gt.N": selects by baseline length.</P>

  <P>"P-Q" or "P:Q" or "PQ": selects interferometer P-Q. (The third form is useful when
  station names are a single character, as for WSRT.) Q may be "*", in which case all
  interferometers formed with antenna P are selected.</P>

  <P>"-xxx": excludes baselines specified by token xxx. E.g. "all -4* -&lt;144 0-1 AB" selects all
  baselines, then throws out those with antenna 4 and those shorter than 144m, then adds 0-1
  and A-B. Note also that as a shortcut, if the first token in the list is negated, then a "*"
  is implicitly prepended. That is, a leading "-AB" is equivalent to "all -AB".</P>
  """;


class IfrSet (object):
  """IfrSet is an object representing a set of interferometers.
  """
  ifr_spec_syntax = ifr_spec_syntax;
  def __init__ (self,station_list=None,ifr_index=None,
                positions=None,observatory=None,label_sep=None,parent=None):
    """Creates a set of interferometers from either
      * 'station_list': a list of station IDs
      * 'station_list': a list of (ip,p): station numbers and ID tuples.
      * 'ifr_index': a list of (ip,p),(iq,q) pairs,
      'positions' is a [nant,3] array of antenna positions, if known. If not supplied,
          the baseline() method will return None, and subset() selection on baseline
          length will not be available.
      'observatory' is an optional observatory name, if known. Some observatories
          define standard baseline sets, for which aliases will be recognized in
          subset() (e.g. "FM" for WSRT fixed-movable set.)
      'label_sep' is a separator string for IFR labels. If None, then the default
          is to use "-", or "" if all station names are single-character.
      'parent' is a parent IfrSet object. If set, then information such as the above
      (position, observatory, etc.) is copied over.
    """
    # init these from args or from parent IfrSet
    self.positions = positions = positions if not parent else parent.positions;
    self.observatory = observatory = observatory if not parent else parent.observatory;
    self.label_sep = label_sep = label_sep if not parent else parent.label_sep;
    # init station list
    if station_list is not None:
      # make _station_index: list of (ip,p) pairs.
      if isinstance(station_list[0],(list,tuple)):
        self._station_index = list(station_list);
      else:
        self._station_index = list(enumerate(station_list));
      # _ifr_index is a list of all IFRs, as (ip,p),(iq,q) pairs.
      self._ifr_index = [ ((ip,p),(iq,q)) for ip,p in self._station_index
                                          for iq,q in self._station_index if ip<iq ];
    elif ifr_index is not None:
      self._ifr_index = list(sorted(ifr_index));
      # rebuild station index from that
      stations = set([ ifrs[0] for ifrs in self._ifr_index ]);
      stations.update([ ifrs[1] for ifrs in self._ifr_index ]);
      self._station_index = list(sorted(stations));
    # _number_to_name: maps ip to p (since ip's are not necessarily contiguous)
    self._number_to_name = dict(self._station_index);
    # make reverse dictionary: from name (also in uppercase) to number.
    self._name_to_number = dict([(str(p).upper(),ip) for ip,p in self._station_index]);
    self._name_to_number.update([(p,ip) for ip,p in self._station_index]);
    # make a dict of positions and baseline lengths
    # _positions[ip] gives the position of antenna #ip
    # _baselines[ip,iq] gives the baseline length for ip,iq
    if self.positions is not None:
      self._baseline_vectors = dict([
        ((ip,iq),(self.positions[ip,:]-self.positions[iq,:]))
        for ip,p in self._station_index for iq,q in self._station_index
      ]);
      self._baselines = dict([
        ((ip,iq),math.sqrt((self._baseline_vectors[ip,iq]**2).sum()))
        for ip,p in self._station_index for iq,q in self._station_index
      ]);
    else:
      self._baselines = self._baseline_vectors = None;
    # make a few more useful lists
    # _stations is just a list of all station anmes
    self._stations = [ p for ip,p in self._station_index ];
    # _ifrs is a list of all IFR names, as (p,q) pairs
    self._ifrs = [ (px[1],qx[1]) for px,qx in self._ifr_index ];
    # _ifr_labels is a dict of ip,iq -> IFR labels
    # default label separator is "-", or "" for single-char stations
    if self.label_sep is None:
      self.label_sep = "-" if max(map(len,map(str,self._stations)))>1 else "";
    self._ifr_labels = dict([((ip,iq),"%s%s%s"%(p,self.label_sep,q)) for (ip,p),(iq,q) in self._ifr_index ]);
    # fill in dictionary of aliases for baseline specifications
    if parent:
      # if we're a subset of parent, make copy of parent dict (need copy since "*" item is different)
      self._ifr_spec_aliases = parent._ifr_spec_aliases;
      self.subset_doc = parent.subset_doc;
    else:
      self._ifr_spec_aliases = dict();
      # form up doc string
      self.subset_doc = ifr_spec_syntax;
      # special aliases for WSRT
      if self.observatory and self.observatory.upper() == "WSRT":
        fixed = set([ (ip,p) for ip,p in self._station_index if ip < 10 ]);
        movable = set([ (ip,p) for ip,p in self._station_index if ip >= 10 ]);
        self._ifr_spec_aliases["FF"] = self._ifr_spec_aliases["F-F"] = \
          set([(px,qx) for px,qx in self._ifr_index if px in fixed and qx in fixed]);
        self._ifr_spec_aliases["FM"] = self._ifr_spec_aliases["F-M"] = \
          set([(px,qx) for px,qx in self._ifr_index if px in fixed and qx in movable]);
        self._ifr_spec_aliases["MM"] = self._ifr_spec_aliases["M-M"] = \
          set([(px,qx) for px,qx in self._ifr_index if px in movable and qx in movable]);
        self._ifr_spec_aliases["S85"] = self.ifr_index_subset("-45 -56 -67 -9A -AB -CD");
        self._ifr_spec_aliases["S83"] = self.ifr_index_subset("-45 -46 -56 -67 -68 -9A -AB -CD");
        self.subset_doc += """
          <P>You appear to have an WSRT MS, so in addition to the above, the following standard
          designations are recognized:</P>

          <P>"FF", "FM", "MM" (or "F-F", "F-M", "M-M"): selects the standard WSRT sets
          of fixed-fixed, fixed-movable and movable-movable baselines. Can also be negated,
          so "* -MM" selects all baselines except the movable-movable ones.
          </P>

          <P>"S85": standard set of 85 baselines, i.e. all except the three shortest spacings
          (9A, AB, CD), and the oft-contaminated 45, 56 and 67.</P>

          <P>"S83": same as above, minus also 46 and 68.</P>
        """;
      # '*' or 'all' selects all baselines in the ifr set
      self._ifr_spec_aliases['*'] = self._ifr_spec_aliases['ALL'] = set(self._ifr_index);

  def stations (self):
    """Returns list of station names.""";
    return self._stations;

  def station_numbers (self):
    """Returns list of ordinal station numbers. These are not necessarily consecutive,
    as they refer to station numbers in the original station list (i.e. the ANTENNA
    subtable of an MS), of which we may be dealing with a subset."""
    return [ ip for ip,p in self._station_index ];

  def ifrs (self):
    """Returns list of ifr names (i.e. p,q pairs, where p and q are station names).""";
    return self._ifrs;

  def ifr_numbers (self):
    """Returns list of ifr numbers (i.e. ip,iq pairs, where ip and iq are station numbers).""";
    return [ (px[0],qx[0]) for px,qx in self._ifr_index ];

  def station_index (self):
    """Returns station index: a list of ip,p pairs, where ip is the station number,
    and p is the name.""";
    return self._station_index;

  def ifr_index (self):
    """Returns ifr index: a list of (ip,p),(iq,q) pairs.""";
    return self._ifr_index;

  def ifr_label (self,ip,iq):
    """Returns ifr label corresponding to IFR ip,iq. This is usually "p-q".""";
    return self._ifr_labels[ip,iq];

  def number_of_station (self,name):
    """Returns the ordinal number of the named station. Raises KeyError if no such name.""";
    return self._name_to_number[name];

  def station_position (self,ip):
    """Returns (x,y,z) position (as an array) of station number ip.""";
    return self.positions[ip,:] if self.positions is not None else array((0,0,0));

  def baseline (self,ip,iq):
    """Returns (x,y,z) baseline for ifr number ip,iq.""";
    return self._baselines[ip,iq] if self._baselines else 0.;

  def baseline_vector (self,ip,iq):
    """Returns (x,y,z) baseline for ifr number ip,iq.""";
    return self._baseline_vectors[ip,iq] if self._baselines else 0.;

  def taql_string (self):
    """Returns TaQL string representing this set of baselines.""";
    return "||".join(["(ANTENNA1==%d&&ANTENNA2==%d)"%(ip,iq) for (ip,p),(iq,q) in self._ifr_index ]);

  _epsilon = 0.5

  _comparison_predicates = {
    '<':  (lambda a,b: a<b-IfrSet._epsilon),
    '<=': (lambda a,b: a<b+IfrSet._epsilon),
    '>':  (lambda a,b: a>b+IfrSet._epsilon),
    '>=': (lambda a,b: a>b-IfrSet._epsilon),
    '=': (lambda a,b: abs(a-b)<IfrSet._epsilon),
    '.lt.':  (lambda a,b: a<b-IfrSet._epsilon),
    '.le.': (lambda a,b: a<b+IfrSet._epsilon),
    '.gt.':  (lambda a,b: a>b+IfrSet._epsilon),
    '.ge.': (lambda a,b: a>b-IfrSet._epsilon),
    '.eq.': (lambda a,b: abs(a-b)<IfrSet._epsilon),
  };

  def subset (self,specstr,strict=False):
    """Parses an interferometer subset specification string, and returns an new IfrSet
    object composed of interferometers found in the string. If strict is False,
    ignores any errors in the string, if true then throws an exception on error.

    %s

    Returns a new IfrSet object.
    """%ifr_spec_syntax;
    return IfrSet(ifr_index=self.ifr_index_subset(specstr,strict=strict),parent=self);

  def ifr_index_subset (self,specstr,strict=False):
    """Parses an interferometer subset specification string, and returns a set of ifr indices
    for interferometers found in the string. If strict is False, ignores any errors in the string,
    if true then throws an exception on error.

    %s

    Returns a set of ((ip,p),(iq,q)) tuples.
    """%ifr_spec_syntax;
    # if None, then return everything
    if not specstr:
      return set(self._ifr_index);
    result = set();
#    print specstr;
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
      match = re.match("^(<|<=|>|>=|=|\\.lt\\.|\\.le\\.|\\.ge\\.|\\.gt\\.|\\.eq\\.)([^=]+)$",spec,re.IGNORECASE);
#      print spec,match;
      if match:
        if self._baselines:
          try:
            length = float(match.group(2));
          except:
            if strict:
              raise ValueError,"invalid ifr specification '%s'"%spec;
            print "Ignoring invalid ifr specification '%s'"%spec;
            continue;
          pred = self._comparison_predicates[match.group(1).lower()];
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
        pq = re.split("[-:]",spec,1);
        if len(pq) == 2:
          p,q = pq;
        else:
          raise ValueError,"illegal interferometer subset: %s"%spec;
      # if first is wildcard, change them around, so X* becomes *X
      if p == "*":
        p,q = q,p;
      p = p.upper();
      q = q.upper();
      psubset = fnmatch.filter(self._stations,p);
      qsubset = fnmatch.filter(self._stations,q);
      print p,psubset,q,qsubset;
      for p in psubset:
        for q in qsubset:
          ip = self._name_to_number.get(p.upper(),None);
          iq = self._name_to_number.get(q.upper(),None);
          # first token must be a valid antenna
          if ip is None:
            if strict:
              raise ValueError,"invalid ifr specification '%s' (station '%s' not found)"%(spec,p);
            print "Ignoring invalid ifr specification '%s' (station '%s' not found)"%(spec,p);
    #        traceback.print_stack();
            continue;
          # second token may be a wildcard
          if q == "*":
            add_or_remove([(px,qx) for px,qx in self._ifr_index if px[0]==ip or qx[0]==ip ]);
          elif iq is None:
            if strict:
              raise ValueError,"invalid ifr specification '%s' (station '%s' not found)"%(spec,q);
            print "Ignoring invalid ifr specification '%s' (station '%s' not found)"%(spec,q);
            traceback.print_stack();
          elif ip<iq:
            add_or_remove([(px,qx) for px,qx in self._ifr_index if (px[0],qx[0])==(ip,iq)]);
          elif ip>iq:
            add_or_remove([(px,qx) for px,qx in self._ifr_index if (px[0],qx[0])==(iq,ip)]);
    # intersect final result with our own set of IFRs: if we're already a subset of some parent set, then the spec
    # string may refer to ifrs from our parent set that are not members of our set. These must be trimmed.
    return result.intersection(self._ifr_index);

def from_ms (ms):
  # import table class
  try:
    from pyrap_tables import table,tablecopy,tableexists,tabledelete
  except:
    try:
      from pyrap.tables import table,tablecopy,tableexists,tabledelete
    except:
      print "Failed to import pyrap_tables or pyrap.tables. Please install the pyrap "
      "package from http://code.google.com/p/pyrap/, or from a MeqTrees binary repository "
      "(see http://www.astron.nl/meqwiki/Downloading.)"
      raise;
  # open MS
  if isinstance(ms,str):
    ms = table(ms);
  elif not isinstance(ms,table):
    raise TypeError,"'ms' argument must be an MS name or a table object";

  # open ANTENNA subtable
  anttab = table(ms.getkeyword("ANTENNA"));
  stations = list(anttab.getcol('NAME'));
  antpos = anttab.getcol('POSITION');
  observatory = anttab.getcol("STATION")[0];
  anttab.close();

  # trim away longest antenna name prefix
  while stations[0]:
    if len([st for st in stations[1:] if st[0] == stations[0][0]]) != len(stations)-1:
      break;
    stations = [ st[1:] for st in stations ];

  # make IfrSet object
  return IfrSet(stations,positions=antpos,observatory=observatory);

