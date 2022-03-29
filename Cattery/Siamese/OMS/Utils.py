# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

import json
import os
import re

import Kittens.utils
_verbosity = Kittens.utils.verbosity(name="vb");
dprint = _verbosity.dprint;
dprintf = _verbosity.dprintf;

class json_beamconfig_reader():
  def __init__(self, beamsets, station_names, verbosity=0):
    self.__verbose = verbosity
    _verbosity.set_verbose(self.__verbose)
    self.__beamsets = beamsets
    self.__station_names = station_names
    self.__station_types = {
      "patterns": {},
      "define-stationtypes": {}
    }
    self.__station_dependent_beams = False
    self.__load_patterns()

  @property
  def is_station_dependent(self):
    return self.__station_dependent_beams

  @property
  def helpstring(self):
    example_use = \
      "Patterns should resemble 'prefix$(stype)infix$(corr)$infix$(reim).fits' to specify " \
      "a set of fits files containing real and imaginary parts for each correlation hand, e.g. " \
      "mybeam_$(stype)_$(corr)_$(reim).fits specifying mybeam_meerkat_xx_re.fits, mybeam_meerkat_xx_im.fits, " \
      "mybeam_meerkat_xy_re.fits... Lists of such patterns with different frequency coverage may be specified.\n"\
      "Station types can be specified using Json configuration files, containing: \n"\
      "{'lband': {\n" \
      "   'patterns': {\n" \
      "       'cmd::default': ['$(stype)_$(corr)_$(reim).fits',...],\n" \
      "   },\n" \
      "  'define-stationtypes': {\n" \
      "      'cmd::default': 'meerkat',\n" \
      "      'ska000': 'ska'\n" \
      "  },\n" \
      "  ...\n" \
      "}\n" \
      "This will substitute 'meerkat' for all antennas but ska000, with 'meerkat_$(corr)_$(reim).fits' " \
      "whereas beams for ska000 will be loaded from 'ska_$(corr)_$(reim).fits' in this example.\n" \
      "The station name may be specified as regex by adding a '~' infront of the pattern to match, e.g " \
      "'~ska[0-9]{3}': 'ska' will assgign all the 'ska' type to all matching names such as ska000, ska001, ..., skaNNN.\n" \
      "Each station type in the pattern section may specify a list of patterns for different frequency ranges.\n" \
      "Multiple keyed dictionaries such as this may be specified within one file. They will be treated as chained " \
      "configurations, adding more patterns and station-types to the first such block.\n" \
      "Warning: Once a station is type-specialized the type applies to **ALL** chained blocks!\n" \
      "Blocks from more than one config file can be loaded by comma separation, e.g. " \
      "'conf1.json,conf2.json,...', however no block may define multiple types for any station.\n" \
      "If patterns for a particular station type already exists more patterns are just appended to the existing list.\n" \
      "Warning: where multiple patterns specify the same frequency range the first such pattern closest to the MS " \
      "SPW frequency coverage will be loaded.\n" \
      "If no configuration file is provided the pattern may not contain $(stype) -- station independence is assumed. " \
      "This is the same as specifing the following config: \n" \
      "{'lband': {\n" \
      "   'patterns': {\n" \
      "       'cmd::default': ['$(corr)_$(reim).fits',...],\n" \
      "   },\n" \
      "  'define-stationtypes': {\n" \
      "      'cmd::default': 'cmd::default',\n" \
      "  }\n" \
      "}\n" \
      "'corr' above may be one of 'corr', 'CORR', 'xy', 'XY' to cast the hands to lower or UPPER case.\n" \
      "'reim' above may be one of 'reim', 'REIM', 'ReIm', 'realimag', 'REALIMAG', 'RealImag' to cast the " \
      "real or imaginary substitutions to lower case, CamelCase or UPPER case respectively.\n" \
      "'stype' above may be one of 'stype' or 'STYPE' to cast the station/antenna type id to the given name " \
      "as defined in the Json file. This is optional depending on whether station types are being used.\n"
    return example_use

  def get_stationtype(self, station):
        return self.station_types["define-stationtypes"].get(
                    station, 
                    self.station_types["define-stationtypes"]["cmd::default"]
        )

  def get_beamsets(self, station):
      return self.station_types["patterns"].get(
          self.__get_stationtype(station),
          self.station_types["patterns"]["cmd::default"]
      )

  def __load_patterns(self):
    for bs in self.__beamsets:
      if os.path.splitext(bs)[1] == ".json":
          if os.path.exists(bs):
              with open(bs) as fbs:
                  vels = json.loads(fbs.read())
              if not isinstance(vels, dict):
                  raise ValueError("Station type config should contain dictionary blocks. " + 
                                   json_beamconfig_reader.example_use())
              for key in vels:
                  if not set(vels[key].keys()).issubset(["patterns", "define-stationtypes"]):
                      raise ValueError("Station type config blocks may only contain keys "
                                      "'patterns' and 'define-stationtypes'. " + self.helpstring)
                  for st_t, patterns in vels[key].get("patterns", {}).items():
                      if not isinstance(patterns, list):
                          patterns = [patterns]
                      self.__station_types["patterns"][st_t] = \
                          self.__station_types["patterns"].get(st_t, []) + patterns

                  for st in filter(lambda key: key.find("~") < 0,
                                  vels[key].get("define-stationtypes", {}).keys()):
                      st_t = vels[key]["define-stationtypes"][st]
                      if st in self.__station_types["define-stationtypes"] and \
                          self.__station_types["define-stationtypes"][st] != st_t:
                          raise ValueError(f"Ambiguous redefinition of station type for '{st}' while "
                                          f"parsing '{bs}'. Check your configuration chain. " + self.helpstring)
                      self.__station_types["define-stationtypes"][st] = st_t

                  for wildcard in map(lambda x: x[1:], 
                                      filter(lambda key: key.find("~") == 0,
                                          vels[key].get("define-stationtypes", {}).keys())):
                      fcmatches = list(filter(lambda st: re.match(wildcard, st), self.ms.StationNames))
                      if len(fcmatches) == 0:
                          raise ValueError(f"No station name matches regex '{wildcard}' while parsing '{bs}'. "
                                            f"Check your station types config!")
                      for st in fcmatches:
                          st_t = vels[key]["define-stationtypes"]["~" + wildcard]
                          if st in self.__station_types["define-stationtypes"] and \
                              self.__station_types["define-stationtypes"][st] != st_t:
                              raise ValueError(f"Ambigous redefinition of station type for '{st}' while "
                                                f"parsing wildcard station name '{wildcard}' in '{bs}'. "
                                                f"Check your configuration chain. " + self.helpstring)
                          self.__station_types["define-stationtypes"][st] = st_t
          else:
              raise FileNotFoundError(f"Station beam pattern map config file '{bs}' does not exist")
      else: # treat as a beam file pattern -- it looks like the code downstream can handle scalar beams
          self.__station_types["patterns"]["cmd::default"] = \
              self.__station_types["patterns"].get("cmd::default", []) + [bs]
      
    self.__station_types["define-stationtypes"].setdefault("cmd::default", "cmd::default")
    self.__station_types["patterns"].setdefault("cmd::default", [])

    for st in self.__station_names:
        beamsets = self.get_beamsets(st)
        this_st_t = self.get_stationtype(st)
        if beamsets == []:
            raise ValueError(f"EJones via FITS beams are enabled in your parset, "
                              f"but no beam patterns are specified for station {st}. "
                              f"Please check your config")

        if this_st_t == "cmd::default" and \
            any(map(lambda x: "$(" + x + ")" in bs,
                    ["stype","STYPE"])):
            raise ValueError(f"One or more patterns for station {st} requests substition on type "
                              f"but you have not specified a type for this station via types configuration file. "
                              f"Did you want to load a configuration file as well? " + _example_use)

    # now we may have many patterns in our little database, but what really determines
    # if we should use station dependent beams is if more than just the default group of
    # stations have been assigned a type to use other than the default
    if len(self.__station_types["define-stationtypes"]) > 1 and \
        any(map(lambda st: self.__station_types["define-stationtypes"][st] != 
                            self.__station_types["define-stationtypes"]["cmd::default"],
                self.__station_types["define-stationtypes"].keys())):
        self.__station_dependent_beams = True
    
    if self.__station_dependent_beams:
        dprint(0, "Using station-dependent E Jones for the array - this "
                  "may take longer to interpolate depending on how many "
                  "unique elements are in the array")
    else:
        dprint(0, "Using station-independent E Jones for the array")

def substitute_pattern (filename_pattern,**substitutions):
  """Substitutes $(key) and $var instances in the given filename pattern, with values from the
  substitutions dict.
  """
  import re
  filename = filename_pattern;
  # loop over substitutions, longest to shortest
  # delimit substitutions with "\n" so that they don't interfere with subsequent word boundaries
  from past.builtins import cmp
  from functools import cmp_to_key
  for key,value in sorted(list(substitutions.items()),key=cmp_to_key(lambda a,b:cmp(len(b[0]),len(a[0])))):
    filename = re.sub("\\$(%s\\b|\\(%s\\))"%(key,key),"\n"+value+"\n",filename);
  filename = re.sub("\\$\\$","$",filename);
  filename = filename.replace("\n","");
  return filename;
