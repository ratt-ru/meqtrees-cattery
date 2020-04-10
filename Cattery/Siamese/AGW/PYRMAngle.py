# A MeqTrees PyNode to de-rotate ionosphee effects from synthesis telescope
# data in Measurement Set format

#  (c) 2019.                             (c) 2019.
#  National Research Council             Conseil national de recherches
#  Ottawa, Canada, K1A 0R6               Ottawa, Canada, K1A 0R6
#
#  This software is free software;       Ce logiciel est libre, vous
#  you can redistribute it and/or        pouvez le redistribuer et/ou le
#  modify it under the terms of          modifier selon les termes de la
#  the GNU General Public License        Licence Publique Generale GNU
#  as published by the Free              publiee par la Free Software
#  Software Foundation; either           Foundation (version 3 ou bien
#  version 2 of the License, or          toute autre version ulterieure
#  (at your option) any later            choisie par vous).
#  version.
#
#  This software is distributed in       Ce logiciel est distribue car
#  the hope that it will be              potentiellement utile, mais
#  useful, but WITHOUT ANY               SANS AUCUNE GARANTIE, ni
#  WARRANTY; without even the            explicite ni implicite, y
#  implied warranty of                   compris les garanties de
#  MERCHANTABILITY or FITNESS FOR        commercialisation ou
#  A PARTICULAR PURPOSE.  See the        d'adaptation dans un but
#  GNU General Public License for        specifique. Reportez-vous a la
#  more details.                         Licence Publique Generale GNU
#                                        pour plus de details.
#
#  You should have received a copy       Vous devez avoir recu une copie
#  of the GNU General Public             de la Licence Publique Generale
#  License along with this               GNU en meme temps que ce
#  software; if not, contact the         logiciel ; si ce n'est pas le
#  Free Software Foundation, Inc.        cas, communiquez avec la Free
#  at http://www.fsf.org.                Software Foundation, Inc. au
#                                                http://www.fsf.org.
#
#  email:                                courriel:
#  business@hia-iha.nrc-cnrc.gc.ca       business@hia-iha.nrc-cnrc.gc.ca
#
#  National Research Council             Conseil national de recherches
#      Canada                               Canada
#  Herzberg Institute of Astrophysics    Institut Herzberg d'astrophysique
#  5071 West Saanich Rd.                 5071 West Saanich Rd.
#  Victoria BC V9E 2E7                   Victoria BC V9E 2E7
#  CANADA                                        CANADA
#

# standard imports
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

from Timba.Meq import meq
from Timba import pynode
import numpy
from scipy.interpolate import interp1d
from string import split, strip

LIGHT_SPEED = 299792458.0     # metres per sec

class PyGetRMAngle (pynode.PyNode):
  """attempts to getFaraday rotation Measure from RMextract/ALBUS output file """
  def __init__ (self,*args):
    pynode.PyNode.__init__(self,*args);
    self.set_symdeps('domain')

    self.albus_filename = ''

    global _albus_time
    global _rm

  def update_state (self,mystate):
    """Standard function to update our state""";
    mystate('albus_filename');
    # Check filename arguments
    if isinstance(self.albus_filename,str):
      self.getdata(self.albus_filename)
      print('Siamese/AGW/PyGetRMAngle albus_filename = ',self.albus_filename)

  def getdata(self,filename):
    global _albus_time
    global _albus_rm

    text = open(filename, 'r').readlines()
    i = 0
    start_time = 0.0
    # skip over all stuff before actual data
    while(text[i].find('time_width') < 0):
       if text[i].find('start and end times') >= 0:
          info = split(strip(text[i]))
          start_time = float(info[4]) 
       i = i+1
    start = i+1
    time=[]
    rm = []
    L = len(text)
    # get actual data
    for i in range(start,L):
      try:
        info = split(strip(text[i]))
#       if int(info[2]) == 0:
        time.append(float(info[3]))
        data = info[8]
        try:
          rm.append(float(data))
        except:
          last = len(data) - 1
          data = data[1:last]
#         print 'adjusted info', data
          rm.append(float(data))

#  for processing of test source with fixed incoming iono RM 
#       rm.append(50.0)
      except:
        pass
    _albus_time = numpy.array(time) + start_time
    _albus_rm  = numpy.array(rm)
#   print '_albus_time ', _albus_time
#   print '_albus_rm ', _albus_rm

  def get_result (self,request):
    global _albus_time
    global _albus_rm

    # get time and frequency for the tile
    cells = request.cells
    cells_shape = meq.shape(cells);
    time = cells.grid.time
    freq = cells.grid.freq

    # use scipy to interpolate and get rm for MS times
    f = interp1d(_albus_time, _albus_rm)
    rotation_measure = f(time) 
#   print 'rotation_measure', rotation_measure

    factor = LIGHT_SPEED / freq
    out_array = meq.vells(cells_shape)
#   print 'out_array.shape', out_array.shape

    metre_sq = factor * factor
    for j in range(cells_shape[0]):
      if cells_shape[0] > 1:
        angle =  rotation_measure[j] * metre_sq
      else:
        angle =  rotation_measure * metre_sq
#     print 'j angle ', j, angle
      out_array[j,:] = angle

    vellsets = [];
    vellsets.append(meq.vellset(out_array))
    res = meq.result(cells = request.cells)
    res.vellsets = vellsets
    return res
