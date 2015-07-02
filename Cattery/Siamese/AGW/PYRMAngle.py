from Timba.Meq import meq
from Timba import pynode
import numpy
from scipy.interpolate import interp1d
from string import split, strip

LIGHT_SPEED = 299792458.0     # metres per sec

class PyGetRMAngle (pynode.PyNode):
  """attempts to getFaraday rotation Measure from ALBUS output file """
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
      print 'Siamese/AGW/PyGetRMAngle albus_filename = ',self.albus_filename

  def getdata(self,filename):
    global _albus_time
    global _albus_rm

    text = open(filename, 'r').readlines()
    i = 0
    start_time = 0.0
    # skip over all stuff before actual data
    while(text[i].find('rel_time') < 0):
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
        if int(info[2]) == 0:
          time.append(float(info[3]))
          rm.append(float(info[8]))
      except:
        pass
    _albus_time = numpy.array(time) + start_time
    _albus_rm  = numpy.array(rm)

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

    factor = LIGHT_SPEED / freq
    out_array = meq.vells(cells_shape)

    metre_sq = factor * factor
    for j in range(cells_shape[0]):
      angle =  rotation_measure[j] * metre_sq
      out_array[j,:] = angle

    vellsets = [];
    vellsets.append(meq.vellset(out_array))
    res = meq.result(cells = request.cells)
    res.vellsets = vellsets
    return res
