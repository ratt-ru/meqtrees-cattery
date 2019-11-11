# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division
from Timba import pynode,mequtils
from Timba.Meq import meq
import numpy
import traceback

from Siamese.OMS.emss_beams.InterpolatedVoltageBeam import unite_multiple_shapes



class ClippedSincBeam(pynode.PyNode):
  """Implements a sin(r)/r beam (where r=scale*sqrt(l^2+m^2)).
  Clips the beam to 0 where l<0 or m <0.
  """
  
  def update_state (self,mystate):
    """Standard function to update our state""";
    mystate('scale',265.667/1.4e+9);
    self._freqaxis = mequtils.get_axis_number("freq");
  
  def get_result (self,request,*children):
    try:
      if len(children) < 1:
        raise TypeError("at least one child is expected");
      # now, figure out the lms 
      # lm may be a 2/3-vector or an Nx2/3 tensor
      lm = children[0];
      dims = getattr(lm,'dims',[len(lm.vellsets)]);
      if len(dims) == 2 and dims[1] in (2,3):
        nsrc,nlm = dims;
        tensor = True;
      elif len(dims) == 1 and dims[0] in (2,3):
        nsrc,nlm = 1,dims[0];
        tensor = False;
      else:
        raise TypeError("expecting a 2/3-vector or an Nx2/3 matrix for child 0 (lm)");
      # pointing offsets (child 1) are optional
      if len(children) > 1:
        dlm = children[1];
        if len(dlm.vellsets) != 2:
          raise TypeError("expecting a 2-vector for child 1 (dlm)");
        dl,dm = dlm.vellsets[0].value,dlm.vellsets[1].value;
      else:
        dl = dm = None;
      # turn freq into an array of the proper shape
      freq = request.cells.grid.freq;
      freqshape = [1]*(self._freqaxis+1);
      freqshape[self._freqaxis] = len(freq);
      freq = freq.reshape(freqshape);
      # now compute the beam per source
      vellsets = [];
      for isrc in range(nsrc):
        # get l/m for this sources
        l,m = lm.vellsets[isrc*nlm].value,lm.vellsets[isrc*nlm+1].value;
        # apply pointing offset, if given
        if dl is not None:
          dl,l,dm,m = unite_multiple_shapes(dl,l,dm,m);
          l = l - dl
          m = m - dm
        l,m,fq = unite_multiple_shapes(l,m,freq);
        # compute distance (scaled by frequency)
        dist = fq*self.scale*numpy.sqrt(l*l+m*m);
        # compute sinc function
        E = meq.vells(dist.shape);
        E[...] = numpy.sin(dist)/dist;
        E[dist==0] = 1;
        E[(l<0)|(m<0)] = 0;
        vellsets.append(meq.vellset(E));
      # form up result and return
      result = meq.result(vellsets[0],cells=request.cells);
      result.vellsets[1:] = vellsets[1:];
      # result.dims = (nsrc,) if tensor else (1.);
      return result;
      #res0 = children[0];
      #res1 = children[1];
      #cells = res0.cells
      #l = res0.vellsets[0].value;
      #m = res1.vellsets[0].value;
      #shape = l.shape
      #for i in range(shape[0]):
        #if l[i,0] < 0.0 or m[i,0] < 0.0:
          #l[i,0] = 0.0
        #else:
          #dist = 265.667 * math.sqrt(l[i,0]*l[i,0] + m[i,0]*m[i,0])
          ##  265.667 is scaling factor to get 0.707 voltage response at 18 arcmin from pointing position (nominal FWHM=36arcmin at 1.4 GHz)
          #l[i,0] = math.sin(dist) / dist
        
      #vellsets = [];
      #vellsets.append(meq.vellset(l))
      #res = meq.result(cells=cells)
      #res.vellsets = vellsets
      #return res
    except:
      traceback.print_exc();
      raise;

