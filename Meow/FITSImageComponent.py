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
from SixpackComponent import *

class FITSImageComponent (SixpackComponent):
  """FitsImageComponent is a sixpack component that uses the FITSImage node 
  to get an image from a FITS file.
  """;
  
  def __init__ (self,ns,name,filename,direction=None,cutoff=1.0,fluxscale=None):
    """Constructor args:
        ns:         the node scope
        name:       the source name
        filename:   FITS file name
        direction:  direction to source (a Direction object). If None, the
                    direction from the FITS file will be used.
        cutoff:     cutoff quotient. FITSImage will find the smallest box
          containing all points with cutoff*100% of the peak flux, and
          return only that box.
        fluxscale:  if specified, image flux will be rescaled by that factor
    """;
    SixpackComponent.__init__(self,ns,name,direction=direction,fluxscale=fluxscale);
    self._filename = str(filename);
    self._cutoff   = float(cutoff);
    
  def sixpack (self):
    image = self.ns.image;
    if not image.initialized():
      image <<= Meq.FITSImage(filename=self._filename,cutoff=self._cutoff);
    return image;
    
    
