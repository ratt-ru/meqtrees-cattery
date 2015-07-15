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
from Timba.Meq import meq
from SkyComponent import *
import Jones
import Context

STOKES = ("I","Q","U","V");

class PointSource(SkyComponent):
  def __init__(self,ns,name,direction,
               I=0.0,Q=None,U=None,V=None,
               spi=None,freq0=None,
               RM=None):
    """Creates a PointSource with the given name; associates with
    direction.
    'direction' is a Direction object or a (ra,dec) tuple
    'I' is I flux (constant, node, or Meow.Parm)
    Optional arguments:
      'Q','U','V' are constants, nodes, or Meow.Parms. If none of the three
          are supplied, an unpolarized source representation is used.
      'spi' and 'freq0' are constants, nodes or Meow.Parms. If both are
          supplied, a spectral index is added, otherwise the fluxes are
          constant in frequency. 'spi' may be a list of such, in which case
          higher-order indices are used.
      'RM' is rotation measure (constants, nodes or Meow.Parms). If None,
          no intrinsic rotation measure is used.
    """;
    SkyComponent.__init__(self,ns,name,direction);
    self._add_parm('I',I,tags="flux");
    # check if polarized
    # NB: disable for now, as Sink can't handle scalar results
    # self._polarized = True;
    self._polarized = Q is not None or U is not None or V is not None or RM is not None;
    self._add_parm('Q',Q or 0.,tags="flux pol");
    self._add_parm('U',U or 0.,tags="flux pol");
    self._add_parm('V',V or 0.,tags="flux pol");
    self._has_rm = RM != 0 and RM is not None;
    if self._has_rm:
      self._add_parm("RM",RM,tags="pol");
    # see if a spectral index is present (freq0!=0 then), create polc
    self._has_spi = spi != 0 and spi is not None;
    if self._has_spi:
      if isinstance(spi,(list,tuple)):
        self._nspi = len(spi);
        for i,x in enumerate(spi):
          self._add_parm('spi' if not i else 'spi_%d'%(i+1),x,tags="spectrum");
      else:  
        self._nspi = 1;
        self._add_parm('spi',spi or 0.,tags="spectrum");
    if self._has_spi or self._has_rm:
      # for freq0, use placeholder node for first MS frequency,
      # unless something else is specified
      self._add_parm('spi_fq0',freq0 or (ns.freq0 ** 0),tags="spectrum");
      self._freq0 = freq0

    # see if intrinsic rotation measure is present
    # flag: flux is fully defined by constants (i.e. no rotation measure).
    # we can still have a spectral index on top of this
    self._constant_flux = not self._has_rm and not \
                          [ flux for flux in STOKES if not self._is_constant(flux) ];

  def is_polarized(self):
    return self._polarized

  def stokes (self,st):
    """Returns flux node for this source. 'st' must be one of 'I','Q','U','V'.
    (This is the flux after RM has been applied, but without spi).
    If flux was defined as a constant, returns constant value, not node!
    """;
    if st not in STOKES:
      raise ValueError,"st: must be one of 'I','Q','U','V'";
    if self._constant_flux:
      return self._get_constant(st);
    # rotation measure rotates Q-U with frequency. So if no RM is given,
    # or if we're just asked for an I or V node, return it as is
    if st == "I" or st == "V" or not self._has_rm:
      return self._parm(st);
    else:
      stokes = self.ns[st+'r'];
      if stokes.initialized():
        return stokes;
      q = self._parm("Q");
      u = self._parm("U");
      # first rotate back to wavelength 0 lambda
      # we assume that reference frequency is equal to that of spectral index
      if self._freq0 is None:
        farot_ref = self.ns.farot_ref << 0.0;
      else:
        lambda_ref2 = self.ns.lambda_ref2 << Meq.Sqr(2.99792458e+8/self._parm('spi_fq0'))
        # remember that polarization position angle = 0.5 * arctan(U/Q) so 
        # factor 2 needed when determining rotation angles for Q and U 
        # coordinate system
        farot_ref = self.ns.farot_ref << self._parm("RM")*lambda_ref2 * -2;
      cosf_ref = self.ns.cos_farot_ref << Meq.Cos(farot_ref);
      sinf_ref = self.ns.sin_farot_ref << Meq.Sin(farot_ref);
      self.ns.Q_ref << cosf_ref*q - sinf_ref*u;
      self.ns.U_ref << sinf_ref*q + cosf_ref*u;
      # rotate forward to requested wavelength
      freq = self.ns0.freq ** Meq.Freq;
      iwl2 = self.ns0.wavelength2 << Meq.Sqr(2.99792458e+8/freq);
      # remember that polarization position angle = 0.5 * arctan(U/Q) so 
      # factor 2 needed when determining rotation angles for Q and U 
      # coordinate system
      farot = self.ns.farot << self._parm("RM")*iwl2 * 2.0;
      cosf = self.ns.cos_farot << Meq.Cos(farot);
      sinf = self.ns.sin_farot << Meq.Sin(farot);
      self.ns.Qr << cosf*self.ns.Q_ref - sinf*self.ns.U_ref;
      self.ns.Ur << sinf*self.ns.Q_ref + cosf*self.ns.U_ref;

      return self.ns[st+'r'];

  def norm_spectrum (self):
    """Returns spectrum normalized to 1 at the reference frequency.
    Flux should be multiplied by this to get the real spectrum""";
    if not self._has_spi:
      return 1;
    nsp = self.ns.norm_spectrum;
    if not nsp.initialized():
      freq = self.ns0.freq ** Meq.Freq;
      fr = self.ns.freqratio << freq/self._parm('spi_fq0');  
      if self._nspi == 1:
        nsp << Meq.Pow(fr,self._parm('spi'));
      else:
        # multi-term spectral index
        logfr = self.ns.logfreqratio << Meq.Log(fr);
        nsp << Meq.Pow(fr,Meq.Add(self._parm('spi'),self._parm('spi_2')*logfr,
          *[ self._parm('spi_%d'%(i+1))*Meq.Pow(logfr,i) for i in range(2,self._nspi)] ));
    return nsp;

  def coherency_elements (self,observation):
    """helper method: returns the four components of the coherency matrix""";
    i,q,u,v = [ self.stokes(st) for st in STOKES ];
    # diagonal = (len(Context.active_correlations) == 2);
    diagonal = False;  # the above was bothersome -- even if we only use 2 corrs, we still want to do all intermediate computations in 2x2
    if observation.circular():
      if self._constant_flux:
        return (i+v,0,0,i-v) if diagonal else (i+v,complex(q,u),complex(q,-u),i-v);
      rr = self.ns.rr ** (self.stokes("I") + self.stokes("V"));
      if diagonal:
        rl = lr = 0;
      else:
        rl = self.ns.rl ** Meq.ToComplex(self.stokes("Q"),self.stokes("U"));
        lr = self.ns.lr ** Meq.Conj(rl);
      ll = self.ns.ll ** (self.stokes("I") - self.stokes("V"));
      return rr,rl,lr,ll;
    else:
      if self._constant_flux:
        return (i+q,0,0,i-q) if diagonal else (i+q,complex(u,v),complex(u,-v),i-q);
      xx = self.ns.xx ** (self.stokes("I") + self.stokes("Q"));
      if diagonal:
        xy = yx = 0;
      else:
        xy = self.ns.xy ** Meq.ToComplex(self.stokes("U"),self.stokes("V"));
        yx = self.ns.yx ** Meq.Conj(xy);
      yy = self.ns.yy ** (self.stokes("I") - self.stokes("Q"));
      return xx,xy,yx,yy;

  def brightness_static (self,observation=None):
    """If the brightness matrix is constant (i.e. no time-freq dependence), returns it.
    Else returns None.""";
    if not self._constant_flux or self._has_spi:
      return None;
    observation = observation or Context.observation;
    if not observation:
      raise ValueError,"observation not specified in global Meow.Context, or in this function call";
    # get a list of the coherency elements
    xx,xy,yx,yy = self.coherency_elements(observation);
    return [[xx,xy],[yx,yy]];

  def brightness (self,observation=None,nodes=None,always_matrix=False):
    """Returns the brightness matrix for a point source.
    'observation' argument is used to select a linear or circular basis;
    if not supplied, the global context is used.
    If always_matrix=True, returns matrix even if source is unpolarized.
    """;
    observation = observation or Context.observation;
    if not observation:
      raise ValueError,"observation not specified in global Meow.Context, or in this function call";
    coh = nodes or self.ns.brightness;
    if not coh.initialized():
      # if not polarized, just return I
      if not always_matrix and not self._polarized:
        if self._has_spi:
          if self._constant_flux:
            coh << Context.unitCoherency(self.stokes("I"))*self.norm_spectrum();
          else:
            coh << Meq.Multiply(self.stokes("I"),self.norm_spectrum(),Context.unitCoherency);
        else:
          coh = Context.unitCoherency(self.stokes("I"));
      else:
        coh_els = self.coherency_elements(observation);
        if self._constant_flux:
          if self._has_spi:
            coh1 = self.ns.coh1 << Meq.Matrix22(*[Context.unitCoherency(x) for x in coh_els]);
            coh << coh1*self.norm_spectrum();
          else:
            coh << Meq.Matrix22(*[Context.unitCoherency(x) for x in coh_els]);
        else:
          coh1 = self.ns.coh1 << Meq.Matrix22(*coh_els);
          if self._has_spi:
            coh << Meq.Multiply(coh1,self.norm_spectrum(),Context.unit_coherency);
          else:
            coh = Context.unitCoherency(coh1);
    return coh;

  def coherency (self,array=None,observation=None,nodes=None,**kw):
    brightness = self.brightness(observation,nodes=nodes,**kw);
    return lambda p,q:brightness;

  def is_station_decomposable (self):
    """Check the type -- subclasses are not necessarily decomposable.""";
    return type(self) == PointSource;

  def sqrt_coherency (self,observation):
    # Cholesky decomposition of the coherency matrix into
    # UU*, where U is lower-triangular
    coh = self.ns.sqrt_coh;
    if not coh.initialized():
      # if unpolarized, then matrix is scalar -- simply take the square root
      if not self._polarized:
        flux = self.stokes("I")*Context.unit_coherency*self.norm_spectrum();
        if isinstance(flux,(int,float,complex)):
          coh << math.sqrt(flux);
        else:
          coh << Meq.Sqrt(flux);
      else:
        # circular and linear matrices have the same form, only QUV is swapped around.
        # So the code below forms up the matrix assuming linear polarization, while
        # here for the circular case we simply swap the quv variables around
        if observation.circular():
          i,q,u,v = self.stokes("I"),self.stokes("V"),self.stokes("Q"),self.stokes("U");
        else:
          i,q,u,v = self.stokes("I"),self.stokes("Q"),self.stokes("U"),self.stokes("V");
        # form up matrix
        if self._constant_flux:
          norm = math.sqrt(.5/(i+q));
          c11 = (i+q)*norm;
          c12 = 0;
          c21 = complex(u,-v)*norm;
          c22 = math.sqrt(i**2-q**2-u**2-v**2)*norm;
          if self._has_spi:
            coh1 = self.ns.sqrt_coh1 << Meq.Matrix22(c11,c12,c21,c22);
            coh << coh1*Meq.Sqrt(self.norm_spectrum());
          else:
            coh << Meq.Matrix22(c11,c12,c21,c22);
        else:
          c11 = self.ns.sqrt_coh(1,1) << i+q;
          c12 = 0;
          c21 = self.ns.sqrt_coh(2,1) << Meq.ToComplex(u,-v);
          i2 = self.ns.I2 << Meq.Sqr(i);
          q2 = self.ns.Q2 << Meq.Sqr(q);
          u2 = self.ns.U2 << Meq.Sqr(u);
          v2 = self.ns.V2 << Meq.Sqr(v);
          c22 = self.ns.sqrt_coh(2,2) << Meq.Sqrt(
            self.ns.I2_QUV2 << self.ns.I2 - (self.ns.QUV2 << Meq.Add(q2,u2,v2))
          );
          # create unnormalized matrix + normalization term
          self.ns.sqrt_coh1 << Meq.Matrix22(c11,c12,c21,c22);
          norm = self.ns.sqrt_coh_norm << Meq.Sqrt(.5*self.norm_spectrum()/c11);
          # and finally form up product
          coh << norm * self.ns.sqrt_coh1;
    return coh;

  def sqrt_visibilities (self,array=None,observation=None,nodes=None):
    self.using_station_decomposition = True;
    observation = Context.get_observation(observation);
    array = Context.get_array(array);
    if nodes is None:
      nodes = self.ns.sqrt_visibility.qadd(observation.radec0());
    stations = array.stations();
    if nodes(stations[0]).initialized():
      return nodes;
    # create decomposition nodes
    sqrtcoh = self.sqrt_coherency(observation);
    # get K jones per station
    if self.direction is observation.phase_center:
      for p in array.stations():
        nodes(p) << Meq.Identity(sqrtcoh);
    # else apply KJones
    else:
      Kj = self.direction.KJones(array=array,dir0=observation.phase_center);
      for p in array.stations():
        nodes(p) << Meq.MatrixMultiply(Kj(p),sqrtcoh);
    return nodes;
