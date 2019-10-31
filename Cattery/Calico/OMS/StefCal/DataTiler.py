from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

import numpy
import operator
from .MatrixOps import *
from Timba.Meq import meq
from functools import reduce


class DataTiler (object):
  """Support class to handle subtiling of data, i.e. covering every axis of length N with K subtiles of length M=N/K.

  The following attributes are defined:
    
      datashape   = N1,N2,...                # original data shape
      subtiling   = M1,M2,...                # tilesizes
      subshape    = K1=N1/M1,K2=N2/M2,...    # tiling shape (number of tiles per axis)

      tiled_shape =  K1,M1,K2,M2,...         # intermediate shape into which a datashape may be reshaped

      # slice to convert from a subshape to a tiled_shape compatible expression
      tiling_slice  = :,numpy.newaxis,:,numpy.newaxis,... 
      
    The following methods are defined:
    
      tile_data(x):        reshapes x of shape datashape into subtiled_shape
      untile_data(x):      reshapes x of shape subtiled_shape into datashape
      tile_gain(x):        reshapes gain into subtiled shape (by inserting a numpy.newaxis for every second axis)
      reduce_subtiles(x):  reduces subtiled_shape to gain_shape by summing every second axis

    If the common case of a 1,1,... subtiling, all these can be an identity.
  """

  def __init__ (self,datashape,subtiling,original_datashape=None,force_subtiling=False):
    """Prepares a DataTiler for tiling 'datashape' with subtiles of size 'subtiling'.
    
    If force_subtiling is True, also subtiles axes where subtiling=1. 
    """
    self.datashape   = datashape;  
    self.subtiling   = subtiling;
    self.subshape    = tuple([ nd//nt for nd,nt in zip(datashape,subtiling) ]);
    # work out various stats
    # total number of slots in subshape
    self.total_slots  = reduce(operator.mul,self.subshape);
    unpadded_subshape = tuple([ int(math.ceil(nd/float(nt))) for nd,nt in zip(original_datashape,subtiling) ]);
    # number of "real" slots given the padding
    self.real_slots = reduce(operator.mul,unpadded_subshape);
    # number of padded slots
    self.padded_slots = self.total_slots - self.real_slots;
    # if subtiling is 1,1,... then override methods with identity relations
    if max(subtiling) == 1 and not force_subtiling:
      self.tiled_shape = self.datashape;
      self.tile_data = lambda x, dtype=None: x.astype(dtype) if ( not numpy.isscalar(x) and dtype and dtype != x.dtype ) else x
      self.untile_data = self.tile_tiling = self.reduce_tiles = identity_function;
      self.expand_tiling = self._expand_trivial_subshape;
      self.subtiled_axes = ();
      self.tiling_slice = ();
    else:
      self.tiled_shape = [];
      self.tiling_slice = [];
      self.subtiled_axes = [];
      for i,(ng,nt) in enumerate(zip(self.subshape,subtiling)):
        if nt>1 or force_subtiling:
          self.tiled_shape += [ng,nt];
          self.tiling_slice += [slice(None),numpy.newaxis];
          self.subtiled_axes.append(len(self.tiling_slice)-1);
        else:
          self.tiled_shape.append(ng);
          self.tiling_slice.append(slice(None));
      self.subtiled_axes = self.subtiled_axes[-1::-1];

  # define methods
  def tile_data (self,x,dtype=None):
    """Converts something of shape datashape into a subtiled_shape"""
    try:
      if numpy.isscalar(x) or x.size == 1:
          return x;
      x = x.reshape(self.tiled_shape);
      return x.astype(dtype) if ( dtype and dtype != x.dtype ) else x
    except:
      print('tile_data exception, tiled_shape:',self.tiled_shape,', arg:',getattr(x,'shape',()));
      raise;
    
  def untile_data (self,x):
    """Converts something of shape subtiled_shape back into datashape"""
    try:
      return x.reshape(self.datashape) if not (numpy.isscalar(x) or x.size == 1) else x;
    except:
      print('untile_data exception, datashape:',self.datashape,', arg:',getattr(x,'shape',()));
      raise;
    
  def tile_subshape (self,x):
    """Converts something of shape subshape into a subtiled_shape"""
    return x[self.tiling_slice] if not (numpy.isscalar(x) or x.size == 1) else x;
    
  def reduce_tiles (self,x,method='sum'):
    """reduces something of shape tiled_shape into a subshape by collapsing the M-axes""";
    try:
      if not (numpy.isscalar(x) or x.size == 1):
        for ax in self.subtiled_axes:
          x = getattr(x,method)(ax);
      return x;
    except:
      print('reduce_tiles exception, axes:',self.subtiled_axes,', arg:',getattr(x,'shape',()));
      raise;
    
  def expand_subshape (self,x,datashape=None,data_subset=None):
    """expands subshape to original data shape"""
    if numpy.isscalar(x):
      return x;
    a = numpy.empty(self.tiled_shape,dtype=x.dtype);
    a[...] = self.tile_subshape(x);
    if x.dtype == bool:
      b = numpy.zeros(datashape or self.datashape,bool);
    else:
      b = meq.complex_vells(datashape or self.datashape);
    b[...] = self.untile_data(a);
    return b[data_subset or ()];
    
  def _expand_trivial_subshape (self,x,datashape=None,data_subset=None):
    a = meq.complex_vells(x.shape);
    a[...] = x;
    return a[data_subset];
      