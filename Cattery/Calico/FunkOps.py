# -*- coding: utf-8 -*-


"""This is a library of "funklet operations".
Each function here is compatible with ParmTab.apply().
Input argument is a ParmTab.FunkSlie object.
Return value is a list of output funklets.
""";
from __future__ import absolute_import
from __future__ import print_function
from __future__ import division

import copy

from Timba.Meq import meq
from Timba import dmi
from Timba import mequtils


def average (funkslice):
  """Reduction function to replace all funklets in a slice with their mean.
  This is the canonical example of a reduction function.
  """
  outfunk = funkslice[0];
  shape = getattr(outfunk.coeff,'shape',None);   # shape will be None when coeff is just a float
  # accumulate sum of coefficients and take mean
  for funk in funkslice[1:]:
    outfunk.coeff += funk.coeff;
    # adjust domain to envelope
    for axis in funkslice.slice_axes:
      outfunk.domain[axis] = (min(outfunk.domain[axis][0],funk.domain[axis][0]),
                              max(outfunk.domain[axis][1],funk.domain[axis][1]));
  outfunk.coeff /= len(funkslice); #mean, float
  return [ outfunk ];

def linear_interpol (funkslice):
  """Reduction function to replace a set of funklets with piecewise linear interpolations. If the original
  N funklet c00 coefficients C0,C1,... are defined for domains [a0,b0] [a1,b1] ..., then the output N-1 funklets
  will be defined on domains [a0,(a1+b1)/2],[(a1+b1)/2,(a2+b2)/2],..., and will correspond to linear slopes
  passing through C0 at (a0+b0)/2, C1 at (a1+b1)/2, etc.
  """
  # transform only available with more than two funklets
  if len(funkslice) < 2:
    return funkslice;
  if funkslice.rank > 1:
    raise TypeError("linear interpolation only available for rank-1 slices");
  iaxis0 = funkslice.slice_iaxes[0];
  axis0 = funkslice.slice_axes[0];
  output = [];
  for ifunk0,funk0 in enumerate(funkslice[:-1]):
    # funklet for next domain
    funk1 = funkslice[ifunk0+1];
    dom0,dom1 = funk0.domain,funk1.domain;
    # centerpoint and c00 of ourselves
    x0 = sum(dom0[axis0])/2;
    c0 = funk0.coeff if isinstance(funk0.coeff,float) else funk0.coeff.ravel()[0];
    # centerpoint and c00 of next domain
    x1 = sum(dom1[axis0])/2;
    c1 = funk1.coeff if isinstance(funk1.coeff,float) else funk1.coeff.ravel()[0];
    # make output domain
    out_domain = copy.copy(dom0);
    # output domain boundaries are original subdomains' centers, with the exception of the first
    # and the last funklet, in which case the domain needs to extend to the edge of the first/last
    # subdomain
    a0 = dom0[axis0][0] if funk0 is funkslice[0] else x0;
    a1 = dom1[axis0][1] if funk1 is funkslice[-1] else x1;
    out_domain[axis0] = (a0,a1);
    # now make output polc
    output.append(meq.polc(coeff=[c0,c1-c0],domain=out_domain,offset=x0,scale=x1-x0,axis_index=iaxis0))
  return output;

def force_rank0 (funkslice):
  """Reduction function to reduce the polynomial rank of a set of funklets"""
  output = [];
  for funk in funkslice:
    c00 = funk.coeff if isinstance(funk.coeff,float) else funk.coeff.ravel()[0];
    output.append(meq.polc(coeff=c00,domain=funk.domain));
  return output;

_sub = dict([(a+b+c,b+c+':'+a) for a in 'ri' for b in 'xy' for c in 'xy' ]);

def rename (funkslice):
  fields = funkslice.name.rsplit(':',1);
  if len(fields) < 2:
    return None;
  global _sub;
  rep = _sub.get(fields[-1]);
  if rep:
    fields[-1] = rep;
    return [':'.join(fields)] + list(funkslice);
  return None;

def make_infinite_domain (funkslice):
  outlist = [];
  for funk in funkslice:
    for axis in funkslice.slice_axes:
      funk.domain[axis] = (-1e+99,1e+99);
  return list(funkslice);
