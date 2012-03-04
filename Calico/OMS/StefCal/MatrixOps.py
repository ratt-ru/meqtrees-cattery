# -*- coding: utf-8 -*-
import math
import numpy
from Timba.Meq import meq


identity_function = lambda x:x;

# returns 0 if an element is 0
def is_null (x):
  return (numpy.isscalar(x) and x == 0) or (x.size == 1 and x.flat[0] == 0);

# list of matrix indices in a flat 4-list representation
IJ2x2 = [ (i,j) for i in range(2) for j in range(2) ];

NULL_MATRIX = [0,0,0,0];

# extracts one matrix for key0 from a dict with keys of key0,i,j (such as the data, model and parm dicts),
# and converts it into a flat 4-list
def make_matrix (datadict,key0):
  return [ datadict.get((key0,i,j),0) for i,j in IJ2x2 ];


def matrix_multiply (A,B):
  """Multiplies two matrices given as flat 4-lists""";
  a11,a12,a21,a22 = A;
  b11,b12,b21,b22 = B;
  return a11*b11+a12*b21,a11*b12+a12*b22,a21*b11+a22*b21,a21*b12+a22*b22;

def matrix_conj (A):
  """Conjugates a matrix given as a flat 4-list""";
  return [ numpy.conj(A[i]) for i in 0,2,1,3 ];

def matrix_transpose (A):
  """Transposes a matrix given as a flat 4-list""";
  return A[0],A[2],A[1],A[3];

def matrix_add (A,B):
  """Sums two matrices given as a flat 4-list""";
  return [ a+b for a,b in zip(A,B) ];

def matrix_sub (A,B):
  """Subtracts two matrices given as a flat 4-list""";
  return [ a-b for a,b in zip(A,B) ];

def matrix_invert (A):
  """Inverts a matrix given as a flat 4-list""";
  a,b,c,d = A;
  idet = 1/(a*d-b*c);
  return d*idet,-b*idet,-c*idet,a*idet;

def array_to_vells (x,vellshape=None,expanded_slice=None):
  if expanded_slice is not None:
    a = meq.complex_vells(vellshape);
    a[...] = x[expanded_slice];
  else:
    a = meq.complex_vells(x.shape);
    a[...] = x;
  return a;

