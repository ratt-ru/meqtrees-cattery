# -*- coding: utf-8 -*-
import math
import numpy
from Timba.Meq import meq


identity_function = lambda x:x;

# returns 0 if an element is 0
def is_null (x):
  return x is None or ( not x if numpy.isscalar(x) else numpy.ma.is_masked(x) or (x.size == 1 and x.flat[0] == 0) );

# list of matrix indices in a flat 4-list representation
IJ2x2 = [ (i,j) for i in range(2) for j in range(2) ];

def NULL_MATRIX():
  return [0,0,0,0];

# extracts one matrix for key0 from a dict with keys of key0,i,j (such as the data, model and parm dicts),
# and converts it into a flat 4-list
def make_matrix (datadict,key0):
  return [ datadict.get((key0,i,j),0) for i,j in IJ2x2 ];

def matrix_copy (A):
  return [ x if numpy.isscalar(x) else x.copy() for x in A ];

def _mul (x,y):
  return 0 if is_null(x) or is_null(y) else x*y;
  
def matrix_multiply (A,B):
  """Multiplies two matrices given as flat 4-lists""";
  a11,a12,a21,a22 = A;
  b11,b12,b21,b22 = B;
  return [ _mul(a11,b11)+_mul(a12,b21),_mul(a11,b12)+_mul(a12,b22),_mul(a21,b11)+_mul(a22,b21),_mul(a21,b12)+_mul(a22,b22) ];

def matrix_conj (A):
  """Conjugates a matrix given as a flat 4-list""";
  return [ numpy.conj(A[i]) for i in 0,2,1,3 ];

def matrix_transpose (A):
  """Transposes a matrix given as a flat 4-list""";
  return A[0],A[2],A[1],A[3];

def matrix_add (A,B):
  """Sums two matrices given as a flat 4-list""";
  return [ a+b for a,b in zip(A,B) ];

def matrix_add1 (A,B):
  """Sums two matrices in place""";
  for i,b in enumerate(B):
    A[i] += b;
  return A;

def matrix_scale (A,c):
  """Scales a matrix by scalar c""";
  return [ _mul(a,c) for a in A ];

def matrix_scale1 (A,c):
  """Scales a matrix by scalar c""";
  for i,a in enumerate(A):
    if is_null(a) or is_null(c):
      A[i] = 0;
    else:
      A[i] *= c;
  return A;

def matrix_sub (A,B):
  """Subtracts two matrices given as a flat 4-list""";
  return [ a-b for a,b in zip(A,B) ];

def matrix_negate (A):
  """Negates a matrix""";
  return [ -a for a in A ];

def matrix_invert (A,reg=0):
  """Inverts a matrix given as a flat 4-list""";
  a,b,c,d = A;
  if reg:
    a = a+reg;
    d = d+reg;
  if not (is_null(b) and is_null(c)):
    det = a*d-b*c;
    return d/det,-b/det,-c/det,a/det;
  else:
    return 1/a,0,0,1/d;

def matrix_maskinf (A):
  """Replaces infinities with unity matrix""";
  mask = ~numpy.isfinite(A[0]);
  for a in A[1:]:
    mask |= ~numpy.isfinite(a);
  for a,defval in zip(A,(1,0,0,1)):
    if not numpy.isscalar(a):
      a[mask] = defval;
  return A
  
def array_to_vells (x,vellshape=None,expanded_slice=None):
  if expanded_slice is not None:
    a = meq.complex_vells(vellshape);
    a[...] = x[expanded_slice];
  else:
    a = meq.complex_vells(x.shape);
    a[...] = x;
  return a;

def mask_to_flags (x,vellshape=None,expanded_slice=None):
  if expanded_slice is not None:
    a = meq.flags(vellshape);
    a[...] = x[expanded_slice];
  else:
    a = meq.flags(x.shape);
    a[...] = x;
  return a;

def matrix_sqrt (A):
  """Returns the matrix square root of A""";
  a,b,c,d = A;
  # find eigenvalues
  r1,r2 = matrix_eigenval(A);
  sqrt_r1,sqrt_r2 = numpy.sqrt(r1),numpy.sqrt(r2);
  r2_r1 = r2-r1;
#  print r1,r2;
  if numpy.isscalar(r2_r1):
    if r1 != r2:
      # m:= (sqrt(r_2)-sqrt(r_1))/(r_2-r_1)
      # p:=(r_2*sqrt(r_1) - r_1*sqrt(r_2))/(r_2-r_1) .
      # sqrt(A) = m*A + p*I
      m = (sqrt_r2 - sqrt_r1)/r2_r1;
      p = (r2*sqrt_r1 - r1*sqrt_r2)/r2_r1;
    else:
      m = 1/(4*r1);
      p = sqrt_r1/2;
  else:
    m  = (sqrt_r2 - sqrt_r1)/r2_r1;
    p  = (r2*sqrt_r1 - r1*sqrt_r2)/r2_r1;
    m1 = 1/(4*r1);
    p1 = numpy.sqrt(r1)/2;
    wh = (r2_r1 == 0);
    m[wh] = m1[wh];
    p[wh] = p1[wh];

  return [a*m+p,b*m,c*m,d*m+p];

def matrix_qrd (A):
  """Returns the QR-decomposition of A.""";
  a,b,c,d = A;
  D  = a*b + c*d;
  x = numpy.sqrt(a*a+c*c);
  z = D/x;
  y = (a*d-b*c)/x;
  q11 = a/x;
  q21 = c/x;
  q12 = (b-q11*z)/y;
  q22 = (d-q21*z)/y;
  return (q11,q12,q21,q22),(x,z,0,y);

def matrix_eigenval (A):
  # eigenvalue decomposition, see
  # http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html
  a,b,c,d = A;
  D = a*b + c*d;
  T = a+d;
  T_2 = T/2;
  X = numpy.sqrt(T_2*T_2-D); #
  L1 = T_2 + X;
  L2 = T_2 - X;
  return L1,L2;

def matrix_eigenval_vec (A):
  # eigenvalue decomposition, see
  # http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html
  # see also: http://ieeexplore.ieee.org/xpl/freeabs_all.jsp?arnumber=486688
  L1,L2 = matrix_eigenval(A);
  a,b,c,d = A;
  if numpy.isscalar(L1):
    # scalar form
    if c<1e-20:
      if b<1e-20:
        return L1,L2,[1,0],[0,1];
      else:
        return L1,L2,[b,L1-a],[c,L1-a];
    else:
      return L1,L2,[L1-d,c],[L2-d,c];
  else:
    # array form
    wh1 = abs(b)<1e-20;
    wh2 = abs(c)<1e-20;
    wh12 = wh1&wh2;
    e1 = [L1-d,c];
    e2 = [L2-d,c];
    if wh2.sum():
      e1[0][wh2] = e2[0][wh2] = b[wh2];
      e1[1][wh2] = (L1-a)[wh2];
      e2[1][wh2] = (L2-a)[wh2];
    e1[0][wh12] = e2[1][wh12] = 1;
    e1[1][wh12] = e2[0][wh12] = 0;
  return L1,L2,e1,e2;
