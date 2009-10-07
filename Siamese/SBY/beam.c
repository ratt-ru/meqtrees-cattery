#include <math.h>
#include <complex.h>
double sqr(double x){
				  return x*x;
}

/* 
 * equation
 * (1-sin(th)^2 sin(phi)^2)sin(2*pi*f/c*h_0*cos(th))^2
 * c: speed of light, f : frequency h_0: height from ground
 * th: pi/2-elevation
 * phi: phi_0+azimuth, phi_0: dipole orientation
 * parms: scale, phi_0, h_0
 * axes: time,freq, az, el
 */
double test_double(const double *par,const double *x){
  const double c=3.e8;
  if (x[3]<=0) return 0; /* below horizon */
  const double theta=M_PI/2-x[3];
  const double phi=par[1]+x[2];
  const double result=(1-sqr(sin(theta))*sqr(sin(phi)))*sqr(sin(2*M_PI*x[1]*par[2]*cos(theta)/c));
  return (par[0]*sqrt(result));
}
complex test_complex(const complex *par,const complex *x){
  return 0;
}
int Npar_test=3;
int Nx_test=4;
