#include <math.h>
#include <complex.h>

/* see writeup for the exact formula */
inline complex double 
Gamma1(double tan_al,double sin_al, double sin_th, double cos_th,
   double cos_ph, double k, double h, double L) {

   double A=sin_th*cos_ph-cos_th/tan_al;
   double inv_sin_al=1/sin_al;
   double tmp=k*A*L;
   double Re1=cos(tmp)*inv_sin_al;
   double Im1=sin(tmp)*inv_sin_al;

   Re1 =Re1-inv_sin_al*cos(k*L*inv_sin_al);
   Im1 =Im1-A*sin(k*L*inv_sin_al);

   complex double C1=Re1+Im1*_Complex_I;

   complex double C2=cos(k*h*cos_th)+sin(k*h*cos_th)*_Complex_I;
   
   C1=-C1*C2/(A*A-inv_sin_al*inv_sin_al);
   return (C1);
}

/* see writeup for the exact formula */
inline complex double 
Gamma2(double tan_al,double sin_al, double sin_th, double cos_th,
   double cos_ph, double k, double h, double L) {

   double A=-sin_th*cos_ph-cos_th/tan_al;
   double inv_sin_al=1/sin_al;
   double tmp=k*A*L;
   double Re1=cos(tmp)*inv_sin_al;
   double Im1=-sin(tmp)*inv_sin_al;

   Re1 =Re1-inv_sin_al*cos(k*L*inv_sin_al);
   Im1 =Im1+A*sin(k*L*inv_sin_al);

   complex double C1=Re1+Im1*_Complex_I;

   complex double C2=cos(k*h*cos_th)+sin(k*h*cos_th)*_Complex_I;
   
   C1=-C1*C2/(A*A-inv_sin_al*inv_sin_al);
   return (C1);
}

/* see writeup for the exact formula */
inline complex double 
Gamma3(double tan_al,double sin_al, double sin_th, double cos_th,
   double cos_ph, double k, double h, double L) {

   double A=-sin_th*cos_ph+cos_th/tan_al;
   double inv_sin_al=1/sin_al;
   double tmp=k*A*L;
   double Re1=cos(tmp)*inv_sin_al;
   double Im1=-sin(tmp)*inv_sin_al;

   Re1 =Re1-inv_sin_al*cos(k*L*inv_sin_al);
   Im1 =Im1+A*sin(k*L*inv_sin_al);

   complex double C1=Re1+Im1*_Complex_I;

   complex double C2=cos(k*h*cos_th)-sin(k*h*cos_th)*_Complex_I;
   
   C1=-C1*C2/(A*A-inv_sin_al*inv_sin_al);
   return (C1);
}

/* see writeup for the exact formula */
inline complex double 
Gamma4(double tan_al,double sin_al, double sin_th, double cos_th,
   double cos_ph, double k, double h, double L) {

   double A=sin_th*cos_ph+cos_th/tan_al;
   double inv_sin_al=1/sin_al;
   double tmp=k*A*L;
   double Re1=cos(tmp)*inv_sin_al;
   double Im1=sin(tmp)*inv_sin_al;

   Re1 =Re1-inv_sin_al*cos(k*L*inv_sin_al);
   Im1 =Im1-A*sin(k*L*inv_sin_al);

   complex double C1=Re1+Im1*_Complex_I;

   complex double C2=cos(k*h*cos_th)-sin(k*h*cos_th)*_Complex_I;
   
   C1=-C1*C2/(A*A-inv_sin_al*inv_sin_al);
   return (C1);
}


/* 
 * equation - droopy dipole
 * equation: see writeup
 * c: speed of light, f : frequency
 * th: pi/2-elevation
 * phi: phi_0+azimuth, phi_0: dipole orientation
 * parameters: h,L,alpha,phi_0
 * h: height of center from ground, L: projected arm length
 * alpha: droop angle
 * axes: time,freq, az, el
 */
double test_double(const double *par,const double *x){
  const double c=3.e8;
  if (x[3]<=0.0) return 0; /* below horizon */
  const double theta=M_PI/2-x[3];
  const double phi=par[3]+x[2];

  /* some essential constants */
  double k=2*M_PI*x[1]/c;

  /* calculate needed trig functions */
  double tan_al=tan(par[2]);
  double sin_al=sin(par[2]);
  double cos_al=cos(par[2]);
  double sin_th=sin(theta);
  double cos_th=cos(theta);
  double sin_ph=sin(phi);
  double cos_ph=cos(phi);

  /* mu/4PI=10e-7  x omega/sin(alpha)*/
  const double A=(1e-7)*2*M_PI*x[1]/sin_al;

  complex double tmp=Gamma1(tan_al,sin_al,sin_th,cos_th,cos_ph,k,par[0],par[1]);
  complex double Eph=A*tmp*(sin_al*sin_ph);

  tmp=Gamma2(tan_al,sin_al,sin_th,cos_th,cos_ph,k,par[0],par[1]);
  Eph+=A*tmp*(sin_al*sin_ph);

  tmp=Gamma3(tan_al,sin_al,sin_th,cos_th,cos_ph,k,par[0],par[1]);
  Eph+=A*tmp*(-sin_al*sin_ph);

  tmp=Gamma4(tan_al,sin_al,sin_th,cos_th,cos_ph,k,par[0],par[1]);
  Eph+=A*tmp*(-sin_al*sin_ph);

  return(creal(Eph));
}

complex double test_complex(const complex *par,const complex *x){
  const double c=3.e8;
  const double x1=creal(x[1]);
  const double x2=creal(x[2]);
  const double x3=creal(x[3]);
  const double p0=creal(par[0]);
  const double p1=creal(par[1]);
  const double p2=creal(par[2]);
  const double p3=creal(par[3]);
  if (x3<=0.0) return (0+0*_Complex_I); /* below horizon */
  const double theta=M_PI/2-x3;
  const double phi=p3+x2;

  /* some essential constants */
  double k=2*M_PI*x1/c;

  /* calculate needed trig functions */
  double tan_al=tan(p2);
  double sin_al=sin(p2);
  double cos_al=cos(p2);
  double sin_th=sin(theta);
  double cos_th=cos(theta);
  double sin_ph=sin(phi);
  double cos_ph=cos(phi);

  /* mu/4PI=10e-7  x omega/sin(alpha)*/
  const double A=(1e-7)*2*M_PI*x1/sin_al;

  complex double tmp=Gamma1(tan_al,sin_al,sin_th,cos_th,cos_ph,k,p0,p1);
  complex double Eph=A*tmp*(sin_al*sin_ph);

  tmp=Gamma2(tan_al,sin_al,sin_th,cos_th,cos_ph,k,p0,p1);
  Eph+=A*tmp*(sin_al*sin_ph);

  tmp=Gamma3(tan_al,sin_al,sin_th,cos_th,cos_ph,k,p0,p1);
  Eph+=A*tmp*(-sin_al*sin_ph);

  tmp=Gamma4(tan_al,sin_al,sin_th,cos_th,cos_ph,k,p0,p1);
  Eph+=A*tmp*(-sin_al*sin_ph);

  return(Eph);

}
int Npar_test=4;
int Nx_test=4;
