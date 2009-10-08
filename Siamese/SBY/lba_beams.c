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
   return C1;
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
   return C1;
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
   return C1;
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
   return C1;
}

/* 
 * equation - droopy dipole
 * equation: see writeup
 * c: speed of light, f : frequency
 * th: pi/2-elevation
 * phi: phi_0+azimuth, phi_0: dipole orientation
 * h: height of center from ground, L: projected arm length
 * alpha: droop angle

 * parameters: az,el,phi0,h,L,alpha
 * axes: time (ignored), freq
 */
complex double lba_theta_complex(const complex *par,const complex *x){
  const double c=3.e8;
  const double freq=creal(x[1]);
  const double az=creal(par[0]);
  const double el=creal(par[1]);
  const double phi0=creal(par[2]);
  const double h=creal(par[3]);
  const double L=creal(par[4]);
  const double alpha=creal(par[5]);

  if (el<=0.0) return (0+0*_Complex_I); /* below horizon */
  const double theta=M_PI/2-el;
  const double phi=phi0+az; /* take orientation into account */

  /* some essential constants */
  double k=2*M_PI*freq/c;

  /* calculate needed trig functions */
  double tan_al=tan(alpha);
  double sin_al=sin(alpha);
  double cos_al=cos(alpha);
  double sin_th=sin(theta);
  double cos_th=cos(theta);
  double sin_ph=sin(phi);
  double cos_ph=cos(phi);

  /* mu/4PI=10e-7  x omega/sin(alpha)*/
  const double A=(1e-7)*2*M_PI*freq/sin_al;

  complex double tmp=Gamma1(tan_al,sin_al,sin_th,cos_th,cos_ph,k,h,L);
  complex double Eth=A*tmp*(cos_al*sin_th-sin_al*cos_th*cos_ph);

  tmp=Gamma2(tan_al,sin_al,sin_th,cos_th,cos_ph,k,h,L);
  Eth+=A*tmp*(-cos_al*sin_th-sin_al*cos_th*cos_ph);

  tmp=Gamma3(tan_al,sin_al,sin_th,cos_th,cos_ph,k,h,L);
  Eth+=A*tmp*(-cos_al*sin_th+sin_al*cos_th*cos_ph);

  tmp=Gamma4(tan_al,sin_al,sin_th,cos_th,cos_ph,k,h,L);
  Eth+=A*tmp*(cos_al*sin_th+sin_al*cos_th*cos_ph);

  return (Eth);
}
int Npar_lba_theta=6;
int Nx_lba_theta=2;

complex double lba_phi_complex(const complex *par,const complex *x){
  const double c=3.e8;
  const double freq=creal(x[1]);
  const double az=creal(par[0]);
  const double el=creal(par[1]);
  const double phi0=creal(par[2]);
  const double h=creal(par[3]);
  const double L=creal(par[4]);
  const double alpha=creal(par[5]);

  if (el<=0.0) return (0+0*_Complex_I); /* below horizon */
  const double theta=M_PI/2-el;
  const double phi=phi0+az;

  /* some essential constants */
  double k=2*M_PI*freq/c;

  /* calculate needed trig functions */
  double tan_al=tan(alpha);
  double sin_al=sin(alpha);
  double cos_al=cos(alpha);
  double sin_th=sin(theta);
  double cos_th=cos(theta);
  double sin_ph=sin(phi);
  double cos_ph=cos(phi);

  /* mu/4PI=10e-7  x omega/sin(alpha)*/
  const double A=(1e-7)*2*M_PI*freq/sin_al;

  complex double tmp=Gamma1(tan_al,sin_al,sin_th,cos_th,cos_ph,k,h,L);
  complex double Eph=A*tmp*(sin_al*sin_ph);

  tmp=Gamma2(tan_al,sin_al,sin_th,cos_th,cos_ph,k,h,L);
  Eph+=A*tmp*(sin_al*sin_ph);

  tmp=Gamma3(tan_al,sin_al,sin_th,cos_th,cos_ph,k,h,L);
  Eph+=A*tmp*(-sin_al*sin_ph);

  tmp=Gamma4(tan_al,sin_al,sin_th,cos_th,cos_ph,k,h,L);
  Eph+=A*tmp*(-sin_al*sin_ph);

  return(Eph);

}
int Npar_lba_phi=6;
int Nx_lba_phi=2;

