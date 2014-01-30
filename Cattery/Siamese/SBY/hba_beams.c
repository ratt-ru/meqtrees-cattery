#include <math.h>
#include <complex.h>

/* Note: formulas have been optimized! */
/* alpha1=50/180*pi;
 * alpha2=80/180*pi;
 * L=0.366;
 * h=45/100;
 */
#define palpha1  0.87266462599716
#define palpha2  1.39626340159546
#define sin_alpha1 0.76604444311898
#define sin_alpha2 0.98480775301221
#define cos_alpha1 0.64278760968654
#define cos_alpha2 0.17364817766693
#define tan_alpha1 1.19175359259421 
#define tan_alpha2 5.67128181961771
#define pH 0.450
#define pL 0.366
#define pL1  0.13275660600238
#define pL2  0.23888953394781
#define pC 2.998e8 /* speed of light */
#define TOL 1e-9
#define del1 0.47777906789561
#define del2 0.37164613995018

/* see writeup for the exact formula */
inline complex double
Gammplus(double A, double B, double k, double L, double alpha) {
  double salpha=1/sin(alpha);
  double tmp=k*B;
  complex double z=-(cos(tmp)+_Complex_I*sin(tmp))/(A*A-salpha*salpha+TOL);
  tmp=k*A*L;
  double ang1=k*(L)*salpha;
  complex double part1=cos(tmp)+_Complex_I*sin(tmp);
  part1=part1*salpha-_Complex_I*A*sin(ang1)-salpha*cos(ang1);
  return z*part1;
}
inline complex double
Gammminus(double A, double B, double k, double L, double alpha) {
  double salpha=1/sin(alpha);
  double tmp=k*B;
  complex double z=-(cos(tmp)+_Complex_I*sin(tmp))/(A*A-salpha*salpha+TOL);
  tmp=-k*A*L;
  double ang1=k*(L)*salpha;
  complex double part1=cos(tmp)+_Complex_I*sin(tmp);
  part1=part1*salpha+_Complex_I*A*sin(ang1)-salpha*cos(ang1);
  return z*part1;
}
inline complex double
Gammplus0(double A, double B, double k, double L) {
  double tmp=k*B;
  complex double z=-(cos(tmp)+_Complex_I*sin(tmp))/(A*A-1.0+TOL);
  tmp=k*A*L;
  double ang1=k*(L);
  complex double part1=cos(tmp)+_Complex_I*sin(tmp);
  part1=part1-_Complex_I*A*sin(ang1)-cos(ang1);
  return z*part1;
}
inline complex double
Gammminus0(double A, double B, double k, double L) {
  double tmp=k*B;
  complex double z=-(cos(tmp)+_Complex_I*sin(tmp))/(A*A-1.0+TOL);
  tmp=-k*A*L;
  double ang1=k*(L);
  complex double part1=cos(tmp)+_Complex_I*sin(tmp);
  part1=part1+_Complex_I*A*sin(ang1)-cos(ang1);
  return z*part1;
}

/* 
 * equation - HBA bowtie
 * equation: see writeup

 * parameters: az,el,phi_0
 * axes: time,freq
 */

complex double hba_phi_complex(const complex *par,const complex *x){
  const double freq=creal(x[1]);
  const double az=creal(par[0]);
  const double el=creal(par[1]);
  const double phi0=creal(par[2]);

  if (el<=0.0) return (0+0*_Complex_I); /* below horizon */
  const double theta=M_PI_2-el;
  const double phi=phi0+az; /* take orientation into account */

  /* some essential constants */
  double k=2*M_PI*freq/pC;

  /* calculate needed trig functions */
  double sin_theta=sin(theta);
  double cos_theta=cos(theta);
  double sin_phi=sin(phi);
  double cos_phi=cos(phi);

  /* mu/4PI=10e-7  x omega*/
  const double mop=(1e-7)*2*M_PI*freq;

  complex double Phi11=Gammplus(sin_theta*cos_phi-cos_theta/tan_alpha1,pH*cos_theta,k,pL,palpha1);
  complex double Eph=sin_alpha1*sin_phi*Phi11;
  Phi11=Gammplus(sin_theta*cos_phi+cos_theta/tan_alpha2,pH*cos_theta,k,pL,M_PI-palpha2);
  Eph+=sin_alpha2*sin_phi*Phi11;


  Phi11=Gammminus(-sin_theta*cos_phi+cos_theta/tan_alpha1,pH*cos_theta,k,pL,palpha1);
  Eph+=sin_alpha1*sin_phi*Phi11;
  Phi11=Gammminus(-sin_theta*cos_phi-cos_theta/tan_alpha2,pH*cos_theta,k,pL,M_PI-palpha2);
  Eph+=sin_alpha2*sin_phi*Phi11;


  Phi11=Gammminus(-sin_theta*cos_phi-cos_theta/tan_alpha1,-pH*cos_theta,k,pL,palpha1);
  Eph+=-sin_alpha1*sin_phi*Phi11;
  Phi11=Gammminus(-sin_theta*cos_phi+cos_theta/tan_alpha2,-pH*cos_theta,k,pL,M_PI-palpha2);
  Eph+=-sin_alpha2*sin_phi*Phi11;


  Phi11=Gammplus(sin_theta*cos_phi+cos_theta/tan_alpha1,-pH*cos_theta,k,pL,palpha1);
  Eph+=-sin_alpha1*sin_phi*Phi11;
  Phi11=Gammplus(sin_theta*cos_phi-cos_theta/tan_alpha2,-pH*cos_theta,k,pL,M_PI-palpha2);
  Eph+=-sin_alpha2*sin_phi*Phi11;

  return (Eph*mop);
}
int Npar_hba_phi=3;
int Nx_hba_phi=2;


complex double hba_theta_complex(const complex *par,const complex *x){
  const double freq=creal(x[1]);
  const double az=creal(par[0]);
  const double el=creal(par[1]);
  //const double p0=creal(par[0]);
 // const double p1=creal(par[1]);
  //const double p2=creal(par[2]);
  const double phi0=creal(par[2]);

  if (el<=0.0) return (0+0*_Complex_I); /* below horizon */
  const double theta=M_PI_2-el;
  const double phi=phi0+az; /* take orientation into account */

  /* some essential constants */
  double k=2*M_PI*freq/pC;

  /* extra params for length of current distribution */
  double wt=(240-freq/1e6)/140.0;
  double pLL1=wt*pL1+(1-wt)*pL2;
  double pLL2=wt*pL2+(1-wt)*pL1;

  /* calculate needed trig functions */
  double sin_theta=sin(theta);
  double cos_theta=cos(theta);
  //double sin_phi=sin(phi);
  double cos_phi=cos(phi);

  /* mu/4PI=10e-7  x omega*/
  const double mop=(1e-7)*2*M_PI*freq;

  complex double expjkL1=sin(k*del1)+_Complex_I*cos(k*del1);
  complex double expjkL1_=sin(k*del1)-_Complex_I*cos(k*del1);
  complex double expjkL2=sin(k*del2)+_Complex_I*cos(k*del2);
  complex double expjkL2_=sin(k*del2)-_Complex_I*cos(k*del2);

  complex double Phi11=Gammplus(sin_theta*cos_phi-cos_theta/tan_alpha1,pH*cos_theta,k,pL,palpha1);
  complex double Eth=(-cos_alpha1*sin_theta-sin_alpha1*cos_theta*cos_phi)*Phi11;
  Phi11=Gammplus(sin_theta*cos_phi+cos_theta/tan_alpha2,pH*cos_theta,k,pL,M_PI-palpha2);
  Eth+=(cos_alpha2*sin_theta-sin_alpha2*cos_theta*cos_phi)*Phi11;
  Phi11=expjkL1*Gammplus0(cos_theta,pL*sin_theta*cos_phi+(pH-pL/tan_alpha1)*cos_theta,k,pLL1);
  Eth+=sin_theta*Phi11;
  Phi11=-expjkL2*Gammplus0(-cos_theta,pL*sin_theta*cos_phi+(pH+pL/tan_alpha2)*cos_theta,k,pLL2);
  Eth+=sin_theta*Phi11;


  Phi11=Gammminus(-sin_theta*cos_phi+cos_theta/tan_alpha1,pH*cos_theta,k,pL,palpha1);
  Eth+=(cos_alpha1*sin_theta-sin_alpha1*cos_theta*cos_phi)*Phi11;
  Phi11=Gammminus(-sin_theta*cos_phi-cos_theta/tan_alpha2,pH*cos_theta,k,pL,M_PI-palpha2);
  Eth+=(-cos_alpha2*sin_theta-sin_alpha2*cos_theta*cos_phi)*Phi11;
  Phi11=-expjkL1_*Gammplus0(cos_theta,-pL*sin_theta*cos_phi+(pH-pL/tan_alpha1)*cos_theta,k,pLL1);
  Eth+=sin_theta*Phi11;
  Phi11=expjkL2_*Gammplus0(-cos_theta,-pL*sin_theta*cos_phi+(pH+pL/tan_alpha2)*cos_theta,k,pLL2);
  Eth+=sin_theta*Phi11;


  Phi11=Gammminus(-sin_theta*cos_phi-cos_theta/tan_alpha1,-pH*cos_theta,k,pL,palpha1);
  Eth+=(cos_alpha1*sin_theta+sin_alpha1*cos_theta*cos_phi)*Phi11;
  Phi11=Gammminus(-sin_theta*cos_phi+cos_theta/tan_alpha2,-pH*cos_theta,k,pL,M_PI-palpha2);
  Eth+=(-cos_alpha2*sin_theta+sin_alpha2*cos_theta*cos_phi)*Phi11;
  Phi11=-expjkL1*Gammminus0(cos_theta,-pL*sin_theta*cos_phi-(pH-pL/tan_alpha1)*cos_theta,k,pLL1);
  Eth+=sin_theta*Phi11;
  Phi11=expjkL2*Gammminus0(-cos_theta,-pL*sin_theta*cos_phi-(pH+pL/tan_alpha2)*cos_theta,k,pLL2);
  Eth+=sin_theta*Phi11;


  Phi11=Gammplus(sin_theta*cos_phi+cos_theta/tan_alpha1,-pH*cos_theta,k,pL,palpha1);
  Eth+=(-cos_alpha1*sin_theta+sin_alpha1*cos_theta*cos_phi)*Phi11;
  Phi11=Gammplus(sin_theta*cos_phi-cos_theta/tan_alpha2,-pH*cos_theta,k,pL,M_PI-palpha2);
  Eth+=(cos_alpha2*sin_theta+sin_alpha2*cos_theta*cos_phi)*Phi11;
  Phi11=expjkL1_*Gammminus0(cos_theta,pL*sin_theta*cos_phi-(pH-pL/tan_alpha1)*cos_theta,k,pLL1);
  Eth+=sin_theta*Phi11;
  Phi11=-expjkL2_*Gammminus0(-cos_theta,pL*sin_theta*cos_phi-(pH+pL/tan_alpha2)*cos_theta,k,pLL2);
  Eth+=sin_theta*Phi11;

  return (Eth*mop);
}
int Npar_hba_theta=3;
int Nx_hba_theta=2;
