#include <math.h>
#include <complex.h>
double sqr(double x){
				  return x*x;
}

/* 
 * not a beam, but cuts off source if elevation is negative 
 */
double test_double(const double *par,const double *x){
  if (x[3]<=0.0) return 0; /* below horizon */
  return(1.0);
}
complex test_complex(const complex *par,const complex *x){
  return 0;
}
int Npar_test=3;
int Nx_test=4;
