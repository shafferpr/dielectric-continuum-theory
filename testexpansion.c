#include <stdio.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gamma.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.14159265
int lmax=0;

void calcexpansion();

int main(int argc, char *argv[]){
  sscanf(argv[1], "%d", &lmax);
  printf("%d\n", lmax);
  calcexpansion();
  return 0;
}

void calcexpansion(){

  signed int i,j,k,l,m=0;
  double complex numerator=0;
  double complex denominator=0;
  double complex energy=0;
  double complex exponentials;
  double complex ylmylm=0;
  double plmplm=0;
  double theta, thetap=0;
  double phi, phip=0;
  double complex expansion1=0;
  double complex expansion2=0;
  float dum, fuckthis=0;
  theta = 0.0;
  thetap = 0;
  phi = 0.0;
  phip = PI/2;

  for(l=0; l<=lmax; l++){
    m=0;
    plmplm = gsl_sf_legendre_sphPlm(l, 0, cos(theta))*gsl_sf_legendre_sphPlm(l, 0, cos(thetap));
    exponentials = cexp(I*phi*m)*cexp(-I*phip*m);
    ylmylm = exponentials*plmplm;
    dum=l;
    expansion1 += ylmylm*(1/(2*dum+1));
    expansion2 += ylmylm*(1/(2*dum+1));
    printf("expansion1 %f\n", creal(expansion1));
    //printf("ylms %f %f %d %f\n", creal(ylmylm), cimag(ylmylm), l, fuckthis);
    //printf("e1 %f %f\n", creal(expansion1), cimag(expansion1));
    //printf("e2 %f %f\n", creal(expansion2), cimag(expansion2));
    for(m=1; m<=l; m++){
       plmplm = gsl_sf_legendre_sphPlm(l, m, cos(theta))*gsl_sf_legendre_sphPlm(l, m, cos(thetap));
       exponentials = cexp(I*phi*m)*cexp(-I*phip*m);
       ylmylm = exponentials*plmplm;
       //printf("pos %f %f %f %f %f\n", creal(ylmylm), cimag(ylmylm), creal(exponentials), cimag(exponentials), plmplm);
       expansion1 += ylmylm*(1/(2*dum+1));
       expansion2 += conj(exponentials)*plmplm*(1/(2*dum+1));
       plmplm = plmplm*pow(gsl_sf_fact(l-m), 2)/pow(gsl_sf_fact(l+m), 2);
       exponentials = conj(exponentials);
       ylmylm = conj(ylmylm);
       //printf("neg %f %f %f %f %f\n", creal(ylmylm), cimag(ylmylm), creal(exponentials), cimag(exponentials), plmplm);
       expansion1 += ylmylm*(1/(2*dum+1));
       //printf("%f\n", creal(expansion1));
       expansion2 += conj(exponentials)*plmplm*(1/(2*dum+1));
       //printf("%d %d\n", l, m);
       }
  }

  printf("e1 %f %f\n", 4*PI*creal(expansion1), 4*PI*cimag(expansion1));
  printf("e2 %f %f\n", 4*PI*creal(expansion2), 4*PI*cimag(expansion2));
}
