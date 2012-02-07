#include <stdio.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gamma.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.14159265
double YLMDsmallr[20][20][2] = {0};
double YLMDbigr[20][20][2] = {0};
double YLMAbigr[20][20][2] = {0};
double epsilon = 80;
double radius = 4;
double depth = 10;
void projectontoYLM();
double Abigr(double, double, double, double);

int main(){
  double x = 5.0;
  double y = gsl_sf_legendre_sphPlm (1, 1, .5);
  //printf("J0(%f) = %.18f\n", x, y);
  projectontoYLM();
  return 0;
}


void projectontoYLM(){
  //at the moment this is just using the rectangle method (I'm pathetic, I know)
  int Ntheta=50;
  int Nphi=100;
  double theta, thetap, phi, phip=0;
  double plmplm=0;
  double complex exponentials=0;
  double complex ylmylm=0;
  double h=0;
  int i,j,k,l,m,n=0;
  h = pow(PI, 4)*4/(Ntheta*Ntheta*Nphi*Nphi);
  for(l=0; l<=5-1; l++){
    for(m=0; m<=l; m++){
      for(i=0; i<=Ntheta-1; i++){
	printf("%d\n", i);
	theta = i*(PI)/Ntheta;
	for(j=0; j<=Nphi-1; j++){
	  phi = j*(2*PI)/Nphi;
	  for(k=0; k<=Ntheta-1; k++){
	    thetap = k*(PI)/Ntheta;
	    for(n=0; n<=Nphi-1; n++){
	      phip = n*(2*PI)/Nphi;
	      plmplm = sin(theta)*sin(thetap)*gsl_sf_legendre_sphPlm(l, m, cos(theta))*gsl_sf_legendre_sphPlm(l, m, cos(thetap));
	      exponentials = cexp(-I*phi*m)*cexp(I*phip*m);
	      ylmylm = plmplm*exponentials;
	      YLMAbigr[l][m][0] += ylmylm*Abigr(theta, phi, thetap, phip);
	      YLMDbigr[l][m][0] += ylmylm*Dbigr(theta, phi, thetap, phip);
	      YLMDsmallr[l][m][0] += ylmylm*Dsmallr(theta, phi, thetap, phip);
	      plmplm = plmplm*pow(gsl_sf_fact(l-m), 2)/pow(gsl_sf_fact(l+m), 2);
	      exponentials = cexp(I*phi*m)*cexp(-I*phip*m);
	      ylmylm = plmplm*exponentials;
	      YLMAbigr[l][m][1] += ylmylm*Abigr(theta, phi, thetap, phip);
	      YLMDbigr[l][m][1] += ylmylm*Dbigr(theta, phi, thetap, phip);
	      YLMDsmallr[l][m][1] += ylmylm*Dsmallr(theta, phi, thetap, phip);
	    }
	  }
	}
      }
      YLMAbigr[l][m][0] = YLMAbigr[l][m][0]*h;
      YLMAbigr[l][m][1] = YLMAbigr[l][m][1]*h;
      YLMDbigr[l][m][0] = YLMDbigr[l][m][0]*h;
      YLMDbigr[l][m][1] = YLMDbigr[l][m][1]*h;
      YLMDsmallr[l][m][0] = YLMDsmallr[l][m][0]*h;
      YLMDsmallr[l][m][1] = YLMDsmallr[l][m][1]*h;
    }
  }
  printf("%f %f\n", YLMAbigr[0][0], h);
}


double Abigr(double theta, double phi, double thetap, double phip){
  double result=0;
  result = pow(cos(phi)*sin(theta) - cos(phip)*sin(thetap), 2) + pow(sin(phi)*sin(theta)-sin(phip)*sin(thetap), 2) + pow(cos(theta) + cos(thetap)-2*depth/radius, 2);
  result = (epsilon-1)/(sqrt(result)*radius*(epsilon+1));
  return result;
}
