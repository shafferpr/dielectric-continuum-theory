#include <stdio.h>
#include <gsl/gsl_sf_legendre.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.14159265
double YLMDsmallr[20][40] = {0};
double YLMDbigr[20][40] = {0};
double YLMAbigr[20][40] = {0};
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
  int i,j,k,l,m,n=0;
  for(i=0; i<=Ntheta-1; i++){
    printf("%d\n", i);
    theta = i*(PI)/Ntheta;
    for(j=0; j<=Nphi-1; j++){
      phi = j*(2*PI)/Nphi;
      for(k=0; k<=Ntheta-1; k++){
	thetap = k*(PI)/Ntheta;
	for(l=0; l<=Nphi-1; l++){
	  phip = l*(2*PI)/Nphi;
	  YLMAbigr[0][0] += sin(theta)*sin(thetap)*gsl_sf_legendre_sphPlm (0, 0, cos(theta))*gsl_sf_legendre_sphPlm (0, 0, cos(thetap))* Abigr(theta, phi, thetap, phip);
	}
      }
    }
  }
  printf("%f\n", YLMAbigr[0][0]);
}

double Abigr(double theta, double phi, double thetap, double phip){
  double result=0;
  result = pow(cos(phi)*sin(theta) - cos(phip)*sin(thetap), 2) + pow(sin(phi)*sin(theta)-sin(phip)*sin(thetap), 2) + pow(cos(theta) + cos(thetap)-2*depth/radius, 2);
  result = (epsilon-1)/(sqrt(result)*radius*(epsilon+1));
  return result;
}
