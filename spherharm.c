#include <stdio.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gamma.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.14159265

double complex YLMDsmallrpos[20][20] = {0};
double complex YLMDbigrpos[20][20] = {0};
double complex YLMAbigrpos[20][20] = {0};
double complex YLMDsmallrneg[20][20] = {0};
double complex YLMDbigrneg[20][20] = {0};
double complex YLMAbigrneg[20][20] = {0};
double epsilon = 80;
double radius = 4;
double depth = 10;
double smallr = 2;
double Dsmallr, Dbigr, Abigr=0;
void projectontoYLM();
void imagecharges(double, double, double, double);

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
  double complex value=0;
  double h=0;
  int i,j,k,l,m,n=0;
  h = pow(PI, 4)*4/(Ntheta*Ntheta*Nphi*Nphi);
  printf("%f\n", h);
  for(l=0; l<=4-1; l++){
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
	      imagecharges(theta, phi, thetap, phip);
	      YLMAbigrpos[l][m] += ylmylm*Abigr;
	      YLMDbigrpos[l][m] += ylmylm*Dbigr;
	      YLMDsmallrpos[l][m] += ylmylm*Dsmallr;
	      plmplm = plmplm*pow(gsl_sf_fact(l-m), 2)/pow(gsl_sf_fact(l+m), 2);
	      exponentials = conj(exponentials);
	      ylmylm = plmplm*exponentials;
	      YLMAbigrneg[l][m] += ylmylm*Abigr;
	      YLMDbigrneg[l][m] += ylmylm*Dbigr;
	      YLMDsmallrneg[l][m] += ylmylm*Dsmallr;
	    }
	  }
	}
      }
      YLMAbigrpos[l][m] = YLMAbigrpos[l][m]*h;
      YLMAbigrneg[l][m] = YLMAbigrneg[l][m]*h;
      YLMDbigrpos[l][m] = YLMDbigrpos[l][m]*h;
      YLMDbigrneg[l][m] = YLMDbigrneg[l][m]*h;
      YLMDsmallrpos[l][m] = YLMDsmallrpos[l][m]*h;
      YLMDsmallrneg[l][m] = YLMDsmallrneg[l][m]*h;
    }
  }

  printf("%f %f\n", cimag(YLMAbigrpos[3][0]), h);
}

/*
double Abigr(double theta, double phi, double thetap, double phip){
  double result=0;
  result = pow(cos(phi)*sin(theta) - cos(phip)*sin(thetap), 2) + pow(sin(phi)*sin(theta)-sin(phip)*sin(thetap), 2) + pow(cos(theta) + cos(thetap)-2*depth/radius, 2);
  result = (epsilon-1)/(sqrt(result)*radius*(epsilon+1));
  return result;
}
*/
/*
double Dbigr(double theta, double phi, double thetap, double phip){
  double result=0;
  result = pow(cos(phi)*sin(theta) - cos(phip)*sin(thetap), 2) + pow(sin(phi)*sin(theta)-sin(phip)*sin(thetap), 2) + pow(cos(theta) + cos(thetap)-2*depth/radius, 2);
  result = (epsilon-1)/(sqrt(result)*radius*(epsilon+1));
  return result;
}
*/
void imagecharges(double theta, double phi, double thetap, double phip){ 
  double result=0;
  double inverser1, inverser2=0;
  double cost, costp, cosp, cospp=0;
  double sint, sintp, sinp, sinpp=0;
  double numerator=0;
  cost=cos(theta);
  cosp=cos(phi);
  costp=cos(thetap);
  cospp=cos(phip);
  sint=sin(theta);
  sinp=sin(phi);
  sintp=sin(thetap);
  sinpp=sin(phip);

  inverser1 = pow(cosp*sint - cospp*sintp, 2) + pow(sinp*sint-sinpp*sintp, 2) + pow(cost + costp - 2*depth/radius, 2);
  inverser1 = 1/(sqrt(inverser1)*radius);
  Abigr = inverser1*(epsilon-1)/(epsilon+1);

  inverser2 = (pow(smallr*cosp*sint - radius*cospp*sintp, 2) + pow(smallr*sinp*sint - radius*sinpp*sintp,2) + pow(smallr*cost-radius*costp-2*depth,2));
  inverser2 = 1/sqrt(inverser2);
  numerator = (smallr*cosp*sint-radius*cospp*sintp)*cospp*sintp + (radius*sinp*sint-radius*sinpp*sintp)*sinpp*sintp - (smallr*cost+radius*costp-2*depth)*costp;
  Dsmallr = numerator*pow(inverser2, 3)*(epsilon-1)/(epsilon+1);

  numerator = -radius*(cosp*sint - cospp*sintp)*cosp*sint - radius*(sinp*sint - sinpp*sintp)*sinp*sint - (radius*cost + radius*costp - 2*depth)*cost;
  Dbigr = numerator*pow(inverser1, 3)*(epsilon-1)/(epsilon+1);
}
