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
double complex YLMDsmallrppos[20][20] = {0};
double complex YLMAsmallrppos[20][20] = {0};
double complex YLMDsmallrpneg[20][20] = {0};
double complex YLMAsmallrpneg[20][20] = {0};
double complex FLMpos[20][20] = {0};
double complex FLMneg[20][20] = {0};
double complex ffunction[3][3]={0};
double epsilon = 80;
double radius = 4;
double depth = 10;
double smallr = 2;
double smallrp = 2;
int lmax=1;
double Dsmallr, Dbigr, Abigr=0;
double Dsmallrp, Asmallrp=0;  //Dsmallrp involves r-prime, and the small r is also the first argument of the image charge

void projectontoYLM();
void imagecharges(double, double, double, double);
void calculateFLM(void);

int main(){

  projectontoYLM();
  calculateFLM();
  return 0;
}


void projectontoYLM(){
  //at the moment this is just using the rectangle method (I'm pathetic, I know)
  int Ntheta=30;
  int Nphi=60;
  double theta, thetap, phi, phip=0;
  double plmplm=0;
  double complex exponentials=0;
  double complex ylmylm=0;
  double complex value=0;
  double h=0;
  int i,j,k,l,m,n=0;
  h = pow(PI, 4)*4/(Ntheta*Ntheta*Nphi*Nphi);
  printf("%f\n", h);
  for(l=0; l<=lmax; l++){
    for(m=0; m<=l; m++){
      for(i=0; i<=Ntheta-1; i++){
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
	      YLMDsmallrppos[l][m] += ylmylm*Dsmallrp;
	      YLMAsmallrppos[l][m] += ylmylm*Asmallrp;
	      plmplm = plmplm*pow(gsl_sf_fact(l-m), 2)/pow(gsl_sf_fact(l+m), 2);
	      exponentials = conj(exponentials);
	      ylmylm = plmplm*exponentials;
	      YLMAbigrneg[l][m] += ylmylm*Abigr;
	      YLMDbigrneg[l][m] += ylmylm*Dbigr;
	      YLMDsmallrneg[l][m] += ylmylm*Dsmallr;
	      YLMDsmallrpneg[l][m] += ylmylm*Dsmallrp;
	      YLMAsmallrpneg[l][m] += ylmylm*Asmallrp;
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
      YLMAsmallrppos[l][m] = YLMAsmallrppos[l][m]*h;
      YLMAsmallrpneg[l][m] = YLMAsmallrpneg[l][m]*h;
      YLMDsmallrppos[l][m] = YLMDsmallrppos[l][m]*h;
      YLMDsmallrpneg[l][m] = YLMDsmallrpneg[l][m]*h;
      
    }
  }

  printf("%f %f\n", creal(YLMAbigrpos[0][0]), h);
}

void calculateFLM(){
  int i,j,k,l,m=0;
  double complex numerator=0;
  double complex denominator=0;
  double complex energy=0;
  double complex exponentials;
  double complex ylmylm=0;
  double plmplm=0;
  double p4 =0;
  double lp21, eps = 0;
  eps=epsilon;
  p4 = 4*PI;
  for(l=0; l<=lmax; l++){
    lp21 = 2*l+1;
    for(m=0; m<=l; m++){
      numerator = (pow(eps-1, 3)/(p4*eps))*((pow(smallr/radius, l)*p4*(l+1)/(lp21))*(-p4/(radius*(lp21)) - YLMAbigrpos[l][m]) + radius*YLMDsmallrpos[l][m]*(p4/(lp21) + radius*YLMAbigrpos[l][m]));
      printf("num1 %f\n", creal(numerator));
      denominator = (eps-(eps-1))*(pow(radius, 2)*YLMDsmallrpos[l][m]-p4*(l+1)/(lp21));
      numerator = numerator*((-p4*(l+1)/(lp21))*pow(smallrp/radius, l) + pow(radius, 2)*YLMDsmallrppos[l][m]);
      printf("num2 %f\n", creal(numerator));
      FLMpos[l][m] = numerator/denominator;
      printf("1 %d %d %f %f %f %f\n", l, m, creal(FLMpos[l][m]), cimag(FLMpos[l][m]), creal(numerator), creal(denominator));
      numerator = (pow(eps-1, 2)/(p4*eps))*(-pow(smallrp/radius, l)*pow(smallr, l)*pow(radius,-(l+1))*pow(p4/lp21, 2)*(l+1));
      numerator +=  (pow(eps-1, 2)/(p4*eps))*pow(smallrp, l)*pow(radius, -l+1)*YLMDsmallrpos[l][m]*p4/lp21;
      numerator += (pow(eps-1, 2)/(p4*eps))*pow(radius, 2)* YLMAsmallrppos[l][m]*YLMDsmallrpos[l][m];
      FLMpos[l][m] += numerator;
      printf("2 %d %d %f %f\n", l, m, creal(FLMpos[l][m]), cimag(FLMpos[l][m]));
      numerator = (pow(eps-1, 3)/(p4*eps))*((pow(smallr/radius, l)*p4*(l+1)/(lp21))*(-p4/(radius*(lp21)) - YLMAbigrneg[l][m]) + radius*YLMDsmallrneg[l][m]*(p4/(lp21) + radius*YLMAbigrneg[l][m]));
      denominator = (eps-(eps-1))*(pow(radius, 2)*YLMDsmallrneg[l][m]-p4*(l+1)/(lp21));
      numerator = numerator*((-p4*(l+1)/(lp21))*pow(smallrp/radius, l) + pow(radius, 2)*YLMDsmallrpneg[l][m]);
      FLMneg[l][m] = numerator/denominator;
      numerator = (pow(eps-1, 2)/(p4*eps))*(-pow(smallrp/radius, l)*pow(smallr, l)*pow(radius,-(l+1))*pow(p4/lp21, 2)*(l+1));
      numerator +=  (pow(eps-1, 2)/(p4*eps))*pow(smallrp, l)*pow(radius, -l+1)*YLMDsmallrneg[l][m]*p4/lp21;
      numerator += (pow(eps-1, 2)/(p4*eps))*pow(radius, 2)* YLMAsmallrpneg[l][m]*YLMDsmallrneg[l][m];
      FLMneg[l][m] += numerator;

    }
  }

  for(l=0; l<=lmax; l++){
    m=0;
    plmplm = gsl_sf_legendre_sphPlm(l, 0, cos(0))*gsl_sf_legendre_sphPlm(l, 0, cos(0));
    exponentials = cexp(I*0*m)*cexp(-I*0*m);
    ylmylm = exponentials*plmplm;
    ffunction[0][0] += FLMpos[l][0]*ylmylm;
    printf("%f %f\n", creal(ffunction[0][0]), cimag(ffunction[0][0]));
    printf("%f %f %f %f\n", creal(FLMpos[l][0]), cimag(FLMpos[l][0]), creal(ylmylm), cimag(ylmylm));
    for(m=1; m<=l; m++){
       plmplm = gsl_sf_legendre_sphPlm(l, m, cos(0))*gsl_sf_legendre_sphPlm(l, m, cos(0));
       exponentials = cexp(I*0*m)*cexp(-I*0*m);
       ylmylm = exponentials*plmplm;
       ffunction[0][0] += FLMpos[l][0]*ylmylm;
       printf("%f %f\n", creal(ffunction[0][0]), cimag(ffunction[0][0]));
       plmplm = plmplm*pow(gsl_sf_fact(l-m), 2)/pow(gsl_sf_fact(l+m), 2);
       exponentials = conj(exponentials);
       ylmylm = exponentials*plmplm;
       ffunction[0][0] += FLMneg[l][m]*ylmylm;
       printf("%f %f\n", creal(ffunction[0][0]), cimag(ffunction[0][0]));
    }
  }
  printf("%f %f\n", creal(ffunction[0][0]), cimag(ffunction[0][0]));
}

void imagecharges(double theta, double phi, double thetap, double phip){ 
  double result=0;
  double inverser1, inverser2, inverser3=0;
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
  numerator = (smallr*cosp*sint-radius*cospp*sintp)*cospp*sintp + (smallr*sinp*sint-radius*sinpp*sintp)*sinpp*sintp - (smallr*cost+radius*costp-2*depth)*costp;
  Dsmallr = numerator*pow(inverser2, 3)*(epsilon-1)/(epsilon+1);
  numerator = -radius*(cosp*sint - cospp*sintp)*cosp*sint - radius*(sinp*sint - sinpp*sintp)*sinp*sint - (radius*cost + radius*costp - 2*depth)*cost;
  Dbigr = numerator*pow(inverser1, 3)*(epsilon-1)/(epsilon+1);

  inverser3 =  (pow(smallrp*cosp*sint - radius*cospp*sintp, 2) + pow(smallrp*sinp*sint - radius*sinpp*sintp,2) + pow(smallrp*cost-radius*costp-2*depth,2));
  inverser3 = 1/sqrt(inverser3);
  Asmallrp = inverser3*(epsilon-1)/(epsilon+1);
  numerator = -(smallrp*cosp*sint-radius*cospp*sintp)*cospp*sintp - (smallrp*sinp*sint-radius*sinpp*sintp)*sinpp*sintp - (smallrp*cost+radius*costp-2*depth)*costp;
  Dsmallrp = numerator*pow(inverser3, 3)*(epsilon-1)/(epsilon+1);
}
