#include <stdio.h>
#include <gsl/gsl_sf_ellint.h>
#include <math.h>
#include <stdlib.h>
#include <armadillo>
#include <complex.h>
using namespace arma;
#define PI 3.14159265
//gcc -O2 -lgsl -lgslcblas shwatermolecule.c -o shwatermolecule
//g++ -O2 -lgsl -lgslcblas -llapack ellipticintegrals.c -o ellipticint
double result=0;
double B[200][200]={0};
doulbe K[200]={0};
double delta=0.01;
double epsilon=80;
double L, d=0;
int gridpoints=10;
int i,j,k=0;
void computeb();
void computek();
double Integral(double, double);


int main(){
  double input=0;
  computeb();
  computek();
  mat A(gridpoints, gridpoints);
  for(i=0; i<=gridpoints-1; i++){
    for(j=0; j<=gridpoints-1; j++){
      A(i,j)=B[i][j];
    }
  }
  mat C = inv(A);
 
  return 0;
}


void computeb(){
  double thetai, thetaj=0;
  double prefactor=0;
  double B[200][200]={0};
  prefactor = (epsilon-1)/(gridpoints*2*sqrt(2)*(epsilon+1));
  for(i=0; i<=gridpoints-1; i++){
    for(j=0; j<=gridpoints-1; j++){
      thetai=3.14159*i/gridpoints;
      thetaj=3.14159*j/gridpoints;
      if(i==j)
	B[i][j]=1;
      B[i][j] += prefactor*Integral(thetai, thetaj)*sin(thetai);
      printf("%d %d %f\n", i, j, B[i][j]);
    }
  }

}

void computek(){
  double spacing=0;
  double theta=0;
  double thetap=0;
  double sum=0;
  for(i=0; i<=gridpoints-1; i++){
    spacing=pi/gridpoints;
    theta=i*spacing;
    K[i] = (epsilon-1)*(epsilon-1)/(4*PI*epsilon*(epsilon+1)*(L*L-4*d*L*cos(theta)+4*d*d));
    sum=0;
    for(j=0; j<=gridpoints-1; j++){
      thetap = j*spacing;
      sum += Integral(theta, thetap)*(sin(thetap)/L + (epsilon-1)*(L-2*d*cos(thetap))*pow(L*L-4*d*l+4*d*d, -1.5));
      sum = sum*(epsilon-1/(4*PI),2)/(epsilon*sqrt(2));
    }
    K[i] -=sum;
  }

}

double Integral(double thetai, double thetaj){
  double A, B=0;
  double argument=0;
  double integral=0;
  A=1-cos(thetai)*cos(thetaj)+delta;
  B=sin(thetai)*sin(thetaj);
  argument = -2*B/(A-B);
  integral = 2*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A-B));
  argument = 2*B/(A+B);
  integral += 2*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A+B));

  return integral;
}

