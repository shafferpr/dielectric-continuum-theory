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
//this file computes solvation energy for d<L(the sphere penetrating the boundary)
double result=0;
double B[400][400]={0};
double K[400]={0};
double ftheta[400]={0};
double delta=0.01;
double epsilon=80;
double L, d=0;
double startingtheta=0;
int gridpoints=10;
int i,j,k=0;
void computeb();
void computek();
void computef();
double Integral(double, double);

int main(int argc, char *argv[]){
  sscanf(argv[1], "%lf", &L);
  sscanf(argv[2], "%lf", &d);
  sscanf(argv[3], "%d", &gridpoints);
  sscanf(argv[4], "%lf", &delta);
  double spacing=0;
  computeb();
  computek();
  mat A(gridpoints, gridpoints);
  startingtheta=0;
  spacing=(PI-startingtheta)/gridpoints;
  for(i=0; i<=gridpoints-1; i++){
    for(j=0; j<=gridpoints-1; j++){
      A(i,j)=B[i][j];
    }
  }
  vec k(gridpoints);
  for(i=0; i<=gridpoints-1; i++){
    k(i)=K[i];
  }
  //k.print();
  mat C = inv(A);
  vec f = C*k;
  f=f*2/(epsilon+1);
  for(i=0; i<=gridpoints-1; i++){
    ftheta[i]=f(i);
  }

  computef();
  return 0;
}

void computef(){
  double integral=0;
  double fdd=0;
  double spacing=0;
  double thetai=0;
  double energy=0;
  spacing=(PI-startingtheta)/gridpoints;
  fdd = -(pow(epsilon-1,2))/(4*PI*epsilon*L);
  //printf("%f\n", fdd);
  integral=0;
  for(i=0; i<=gridpoints-1; i++){
    thetai=startingtheta + (3.14159-startingtheta)*i/(gridpoints);
    integral += sin(thetai)*spacing*L*(L-2.0*d*cos(thetai))/pow(L*L-4.0*d*L*cos(thetai)+4.0*d*d,1.5);
  }
  integral = integral*pow(epsilon-1,2)/(8*PI*(epsilon+1));
  fdd -=integral;
  integral=0;
  //printf("%f\n", fdd);
  for(i=0; i<=gridpoints-1; i++){
    thetai=startingtheta + (3.14159-startingtheta)*i/(gridpoints);
    integral += sin(thetai)*spacing/pow(L*L-4.0*d*L*cos(thetai)+4.0*d*d,0.5);
  }
  integral=integral*pow(epsilon-1,2)/(8*PI*epsilon*(epsilon+1));
  fdd += integral;
  //printf("%f\n", fdd);
  integral=0;
  for(i=0; i<=gridpoints-1; i++){
    thetai = startingtheta + (3.14159-startingtheta)*i/(gridpoints);
    integral += sin(thetai)*spacing*ftheta[i];
  }
  integral=integral*(epsilon-1)/2;
  fdd -= integral;
  //printf("%f\n", fdd);
  energy = fdd*epsilon*2*PI/(epsilon-1);
  printf("%f %f %f %d\n", d, energy, delta, gridpoints);
}

void computeb(){
  double thetai, thetaj=0;
  double prefactor=0;
  double spacing=0;
  spacing=(PI-startingtheta)/gridpoints;
  prefactor = spacing*(epsilon-1)/(2*sqrt(2)*PI*(epsilon+1));
  for(i=0; i<=gridpoints-1; i++){
    for(j=0; j<=gridpoints-1; j++){
      thetai=startingtheta + (3.14159-startingtheta)*i/(gridpoints);
      thetaj=startingtheta + (3.14159-startingtheta)*j/(gridpoints);
      if(i==j){
	B[i][j]=1;
      }
      B[i][j] += prefactor*Integral(thetai, thetaj)*sin(thetai);
      //printf("%d %d %f\n", i, j, B[i][j]);
    }
  }

}

void computek(){
  double spacing=0;
  double theta=0;
  double thetap=0;
  double sum=0;
  spacing=(PI-startingtheta)/gridpoints;
  for(i=0; i<=gridpoints-1; i++){
    theta=startingtheta + (3.14159-startingtheta)*i/(gridpoints);
    K[i] = (epsilon-1)*(epsilon-1)/(4*PI*epsilon*(epsilon+1)*sqrt(L*L-4*d*L*cos(theta)+4*d*d));
    sum=0;
    for(j=0; j<=gridpoints-1; j++){
      thetap = startingtheta + j*spacing;
      sum += Integral(theta, thetap)*sin(thetap)*(1/L + (epsilon-1)*(L-2*d*cos(thetap))*pow(L*L-4*d*L*cos(thetap)+4*d*d, -1.5)*L/(epsilon+1));
    }
    sum = sum*spacing*pow((epsilon-1)/(4*PI),2)/(epsilon*sqrt(2));
    K[i] -= sum;
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

