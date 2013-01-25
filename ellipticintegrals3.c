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
  double thetai=0;
  computeb();
  computek();
  mat A(gridpoints, gridpoints);
  //startingtheta=0;
  startingtheta=acos(d/L);
  spacing=PI/gridpoints;
  for(i=0; i<=gridpoints-1; i++){
    for(j=0; j<=gridpoints-1; j++){
      A(i,j)=B[i][j];
    }
  }
  vec k(gridpoints);
  for(i=0; i<=gridpoints-1; i++){
    thetai=spacing*i;
    if(thetai<startingtheta){
      k(i)=K[i]*(epsilon+1)/(2*epsilon);
    }
    else{
      k(i)=K[i]*2/(epsilon+1);
    }
  }
  //k.print();
  mat C = inv(A);
  vec f = C*k;
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
  spacing=PI/gridpoints;
  //printf("%f\n", fdd);
  integral=0;
  
  fdd=pow(epsilon-1,2)/(8*PI*d*(epsilon+1));
  for(i=0; i<=gridpoints-1; i++){
    thetai = spacing*i;
    if(thetai<startingtheta){
      integral += tan(thetai)*2*epsilon/(d*(epsilon+1));
    }
    else{
      integral +=sin(thetai)*(1/L + (epsilon-1)/((epsilon+1)*sqrt(L*L-4*d*L*cos(thetai)+4*d*d)));
    }
  }
  integral = integral*spacing*(epsilon-1)*(epsilon-1)/(8*PI*epsilon);
  fdd -= integral;

  integral=0;
  //printf("%f\n", fdd);
  for(i=0; i<=gridpoints-1; i++){
    thetai = spacing*i;
    if(thetai<startingtheta){
      integral -= tan(thetai)*ftheta[i]*2*epsilon/(cos(thetai)*(epsilon+1));
    }
    else{
      integral += sin(thetai)*(-1/L*L + (epsilon-1)*(2*d*cos(theta)-L)/((epsilon+1)*pow(L*L-4*d*L*cos(thetai+4*d*d,1.5)));
    }
  }
  integral = integral*spacing*(epsilon-1)/2;
  fdd += integral;
  //printf("%f\n", fdd);
  energy = fdd*epsilon*2*PI/(epsilon-1);
  printf("%f %f %f %d\n", d, energy, delta, gridpoints);

}


void computeb(){
  double thetai, thetaj=0;
  double prefactor=0;
  double spacing=0;

  spacing=PI/gridpoints;
  for(j=0; j<=gridpoints-1; j++){
    thetaj=3.14159*j/gridpoints;
    if(thetaj<startingtheta){
      prefactor = spacing*(epsilon-1)*(epsilon+1)/(8*PI*epsilon);
    }
    else{
      prefactor = spacing*(epsilon-1)/(2*PI*(epsilon+1));
    }
    for(i=0; i<=gridpoints-1; i++){
      thetai=3.14159*i/gridpoints;
      if(i==j){
	B[i][j]=1;
      }
      B[i][j] -= prefactor*Integral(thetai, thetaj);
      //printf("%d %d %f\n", i, j, B[i][j]);
    }
  }

}


void computek(){
  double spacing=0;
  double theta=0;
  double thetap=0;
  double sum=0;
  double fullintegral=0;
  spacing=PI/gridpoints;
  for(i=0; i<=gridpoints-1; i++){
    theta = i*spacing;
    if(thetai<startingtheta){
      K[i] = (pow(epsilon,2)-1)*(epsilon-1)*cos(theta)/(d*8.0*PI*epsilon*(epsilon+1));
      K[i] -= cos(theta)*pow(epsilon-1,2)/(d*8.0*PI*epsilon);
    }
    else{
      K[i] = (pow(epsilon,2)-1)*(epsilon-1)/(8.0*PI*epsilon*(epsilon+1)*sqrt(L*L-4*L*d*cos(theta)+4*d*d));
      K[i] -= pow(epsilon-1,2)/(L*8.0*PI*epsilon);
    }
    sum=0;
    for(j=0; j<=gridpoints-1; j++){
      thetap = j*spacing;
      fullintegral = Integral2(theta, thetap);
      sum += fullintegral;
    }
    sum = sum*spacing*pow((epsilon-1)/(4*PI),2)/(epsilon);
    K[i] += sum;
  }

}

double Integral2(double thetai, double thetaj){
  double A, B=0;
  double argument=0;
  double integral=0;
  double Jac=0;
  double factor1, factor2=0;
  double tani, tanj=0;
  double element=0;
  tani = tan(thetai);
  tanj = tan(thetaj);
  if(thetai<startingtheta){
    if(thetaj<startingtheta){//r'' and R on boundary
      Jac=tani*d*d/(cos(thetai)*cos(thetai));
      A = d*d*(tani*tani+tanj*tanj) + delta;
      B = 2.0*d*tani*tanj;
      argument = -2*B/(A-B);
      integral = 2*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A-B));
      argument = 2*B/(A+B);
      integral += 2*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A+B));
      integral = integral*Jac*(-2*epsilon*cos(thetaj))/(d*d*(epsilon+1));
    }
    else{//r'' on sphere R on boundary
      Jac=L*L*sin(thetai);
      A = d*d*tani*tani + d*d -2*d*L*cos(thetaj) + L*L;
      B = 2*d*L*tani*sin(thetaj);
      argument = -2*B/(A-B);
      integral = 2*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A-B));
      argument = 2*B/(A+B);
      integral += 2*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A+B));
      integral = integral*Jac*(-2*epsilon/(L*L*(epsilon+1)));      
    }
  }
  else{
    if(thetaj<startingtheta){//r'' on boundary R on sphere
      Jac=tani*d*d/(cos(thetai)*cos(thetai));
      A = d*d*tanj*tanj + d*d -2*d*L*cos(thetai) + L*L;
      B = 2*d*L*tanj*sin(thetai);
      argument = -2*B/(A-B);
      integral = 2*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A-B));
      argument = 2*B/(A+B);
      integral += 2*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A+B));
      integral = integral*Jac*(-2*epsilon*cos(thetaj))/(d*d*(epsilon+1));
    }
    else{// r'' and R on sphere
      Jac=L*L*sin(thetai);
      A = 2*L*L - 2*L*L*cos(thetai)*cos(thetaj) + delta;
      B = 2*L*L*sin(thetai)*sin(thetaj); 
      argument = -2*B/(A-B);
      integral = 2*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A-B));
      argument = 2*B/(A+B);
      integral += 2*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A+B));
      A = 2*L*L - 2*L*L*cos(thetai)*cos(thetaj) -4*L*d*(cos(thetai)+cos(thetaj)) + 4*d*d + delta;
      B = 2*L*L*sin(thetai)*sin(thetaj);
      argument=-2*B/(A-B);
      integral += 2*(epsilon-1)*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A-B)*(epsilon+1));
      argument = 2*B/(A+B);
      integral += 2*(epsilon-1)gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A+B)*(epsilon+1));      
      integral = integral*Jac*(-1/(L*L));
    }
  }

}

double Integral(double thetai, double thetaj){
  double A, B=0;
  double argument=0;
  double integral=0;
  double Jac=0;
  double factor1, factor2=0;
  double factor3, factor4=0;
  double tani, tanj=0;
  //i corresponds to the double prime coordinate
  tani = tan(thetai);
  tanj = tan(thetaj);
  if(thetai<startingtheta){
    Jac=tani*d*d/(cos(thetai)*cos(thetai));
    if(thetaj<startingtheta){//corresponds to r'' and R on interface
      A = d*d*(pow(tani,2)+pow(tanj,2))+delta;
      B = 2*d*tani*tanj;
      factor1 = 4*EllipticE(-2*B/(A-B))/(sqrt(A-B)*(A+B));
      factor1 = factor1*d*(pow(tani,2)+pow(tanj,2));
      factor2 = 4*A*(A-B)*EllipticE(-2*B/(A-B));
      factor2 -= 4*(A*A-B*B)*ellipticK(-2*B/(A-B));
      factor2 = factor2/(B*(A+B)*pow(A-B,1.5));
      factor2 = factor2*tan(thetai)*tan(thetaj);
      integral = factor2 - factor1;
      integral = integral*Jac*2*epsilon/(epsilon+1);
    }
    else{//corresponds to r'' on interface and R on sphere
      A = pow(d*tani,2) + d*d -2*d*L*cos(thetaj)+L*L + delta;
      B = 2*d*L*tani*sin(thetaj);
      factor1 = 4*EllipticE(-2*B/(A-B))/(sqrt(A-B)*(A+B));
      factor1 = factor1*(d*tani*tani + d - L*cos(thetaj));
      factor2 = 4*A*(A-B)*EllipticE(-2*B/(A-B));
      factor2 -= 4*(A*A-B*B)*ellipticK(-2*B/(A-B));
      factor2 = factor2*L*tani*sin(thetaj);
      integral = factor2 - factor1;
      integral = integral*Jac*2*epsilon/(epsilon+1);
    }
  }
  else{
    Jac=L*L*sin(thetai);
    if(thetaj<startingtheta){//corresponds to r'' on sphere and R on interface
      A = pow(d*tanj,2) + d*d - 2*d*L*cos(thetai) + L*L +delta;
      B = 2*d*L*tanj*sin(thetai);
      factor1 = 4*EllipticE(-2*B/(A-B))/(sqrt(A-B)*(A+B));
      factor1 = factor1*(cos(thetai)*d-L);
      factor2 = 4*A*(A-B)*EllipticE(-2*B/(A-B));
      factor2 -= 4*(A*A-B*B)*ellipticK(-2*B/(A-B)); 
      factor2 = factor2*d*tanj*sin(thetai);
      integral = factor1 + factor2;
      integral = integral*Jac*2*epsilon/(epsilon+1);
    }
    else{//corresponds to r'' and R on sphere
      A = 2*L*L-2*L*L*cos(thetai)*cos(thetaj)+delta;
      B = 2*L*L*sin(thetai)*sin(thetaj);
      factor1 = 4*EllipticE(-2*B/(A-B))/(sqrt(A-B)*(A+B));
      factor1 = factor1*(-2*L-2*L*cos(thetai)*cos(thetaj));
      factor2 = 4*A*(A-B)*EllipticE(-2*B/(A-B));
      factor2 -= 4*(A*A-B*B)*ellipticK(-2*B/(A-B));
      factor2 = factor2*2.0*L*sin(thetai)*sin(thetaj);
      A = 2*L*L + 2*L*L*cos(thetai)*cos(thetaj)-4.0*L*d*(cos(thetai)+cos(thetaj)) + 4*d*d + delta;
      //B stays the same
      factor3 = 4*EllipticE(-2*B/(A-B))/(sqrt(A-B)*(A+B));
      factor3 = factor3*(-2*L-2*L*cos(thetai)*cos(thetaj) - 4*d*(cos(thetai)+cos(thetaj)));
      factor4 = 4*A*(A-B)*EllipticE(-2*B/(A-B));
      factor4 -= 4*(A*A-B*B)*ellipticK(-2*B/(A-B));
      factor4 = factor4*2.0*L*sin(thetai)*sin(thetaj);
      integral = factor1+factor2+(epsilon-1)*(factor3 + factor4)/(epsilon+1);
      integral = integral*Jac;
      
    }
  }

  return integral;
}

