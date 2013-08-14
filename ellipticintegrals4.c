#include <stdio.h>
#include <gsl/gsl_sf_ellint.h>
#include <math.h>
#include <stdlib.h>
#include <armadillo>
#include <complex.h>
using namespace arma;
#define PI 3.14159265
//gcc -O2 -lgsl -lgslcblas shwatermolecule.c -o shwatermolecule
//g++ -O2 ellipticintegrals.c -o ellipticint -lgsl -lgslcblas -lblas -llapack
//this file computes solvation energy for d<L(the sphere penetrating the boundary), i have also switch back here to the original ansatz for chi_in^(-1), the one in which the expression for solvation energy involves an image charge
double result=0;
double B[800][800]={0};
double K[800]={0};
double ftheta[800]={0};
double delta=0.01;
double epsilon=80;
double L, d=0;
double startingtheta=0;
int gridpoints=10;
int i,j,k=0;
void computeb();
void computek();
void computef();
void check();
double Integral(double, double);
double Integral2(double, double);
int main(int argc, char *argv[]){
  sscanf(argv[1], "%lf", &L);
  sscanf(argv[2], "%lf", &d);
  sscanf(argv[3], "%d", &gridpoints);
  sscanf(argv[4], "%lf", &delta);
  double spacing=0;
  double thetai=0;
  startingtheta=0;
  //startingtheta=acos(d/L);
  computeb();
  computek();
  mat A(gridpoints, gridpoints);
  //printf("%lf %lf\n", d, startingtheta);
  spacing=PI/gridpoints;
  for(i=0; i<=gridpoints-1; i++){
    for(j=0; j<=gridpoints-1; j++){
      A(i,j)=B[i][j];
    }
    //printf("%lf %lf\n", i*spacing, B[i][0]);
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
    //printf("%lf %lf\n", thetai, K[i]);
  }

  mat C = inv(A);
  //C.print();
  mat D = A*C;
  //D.print();
  vec f = C*k;
  //vec k1 = A*f;
  /*for(i=0; i<=gridpoints-1; i++){
    A(i,i) -=1 ;
  }
  k1 = -A*f;
  k1 = k1 +k;*/
  //k1.print();
  //f.print();
  for(i=0; i<=gridpoints-1; i++){
    ftheta[i]=f(i);
    //printf("%lf %lf\n", spacing*i, ftheta[i]);
  }
  //check();
  computef();
  return 0;
}

void check(){
  double sum, spacing=0;
  spacing = PI/gridpoints;
  for(i=0; i<=gridpoints-1; i++){
    sum=0;
    for(j=0; j<=gridpoints-1; j++){
      if(i==j)
	B[i][j]-=1;
      sum -= B[i][j]*ftheta[j]; 
      //printf("%f\n", B[i][j]);
    }
    sum += K[i]*2/(epsilon+1);
    //printf("%f\n", K[i]);
    //printf("a %f %f\n", ftheta[i], sum);
  }

}

void computef(){
  double integral=0;
  double fdd=0;
  double spacing=0;
  double imageR=0;
  double thetai=0;
  double energy=0;
  double product=0;
  double firstterm=0;
  spacing=PI/gridpoints;
  
  integral=0;
  
  for(i=0; i<=gridpoints-1; i++){
    thetai = spacing*i;
    if(thetai<startingtheta){
      integral -= tan(thetai)*8*PI*epsilon/(d*(epsilon+1)*(epsilon+1));
    }
    else{
      imageR = L*L - 4*d*L*cos(thetai)+4*d*d;
      product = sin(thetai)*L*L*2*PI*(1/L + (epsilon-1)/((epsilon+1)*sqrt(imageR)));
      product = product*(-1/(L*L) + (epsilon-1)*(-L+2*d*cos(thetai))/((epsilon+1)*pow(imageR,1.5)));
      integral += product;
    }
  }
  integral = integral*spacing*pow((epsilon-1)/(4*PI),2)/epsilon;
  firstterm = integral;
  fdd += integral;

  integral=0;
  //printf("%f\n", fdd);
  for(i=0; i<=gridpoints-1; i++){
    thetai = spacing*i;
    if(thetai<startingtheta){
      integral -= tan(thetai)*ftheta[i]*2/(cos(thetai)*(epsilon+1));
    }
    else{
      imageR = L*L - 4*d*L*cos(thetai)+4*d*d;
      integral -= sin(thetai)*ftheta[i];
      //printf("A %f %f %f\n", thetai, ftheta[i], integral);
      integral += sin(thetai)*L*L*(epsilon-1)*(2*d*cos(thetai-L))/((epsilon+1)*pow(imageR,1.5));
      //printf("B %f %f %f\n", thetai, ftheta[i], integral);
    }

  }
  integral = integral*spacing*(epsilon-1)/2;
  fdd += integral;
  //printf("%f\n", fdd);
  energy = fdd*epsilon*2*PI/(epsilon-1) + (epsilon-1)/(2*d*(epsilon+1));
  printf("%f %f %f\n", d, energy, firstterm*epsilon*2*PI/(epsilon-1));

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
      //printf("%f %f %f\n", i,j, thetai, thetaj, startingtheta);
      if(i!=j)
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
    sum=0;
    for(j=0; j<=gridpoints-1; j++){
      thetap = j*spacing;
      fullintegral = Integral2(theta, thetap);
      sum += fullintegral;
    }
    sum = sum*spacing*pow((epsilon-1)/(4*PI),2)/(epsilon);
    K[i] = sum;
    //printf("%f\n", K[i]);
  }
  //printf("a %lf %lf\n", d, fullintegral);
}

double Integral2(double thetaj, double thetai){
  double A, B=0;
  double argument=0;
  double argument1=0;
  double integral=0;
  double Jac=0;
  double factor1, factor2=0;
  double tani, tanj=0;
  double element=0;
  tani = tan(thetai);
  tanj = tan(thetaj);
  if(thetai<startingtheta){
    if(thetaj<startingtheta){//r'' and R on boundary
      Jac=sin(thetai)/(cos(thetai)*cos(thetai));
      A = d*d*(tani*tani+tanj*tanj) + delta;
      B = 2.0*d*d*tani*tanj;
      argument = -2*B/(A-B);
      //printf("%f\n", argument);
      integral = 2*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A-B));
      //printf("%f\n", integral);
      argument = 2*B/(A+B);
      integral += 2*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A+B));
      //printf("%f\n", integral);
      integral = integral*Jac*(-4*epsilon)/(pow(epsilon+1,2));
    }
    else{//r'' on boundary R on sphere
      Jac=sin(thetai)/(cos(thetai)*cos(thetai));
      A = d*d*tani*tani + d*d -2*d*L*cos(thetaj) + L*L + delta;
      B = 2*d*L*tani*sin(thetaj);
      argument = -2*B/(A-B);
      integral = 2*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A-B));
      argument = 2*B/(A+B);
      integral += 2*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A+B));
      integral = integral*Jac*(-4*epsilon)/(pow(epsilon+1,2));
    }
  }
  else{
    if(thetaj<startingtheta){//r'' on sphere R on boundary
      Jac=L*L*sin(thetai)*(epsilon-1)*(-L+2*d*cos(thetai))/((epsilon+1)*pow(L*L-4*d*L*cos(thetai)+4*d*d,1.5));
      Jac -= sin(thetai);
      A = d*d*tanj*tanj + d*d -2*d*L*cos(thetai) + L*L + delta;
      B = 2*d*L*tanj*sin(thetai);
      argument = -2*B/(A-B);
      integral = 2*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A-B));
      argument = 2*B/(A+B);
      integral += 2*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A+B));
      integral = integral*Jac*(2*epsilon)/(epsilon+1);
    }
    else{// r'' on sphere and R sphere
      Jac=L*L*sin(thetai)*(epsilon-1)*(-L+2*d*cos(thetai))/((epsilon+1)*pow(L*L-4*d*L*cos(thetai)+4*d*d,1.5));
      Jac -= sin(thetai);
      A = 2*L*L - 2*L*L*cos(thetai)*cos(thetaj) + delta;
      B = 2*L*L*sin(thetai)*sin(thetaj); 
      argument = -2*B/(A-B);
      integral = 2*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A-B));
      argument = 2*B/(A+B);
      integral += 2*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A+B));
      //printf("1 %f\n", integral);
      A = 2*L*L + 2*L*L*cos(thetai)*cos(thetaj) -4*L*d*(cos(thetai)+cos(thetaj)) + 4*d*d + delta;
      B = 2*L*L*sin(thetai)*sin(thetaj);
      argument=-2*B/(A-B);
      integral += 2*(epsilon-1)*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A-B)*(epsilon+1));
      argument1 = 2*(epsilon-1)*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A-B)*(epsilon+1));
      argument = 2*B/(A+B);
      integral += 2*(epsilon-1)*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A+B)*(epsilon+1));
      //printf("2 %f\n", integral);
      argument1 += 2*(epsilon-1)*gsl_sf_ellint_RF(0, 1-argument, 1, 0.1)/(sqrt(A+B)*(epsilon+1));
      //printf("a1 %f\n", argument1);
      integral = integral*Jac;
    }
  }
  //if(thetai==0 && thetaj==0)
    //printf("%lf %lf\n", d, integral);
  return integral;
}

double Integral(double thetai, double thetaj){
  double A, B=0;
  double argument=0;
  double integral=0;
  double Jac=0;
  double factor1, factor2=0;
  double factor3, factor4=0;
  double tani, tanj=0;
  double EllipticE, EllipticK=0;
  //i corresponds to the double prime coordinate
  tani = tan(thetai);
  tanj = tan(thetaj);
  if(thetai<startingtheta){
    Jac=tani*d*d/(cos(thetai)*cos(thetai));
    if(thetaj<startingtheta){//corresponds to r'' and R on interface
      A = d*d*(pow(tani,2)+pow(tanj,2))+delta;
      B = 2*d*d*tani*tanj;
      argument = -2*B/(A-B);
      //printf("%f %f %f\n", argument, A, B);
      EllipticE = gsl_sf_ellint_RD(0, 1-argument, 1, 0.01)-argument*gsl_sf_ellint_RD(0,1-argument, 1, 0.01)/3;
      factor1 = 4*EllipticE/(sqrt(A-B)*(A+B));
      factor1 = factor1*d*(pow(tani,2)+pow(tanj,2));
      factor2 = 4*A*(A-B)*EllipticE;
      EllipticK = gsl_sf_ellint_RF(0, 1-argument, 1, 0.01);
      factor2 -= 4*(A*A-B*B)*EllipticK;
      //printf("%f\n", factor2);
      factor2 = factor2/(B*(A+B)*pow(A-B,1.5));
      //printf("%f\n", factor2);
      factor2 = factor2*tan(thetai)*tan(thetaj)*2*d;
      if(tan(thetai)==0 || tan(thetaj)==0)
	factor2=0;
      //printf("%f\n", factor2);
      integral = factor2 - factor1;
      integral = integral*Jac*2/(epsilon+1);//agree with everything up to here 07/28/13
      //printf("both less %f %f %f %f %f\n", integral, factor2, factor1, tan(thetai), tan(thetaj));
    }
    else{//corresponds to r'' on interface and R on sphere
      A = pow(d*tani,2) + d*d -2*d*L*cos(thetaj)+L*L + delta;
      B = 2*d*L*tani*sin(thetaj);
      argument = -2*B/(A-B);
      EllipticE = gsl_sf_ellint_RD(0, 1-argument, 1, 0.01)-argument*gsl_sf_ellint_RD(0,1-argument, 1, 0.01)/3;
      factor1 = 4*EllipticE/(sqrt(A-B)*(A+B));
      factor1 = factor1*(d*tani*tani + d - L*cos(thetaj));
      factor2 = 4*A*(A-B)*EllipticE;
      EllipticK = gsl_sf_ellint_RF(0, 1-argument, 1, 0.01);
      factor2 -= 4*(A*A-B*B)*EllipticK;
      factor2 = factor2*L*tani*sin(thetaj);
      integral = factor2 - factor1;
      integral = integral*Jac*2/(epsilon+1);//agree with everything up to here 07/28/13
      //printf("one less %f\n", integral);
    }
  }
  else{
    Jac=L*L*sin(thetai);
    if(thetaj<startingtheta){//corresponds to r'' on sphere and R on interface
      A = pow(d*tanj,2) + d*d - 2*d*L*cos(thetai) + L*L +delta;
      B = 2*d*L*tanj*sin(thetai);
      argument = -2*B/(A-B);
      EllipticE = gsl_sf_ellint_RD(0, 1-argument, 1, 0.01)-argument*gsl_sf_ellint_RD(0,1-argument, 1, 0.01)/3;
      factor1 = 4*EllipticE/(sqrt(A-B)*(A+B));
      factor1 = factor1*(cos(thetai)*d-L);
      factor2 = 4*A*(A-B)*EllipticE;
      EllipticK = gsl_sf_ellint_RF(0, 1-argument, 1, 0.01);
      factor2 -= 4*(A*A-B*B)*EllipticK; 
      B+=delta/2;
      factor2 = factor2*d*tanj*sin(thetai)/(pow(A-B, 1.5)*B*(A+B));
      integral = factor1 + factor2;
      integral = integral*Jac*2*epsilon/(epsilon+1);//still agree 07/28/13
      //printf("one less %f\n", integral);
    }
    else{//corresponds to r'' and R on sphere
      A = 2-2*cos(thetai)*cos(thetaj)+delta;
      B = 2*sin(thetai)*sin(thetaj);
      argument = -2*B/(A-B);
      EllipticE = gsl_sf_ellint_RD(0, 1-argument, 1, 0.01)-argument*gsl_sf_ellint_RD(0,1-argument, 1, 0.01)/3;
      factor1 = 4*EllipticE/(sqrt(A-B)*(A+B));
      factor1 = factor1*(-1/(L*L));
      A = 2*L*L + 2*L*L*cos(thetai)*cos(thetaj)-4.0*L*d*(cos(thetai)+cos(thetaj)) + 4*d*d + delta;
      B = 2*L*L*sin(thetai)*sin(thetaj);
      argument = -2*B/(A-B);
      EllipticE = gsl_sf_ellint_RD(0, 1-argument, 1, 0.01)-argument*gsl_sf_ellint_RD(0,1-argument, 1, 0.01)/3;
      factor3 = 4*EllipticE/(sqrt(A-B)*(A+B));
      factor3 = factor3*(-2*L-2*L*cos(thetai)*cos(thetaj) + 2*d*(cos(thetai)+cos(thetaj)));
      factor4 = 4*A*(A-B)*EllipticE;
      EllipticK = gsl_sf_ellint_RF(0, 1-argument, 1, 0.01);
      factor4 -= 4*(A*A-B*B)*EllipticK;
      //printf("A %lf\n", factor4);
      //if(Jac==0){
	//factor4 = 0;
      //}
      //else{
      B+=delta/2;
      factor4 = factor4*2.0*L*sin(thetai)*sin(thetaj)/(pow(A-B, 1.5)*B*(A+B));
	//}
      integral = factor1+(epsilon-1)*(factor3 + factor4)/(epsilon+1);
      integral = integral*Jac;//still agree 07/28/13

    }
  }

  return integral;
}

