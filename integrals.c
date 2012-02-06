#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <armadillo>
using namespace arma;

double B[5000][5000][3]={0};
double Binverse[5000][5000][3]={0};
double xyspacevalues[5000][3]={0};
double rx, ry, rz=0; //These are the coordinates of little r (the point on the interior of the sphere)
long int ngridpoints=0;

void calculateB();

int main(){
  calculateB();
  mat A(ngridpoints,ngridpoints);
  for(i=0; i<=ngridpoints; i++){
    for(j=0; j<=ngridpoints; j++){
      A(i,j) = B[i][j];
    }
  }

  mat BI(ngridpoints, ngridpoints);
  BI = inv(A, slow = false);  


}


calculateB(){


}
