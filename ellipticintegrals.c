#include <stdio.h>
#include <gsl/gsl_sf_ellint.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>


double result=0;


int main(){
  double input=0;
  input = sqrt(0.5);
  result = gsl_sf_ellint_Kcomp(input, 0.01);
  printf("%.12f\n", result);
  return 0;
}
