#include <stdio.h>
#include <gsl/gsl_sf_legendre.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
double YLMcomponents[20][40]={0};
void projectontoYLM();


int main(void){
  double x = 5.0;
  double y = gsl_sf_legendre_sphPlm (1, 1, .5);
  printf("J0(%f) = %.18f\n", x, y);
  return 0;
}


void projectontoYLM(){

}
