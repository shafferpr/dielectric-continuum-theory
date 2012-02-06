#include <stdio.h>
#include <gsl/gsl_sf_bessel.h>

int main(void){
  double x = 5.0;
  double y = gsl_sf_bessel_J0 (x);
  printf("J0(%f) = %.18f\n", x, y);
  return 0;
}
