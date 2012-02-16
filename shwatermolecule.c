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
double FLMpos[20][20][3][3] = {0};
double FLMneg[20][20][3][3] = {0};
double complex ffunction[3][3]={0};
double radii[3] = {0};
double charges[3] = {0};
double thetas[3] = {0};
double phis[3] = {0};
double epsilon = 80;
double radius = 2.9;
double depth = 10;
double depths[10]={0};
double smallr = 1;
double smallrp = 1;
double thetacoord, phicoord=0;
double thetacoordp, phicoordp=0;
int gridsize=0;
int lmax=4;
double Dsmallr, Dbigr, Abigr=0;
double Dsmallrp, Asmallrp=0;  //Dsmallrp involves r-prime, and the small r is also the first argument of the image charge
char enfile[128];

void initialize(void);
void projectontoYLM();
void imagecharges(double, double, double, double);
void calculateFLM(const int, const int);
void orientationalaverage(void);
void generateorientation(void);

int main(int argc, char *argv[]){
  sscanf(argv[1], "%d", &lmax);
  sscanf(argv[2], "%s", &enfile);
  sscanf(argv[3], "%d", &gridsize);
  int p=0;
  int p1, p2=0;
  initialize();
  for(p=0; p<=64; p++){
    depth = 22 - 0.3*p;
    for(p1=0; p1<=2; p1++){
      for(p2=0; p2<=2; p2++){
	smallr = radii[p1];
	smallrp = radii[p2];
	projectontoYLM();
	calculateFLM(p1, p2);
      }
    }
    orientationalaverage();
  }
  return 0;
}

void initialize(){
  radii[0]=0;
  radii[1]=1;
  radii[2]=1;
  charges[0]=-.8476;
  charges[1]=.4328;
  charges[2]=.4328;
}

void projectontoYLM(){
  //at the moment this is just using the rectangle method (I'm pathetic, I know)
  int Ntheta=0;
  int Nphi=0;
  double theta, thetap, phi, phip=0;
  double plmplm=0;
  double complex exponentials=0;
  double complex ylmylm=0;
  double complex value=0;
  double h=0;
  int i,j,k,l,m,n=0;
  Ntheta = gridsize;
  Nphi = 2*gridsize;
  h = pow(PI, 4)*4/(Ntheta*Ntheta*Nphi*Nphi);
  for(l=0; l<=lmax; l++){
    for(m=0; m<=l; m++){
      YLMAbigrpos[l][m] = 0;
      YLMAbigrneg[l][m] = 0;
      YLMDbigrpos[l][m] = 0;
      YLMDbigrneg[l][m] = 0;
      YLMDsmallrpos[l][m] = 0;
      YLMDsmallrneg[l][m] = 0;
      YLMAsmallrppos[l][m] = 0;
      YLMAsmallrpneg[l][m] = 0;
      YLMDsmallrppos[l][m] = 0;
      YLMDsmallrpneg[l][m] = 0;
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
      //printf("Abigr %f %f\n", creal(YLMAbigrpos[l][m]), cimag(YLMAbigrpos[l][m]));
      //printf("Dbigr %f %f\n", creal(YLMDbigrpos[l][m]), cimag(YLMDbigrpos[l][m]));
      //printf("Dsmallr %f %f\n", creal(YLMDsmallrpos[l][m]), cimag(YLMDsmallrpos[l][m]));
      //printf("Asmallrp %f %f\n", creal(YLMAsmallrppos[l][m]), cimag(YLMAsmallrppos[l][m]));
      //printf("Dsmallrp %f %f\n\n", creal(YLMDsmallrppos[l][m]), cimag(YLMDsmallrppos[l][m]));
    }
  }

  //printf("%f %f\n", cimag(YLMAbigrpos[3][0]), h);
}


void calculateFLM(const int i1, const int i2){
  int i,j,k,l,m=0;
  double complex numerator=0;
  double complex denominator=0;
  double complex energy=0;
  double complex exponentials;
  double complex flmsurfpos=0;
  double complex flmsurfneg=0;
  double ldum=0;
  double lp21, eps = 0;
  eps=epsilon;
  for(l=0; l<=lmax; l++){
    lp21 = 2*l+1;
    ldum=l;
    for(m=0; m<=l; m++){
      flmsurfpos=0;
      flmsurfneg=0;
      FLMpos[l][m][i1][i2]=0;
      FLMneg[l][m][i1][i2]=0;

      flmsurfpos = ((ldum+1)/(pow(lp21,2)))*pow(smallr, l)/pow(radius, l+1);
      flmsurfpos += ((ldum+1)/(lp21))*pow(smallr/radius, l)*YLMAbigrpos[l][m];
      flmsurfpos -= YLMDsmallrpos[l][m]*radius/(lp21);
      flmsurfpos -= YLMDsmallrpos[l][m]*YLMAbigrpos[l][m]*pow(radius, 2);
      flmsurfpos = -pow(eps-1, 2)*flmsurfpos/eps;
      denominator = ((ldum+1)*eps + ldum)/(lp21) - (eps-1)*YLMDbigrpos[l][m]*radius;
      flmsurfpos = flmsurfpos/denominator;
      FLMpos[l][m][i1][i2] = -((ldum+1)/(pow(lp21, 2)))*pow(smallr*smallrp, l)/pow(radius, 2*l+1);
      FLMpos[l][m][i1][i2] += pow(smallrp/radius, l)*YLMDsmallrpos[l][m]*radius/(lp21);
      FLMpos[l][m][i1][i2] += pow(radius, 2)*YLMAsmallrppos[l][m]*YLMDsmallrpos[l][m];
      FLMpos[l][m][i1][i2] -= ((ldum+1)/(lp21))*YLMAsmallrppos[l][m]*pow(smallr/radius, l);
      FLMpos[l][m][i1][i2] = FLMpos[l][m][i1][i2]*pow(eps-1, 2)/eps;
      FLMpos[l][m][i1][i2] += (eps-1)*flmsurfpos*(pow(radius, 2)*YLMDsmallrppos[l][m] - pow(smallrp/radius, l)*(ldum+1)/(lp21));


      flmsurfneg = ((ldum+1)/(pow(lp21,2)))*pow(smallr, l)/pow(radius, l+1);
      flmsurfneg += ((ldum+1)/(lp21))*pow(smallr/radius, l)*YLMAbigrneg[l][m];
      flmsurfneg -= YLMDsmallrneg[l][m]*radius/(lp21);
      flmsurfneg -= YLMDsmallrneg[l][m]*YLMAbigrneg[l][m]*pow(radius, 2);
      flmsurfneg = -pow(eps-1, 2)*flmsurfneg/eps;
      denominator = ((ldum+1)*eps + ldum)/(lp21) - (eps-1)*YLMDbigrneg[l][m]*radius;
      flmsurfneg = flmsurfneg/denominator;
      FLMneg[l][m][i1][i2] = -((ldum+1)/(pow(lp21, 2)))*pow(smallr*smallrp, l)/pow(radius, 2*l+1);
      FLMneg[l][m][i1][i2] += pow(smallrp/radius, l)*YLMDsmallrneg[l][m]*radius/(lp21);
      FLMneg[l][m][i1][i2] += pow(radius, 2)*YLMAsmallrpneg[l][m]*YLMDsmallrneg[l][m];
      FLMneg[l][m][i1][i2] -= ((ldum+1)/(lp21))*YLMAsmallrpneg[l][m]*pow(smallr/radius, l);
      FLMneg[l][m][i1][i2] = FLMneg[l][m][i1][i2]*pow(eps-1, 2)/eps;
      FLMneg[l][m][i1][i2] += (eps-1)*flmsurfneg*(pow(radius, 2)*YLMDsmallrpneg[l][m] - pow(smallrp/radius, l)*(ldum+1)/(lp21));
      //printf("%d %d %f %f\n", l, m, creal(FLMpos[l][m]), cimag(FLMpos[l][m]));
      //printf("neg %d %d %f %f\n", l, m, creal(FLMneg[l][m]), cimag(FLMneg[l][m]));
    }
  }

}

void orientationalaverage(){
  int i,j,k,l,m=0;
  double complex ylmylm=0;
  double complex exponentials=0;
  double complex energy=0;
  double plmplm=0;
  FILE *enptr;
  enptr=fopen(enfile, "a");
  generateorientation();
  for(i=0; i<=2; i++){
    for(j=0; j<=2; j++){
      ffunction[i][j]=0;
      thetacoord = thetas[i];
      thetacoordp = thetas[j];
      phicoord = phis[i];
      phicoordp = phis[j];
      for(l=0; l<=lmax; l++){
	m=0;
	plmplm = gsl_sf_legendre_sphPlm(l, 0, cos(thetacoord))*gsl_sf_legendre_sphPlm(l, 0, cos(thetacoordp));
	exponentials = cexp(I*phicoord*m)*cexp(-I*phicoordp*m);
	ylmylm = exponentials*plmplm;
	ffunction[i][j] += FLMpos[l][0][i][j]*ylmylm;
	for(m=1; m<=l; m++){
	  plmplm = gsl_sf_legendre_sphPlm(l, m, cos(thetacoord))*gsl_sf_legendre_sphPlm(l, m, cos(thetacoordp));
	  exponentials = cexp(I*phicoord*m)*cexp(-I*phicoordp*m);
	  ylmylm = exponentials*plmplm;
	  ffunction[i][j] += FLMpos[l][m][i][j]*ylmylm;
	  plmplm = plmplm*pow(gsl_sf_fact(l-m), 2)/pow(gsl_sf_fact(l+m), 2);
	  exponentials = conj(exponentials);
	  ylmylm = exponentials*plmplm;
	  ffunction[i][j] += FLMneg[l][m][i][j]*ylmylm;
	}
      }
      energy = charges[i]*charges[i]*ffunction[i][j]*epsilon*4*PI/(epsilon-1);
    }
  }

  printf("energies .. %f %f\n", creal(energy), cimag(energy));
  fprintf(enptr, "%f %f\n", depth, creal(energy));
  fclose(enptr);
}


void generateorientation(void){
  thetas[0] = 0.5;
  thetas[1] = 0.0;
  thetas[2] = 1.902;
  phis[0] = 0.5;
  phis[1] = 0.0;
  phis[2] = 0.0;
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
  //printf("%f %f\n", Dsmallr, Dsmallrp);
}
