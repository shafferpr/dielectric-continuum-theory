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
int ij[15000][2]={0};
double fsurf[15000]={0};
double fsurfnew[15000]={0};
double rho[15000][7]={0}; // This is where I will store coordinates of the triangulation
double delta=0.01;
double epsilon=80;
double L, d=0;
double startingtheta=0;
double derivativedelta=0.03;
int gridpoints=10;
int gridsize1, gridsize2=0;
int Npoints=0;
int i,j,k=0;
char enfile[128];
void computef();
void triangulate();
void initialize();
void check();
double Integral2(double, double);
double computedelg1(const int);
double computedelg2(const int, const int);
double computeg(const int, const int);
double computeg1(const int);
void computefsnew();
void printsurf(const int);
void update();

int main(int argc, char *argv[]){
  sscanf(argv[1], "%lf", &L);
  sscanf(argv[2], "%lf", &d);
  sscanf(argv[3], "%d", &Npoints);
  sscanf(argv[4], "%lf", &delta);
  sscanf(argv[5], "%lf", &derivativedelta);
  sscanf(argv[6], "%s", &enfile);
  double spacing=0;
  double thetai=0;
  int iterate, height=0;
  triangulate();
  d=30;
  initialize();
  for(height=30; height>=2; height=height-1){
    d=height;
    for(iterate=0; iterate<=115; iterate++){
      printf("%d %f %f\n", iterate, fsurf[500], d);
      if(iterate%10==0){
	printsurf(iterate);
      }
      computefsnew();
    }
    computef();
    update();
  }
  return 0;
}

void update(){

  for(i=0; i<=Npoints-1; i++){
    fsurf[i] -= pow((epsilon-1)/epsilon,2)/(4*PI)*((epsilon-1)/(2*(d+1)*(epsilon+1)));
    fsurf[i] += pow((epsilon-1)/epsilon,2)/(4*PI)*((epsilon-1)/(2*d*(epsilon+1)));
  }

}

void printsurf(const int marker){
  FILE *fsurfptr;
  char fsfile[128];
  sprintf(fsfile, "fs35_%d.txt", marker);
  fsurfptr=fopen(fsfile, "w");
  for(i=0; i<=Npoints-1; i++){
    fprintf(fsurfptr, "%f %f\n", rho[i][2], fsurf[i]); 
  }
  fclose(fsurfptr);
}


void initialize(){
  double initialvalue=0;
  initialvalue = -pow((epsilon-1)/epsilon,2)/(4*PI)*(1/L + (epsilon-1)/(2*d*(epsilon+1)));

  for(i=0; i<=Npoints-1; i++){
    fsurf[i]=initialvalue;
  }

}

void computefsnew(){
  double x0, y0, z0=0;
  double xp, yp, zp=0;
  double sum1, sum2=0;

  for(i=0; i<=Npoints-1; i++){
    sum1=0;
    sum2=0;
    for(j=0; j<=Npoints-1; j++){
      if(i!=j){
	sum1 += computedelg1(j)*computeg(i,j);
	sum2 +=fsurf[j]*computedelg2(j,i);
      }
    }
    sum1 = sum1*rho[500][6];
    sum2 = sum2*rho[500][6];
    fsurfnew[i] = sum1*2.0*pow((epsilon-1)/(4*PI),2)/(epsilon*(epsilon+1));
    fsurfnew[i] += sum2*2.0*(epsilon-1)/(4*PI*(epsilon+1));
  }

  for(i=0; i<=Npoints-1; i++){
    fsurf[i]=fsurfnew[i];
  }

}

double computedelg1(const int m){
  double x1, y1, z1=0;
  double z0=0;
  double delg=0;
  double g1, g2=0;
  double invr, invrz=0;
  
  x1=rho[m][0];
  y1=rho[m][1];
  z1=rho[m][2];
  invr = sqrt(pow(x1,2) + pow(y1,2) + pow(z1,2)+delta);
  invr = 1/invr;
  invrz = sqrt(pow(x1,2) + pow(y1,2) + pow(z1-2*d,2));
  invrz = 1/invrz;
  g1 = invr + (epsilon-1)*invrz/(epsilon+1);
  //g1 = invr; //remove after check
  x1=rho[m][3];
  y1=rho[m][4];
  z1=rho[m][5];
  invr = sqrt(pow(x1,2) + pow(y1,2) + pow(z1,2)+delta);
  invr = 1/invr;
  invrz = sqrt(pow(x1,2) + pow(y1,2) + pow(z1-2*d,2));
  invrz = 1/invrz;
  g2 = invr + (epsilon-1)*invrz/(epsilon+1);
  //g2 = invr;//me too
  delg = (g1-g2)/derivativedelta;
  return delg;
}

double computedelg2(const int m, const int n){
  double x1, y1, z1=0;
  double x2, y2, z2=0;
  double delg=0;
  double g1, g2=0;
  double invr, invrz=0;
  
  x1 = rho[m][0];
  y1 = rho[m][1];
  z1 = rho[m][2];
  x2 = rho[n][0];
  y2 = rho[n][1];
  z2 = rho[n][2];

  invr = sqrt(pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2)+delta);
  invr = 1/invr;
  invrz = sqrt(pow(x1-x2,2) + pow(y1-y2,2) + pow(z1+z2-2*d,2));
  invrz = 1/invrz;
  g1 = invr + (epsilon-1)*invrz/(epsilon+1);
  //g1 = invr;
  x1 = rho[m][3];
  y1 = rho[m][4];
  z1 = rho[m][5];

  invr = sqrt(pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2)+delta);
  invr = 1/invr;
  invrz = sqrt(pow(x1-x2,2) + pow(y1-y2,2) + pow(z1+z2-2*d,2));
  invrz = 1/invrz;
  g2 = invr + (epsilon-1)*invrz/(epsilon+1);
  //g2 = invr;
  delg = (g1-g2)/derivativedelta;
  return delg;
}

double computeg(const int m, const int n){
  double x1, y1, z1=0;
  double x2, y2, z2=0;
  double g=0;
  double invr, invrz=0;

  x1 = rho[m][0];
  y1 = rho[m][1];
  z1 = rho[m][2];
  x2 = rho[n][0];
  y2 = rho[n][1];
  z2 = rho[n][2];

  invr = sqrt(pow(x1-x2,2) + pow(y1-y2,2) + pow(z1-z2,2)+delta);
  invr = 1/invr;
  invrz = sqrt(pow(x1-x2,2) + pow(y1-y2,2) + pow(z1+z2-2*d,2));
  invrz = 1/invrz;
  g = invr + (epsilon-1)*invrz/(epsilon+1);
  //g = invr;
 
  return g;
}

/*
void triangulate(){
  double phii, thetaj=0;
  double x1, y1, z1=0;
  double L1, L2, L3=0;
  double v1, v2, v3=0;
  double q1, q2, q3=0;
  double area, perimeter=0;
  double AREA=0;
  int place=0;
  int place2, place3=0;
  int j1, i1=0;
  int j2, i2=0;
  int j3, i3=0;
  FILE *xyzptr;
  gridsize2=gridsize1/(2*sin(PI/3));
  //printf("%d %d\n", gridsize1, gridsize2);
  for(i=0; i<=gridsize1-1; i++){
    for(j=0; j<=gridsize2-1; j++){
      if(j%2==0){
	phii = 2*PI*i/gridsize1;
      }
      else{
	phii = 2*PI*i/gridsize1 + PI/gridsize1;
      }
      thetaj=2*PI*sin(PI/3)*(j+1)/gridsize1;
      place = i*gridsize2 + j;
      x1 = L*sin(thetaj)*cos(phii);
      y1 = L*sin(thetaj)*sin(phii);
      z1 = L*cos(thetaj);
      rho[place][0]=x1;
      rho[place][1]=y1;
      rho[place][2]=z1;
      Npoints+=1;
      ij[place][0]=i;
      ij[place][1]=j;
    }    
  }

  for(i=0; i<=Npoints-1; i++){
    i1 = ij[i][0];
    j1 = ij[i][1];
    if(j1%2==0){
      j2 = j1;
      i2 = i1+1;
      j3 = j1+1;
      i3 = i1;
      if(i1==gridsize1-1){
	i2=0;
	j2=j1;
      }
      place2 = i2*gridsize2 + j2;
      place3 = i3*gridsize2 + j3;
      L1 = sqrt(pow(rho[i][0]-rho[place2][0],2)+pow(rho[i][1]-rho[place2][1],2)+pow(rho[i][2]-rho[place2][2],2));
      L2 = sqrt(pow(rho[i][0]-rho[place3][0],2)+pow(rho[i][1]-rho[place3][1],2)+pow(rho[i][2]-rho[place3][2],2));
      L3 = sqrt(pow(rho[place3][0]-rho[place2][0],2)+pow(rho[place3][1]-rho[place2][1],2)+pow(rho[place3][2]-rho[place2][2],2));
      perimeter = (L1 + L2 + L3)/2;
      area = sqrt(perimeter*(perimeter-L1)*(perimeter-L2)*(perimeter-L3));
      //printf("%f %d %d\n", area, i1, j1);
      rho[i][6]=2*area;
      AREA +=2*area;
      v1 = rho[place2][0]-rho[i][0];
      v2 = rho[place2][1]-rho[i][1];
      v3 = rho[place2][2]-rho[i][2];
      q1 = rho[place3][0]-rho[i][0];
      q2 = rho[place3][1]-rho[i][1];
      q3 = rho[place3][2]-rho[i][2];
      L1 = v2*q3-v3*q2;
      L2 = v3*q1-v1*q3;
      L3 = v1*q2-v2*q1;
      area=sqrt(L1*L1+L2*L2+L3*L3);
      L1=L1/area;
      L2=L2/area;
      L3=L3/area;
      rho[i][3]=rho[i][0]+L1*derivativedelta;
      rho[i][4]=rho[i][1]+L2*derivativedelta;
      rho[i][5]=rho[i][2]+L3*derivativedelta;
    }
    else{
      j2 = j1-1;
      i2 = i1+1;
      j3 = j1;
      i3 = i1+1;
      if(i1==gridsize1-1){
	i2=0;
	i3=0;
      }
      place2 = i2*gridsize2 + j2;
      place3 = i3*gridsize2 + j3;
      L1 = sqrt(pow(rho[i][0]-rho[place2][0],2)+pow(rho[i][1]-rho[place2][1],2)+pow(rho[i][2]-rho[place2][2],2));
      L2 = sqrt(pow(rho[i][0]-rho[place3][0],2)+pow(rho[i][1]-rho[place3][1],2)+pow(rho[i][2]-rho[place3][2],2));
      L3 = sqrt(pow(rho[place2][0]-rho[place3][0],2)+pow(rho[place2][1]-rho[place3][1],2)+pow(rho[place2][2]-rho[place3][2],2));
      perimeter = (L1 + L2 + L3)/2;
      area = sqrt(perimeter*(perimeter-L1)*(perimeter-L2)*(perimeter-L3));
      //printf("%f %d %d %f %f %f\n", area, i1, j1, L1, L2, L3);
      rho[i][6]=2*area;
      AREA +=2*area;
      v1 = rho[place2][0]-rho[i][0];
      v2 = rho[place2][1]-rho[i][1];
      v3 = rho[place2][2]-rho[i][2];
      q1 = rho[place3][0]-rho[i][0];
      q2 = rho[place3][1]-rho[i][1];
      q3 = rho[place3][2]-rho[i][2];
      L1 = v2*q3-v3*q2;
      L2 = v3*q1-v1*q3;
      L3 = v1*q2-v2*q1;
      area=sqrt(L1*L1+L2*L2+L3*L3);
      L1=L1/area;
      L2=L2/area;
      L3=L3/area;
      rho[i][3]=rho[i][0]+L1*derivativedelta;
      rho[i][4]=rho[i][1]+L2*derivativedelta;
      rho[i][5]=rho[i][2]+L3*derivativedelta;
    }
  }
  //printf("%f\n", AREA);
  xyzptr=fopen("triangulation.xyz", "w");
  fprintf(xyzptr, "%d\n\n", 2*Npoints);
  for(i=0; i<=Npoints-1; i++){
    fprintf(xyzptr, "C %f %f %f\n", rho[i][0], rho[i][1], rho[i][2]);
    fprintf(xyzptr, "O %f %f %f\n", rho[i][3], rho[i][4], rho[i][5]);
  }

  fclose(xyzptr);

  }*/

void triangulate(){
  double inc, off=0;
  double x1, y1, z1=0;
  double L1, L2, L3=0;
  double v1, v2, v3=0;
  double q1, q2, q3=0;
  double r1, phi=0;
  double area, perimeter=0;
  double AREA=0;
  int place=0;
  int place1, place2=0;
  int j1, i1=0;
  int j2, i2=0;
  int j3, i3=0;
  FILE *xyzptr;
  
  inc = PI*(3-sqrt(5));
  off = 2/(float)Npoints;
  for(i=0; i<=Npoints-1; i++){
    y1 = i*off -1 + off/2;
    r1 = sqrt(1-y1*y1);
    phi = i*inc;
    rho[i][0]=cos(phi)*r1*L;
    rho[i][1]=y1*L;
    rho[i][2]=sin(phi)*r1*L;
  }
  for(i=0; i<=Npoints-1; i++){
    L1=30;
    L2=33;
    place1=0;
    place2=0;
    for(j=0; j<=Npoints-1; j++){
      if(j!=i){
	v1 = sqrt(pow(rho[i][0]-rho[j][0],2)+pow(rho[i][1]-rho[j][1],2)+pow(rho[i][2]-rho[j][2],2));
	if(v1<L1){
	  L1=v1;
	  place1=j;
	}
	else if(v1<L2){
	  L2=v1;
	  place2=j;
	}
      }
    }

    v1 = rho[place1][0]-rho[i][0];
    v2 = rho[place1][1]-rho[i][1];
    v3 = rho[place1][2]-rho[i][2];
    q1 = rho[place2][0]-rho[i][0];
    q2 = rho[place2][1]-rho[i][1];
    q3 = rho[place2][2]-rho[i][2];    
    L1 = rho[i][0];
    L2 = rho[i][1];
    L3 = rho[i][2];
    area=sqrt(L1*L1+L2*L2+L3*L3);
    L1=L1/area;
    L2=L2/area;
    L3=L3/area;
    rho[i][3]=rho[i][0]-L1*derivativedelta;
    rho[i][4]=rho[i][1]-L2*derivativedelta;
    rho[i][5]=rho[i][2]-L3*derivativedelta;
    rho[i][6] = 4*PI*L*L/(float)Npoints;
  }

  //printf("%f\n", AREA);
  xyzptr=fopen("triangulation.xyz", "w");
  fprintf(xyzptr, "%d\n\n", 2*Npoints);
  for(i=0; i<=Npoints-1; i++){
    fprintf(xyzptr, "C %f %f %f\n", rho[i][0], rho[i][1], rho[i][2]);
    fprintf(xyzptr, "O %f %f %f\n", rho[i][3], rho[i][4], rho[i][5]);
  }

  fclose(xyzptr);

}


void computef(){
  double fdd=0;
  double sum1, sum2=0;
  double energy=0;
  double area=0;
  FILE *enptr;
  enptr=fopen(enfile, "a");
  for(i=0; i<=Npoints-1; i++){
    sum1 += computedelg1(i)*computeg1(i)*rho[i][6];
    sum2 += computedelg1(i)*fsurf[i]*rho[i][6];
    area += rho[i][6];
  }
  sum1 = sum1*pow((epsilon-1)/(4*PI),2)/epsilon;
  sum2 = sum2*(epsilon-1)/(4*PI);
  fdd = sum1+sum2;
  energy = fdd*epsilon*2*PI/(epsilon-1) + (epsilon-1)/(4*d*(epsilon+1)); //important change here, this 4 used to be a 2
  //energy = fdd*epsilon*2*PI/(epsilon-1);
  fprintf(enptr, "%f %f %f %f %f %f\n", d, energy, fdd, sum1, sum2, area);
  fclose(enptr);
}


double computeg1(const int m){
  double x1, y1, z1=0;
  double x2, y2, z2=0;
  double g=0;
  double invr, invrz=0;

  x1=rho[m][0];
  y1=rho[m][1];
  z1=rho[m][2];
  invr = sqrt(pow(x1,2) + pow(y1,2) + pow(z1,2)+delta);
  invr = 1/invr;
  invrz = sqrt(pow(x1,2) + pow(y1,2) + pow(z1-2*d,2));
  invrz = 1/invrz;
  g = invr + (epsilon-1)*invrz/(epsilon+1);
  //g = invr;

  return g;
}


/*
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
o    }
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
 
*/
