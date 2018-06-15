#include <stdio.h>
#include <math.h>
#include "nr.h"
#include "ssfunctions2.h"
#include "nrutil.h"
#include "rates2.h"
#include "vecfunc2.h"
#define TOL 1.0e-20
#define SOLARMASS 1.998e33
#define PI 3.14159265359
#define BOLTZMANN 1.38e-16
#define ELECTRONMASS 9.109e-28
#define PROTONMASS 1.67262e-24
#define CLIGHT2 8.9875518e20
#define N 3

double m,mdot,RS,alphapar,betapar,gammapar,rmin,rmax,step;
int mgrid;
double r,radius,eTemp,eDensity,iDensity,magField;

int main()
{
    double *x;
    int check=0;
  
    x=dvector(1,N);
    constants2();
    r=rmin;
    FILE *ftemps;
    ftemps=fopen("temperatures.txt","w");
    x[1]=1.0e12;
    x[2]=1.0e9;
    x[3]=0.9;
    for (int i=1;i<=mgrid;i++) {
        r=r*step;
        radius=r*RS;
        newt(x,N,&check,vecfunc2);
        double f=x[3];
        double qie=qiefunc2(x[1],x[2],x[3]);
        double qemi=qem2(x[2],x[3]);
        double qplus=qp2(x[3])*(1.0-f);
        //double realtempi=PROTONMASS*CLIGHT2*x[1]/BOLTZMANN;
        //double realtempe=ELECTRONMASS*CLIGHT2*x[2]/BOLTZMANN;
        fprintf(ftemps,"%f %f %f %f\n",log10(r),log10(x[1]),log10(x[2]),x[3]);
    }
    fclose(ftemps);
    return 0;
}
