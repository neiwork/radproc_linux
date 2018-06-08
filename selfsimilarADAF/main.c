#include <stdio.h>
#include <math.h>
#include "nr.h"
#include "ssfunctions.h"
#include "nrutil.h"
#include "rates.h"
#include "vecfunc.h"
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
    constants();
    r=rmin;
    FILE *ftemps;
    ftemps=fopen("temperatures.txt","w");
    x[1]=0.1;
    x[2]=0.1;
    x[3]=0.99;
    for (int i=1;i<=mgrid;i++) {
        r=r*step;
        radius=r*RS;
        newt(x,N,&check,vecfunc);
        double qie=qiefunc(x[1],x[2],x[3]);
        double qemi=qem(x[2],x[3]);
        double qplus=qp(x[3])*(1.0-x[3]);
        double realtempi=PROTONMASS*CLIGHT2*x[1]/BOLTZMANN;
        double realtempe=ELECTRONMASS*CLIGHT2*x[2]/BOLTZMANN;
        fprintf(ftemps,"%f %f %f %f\n",log10(r),log10(realtempi),log10(realtempe),log10(mdot));
    }
    fclose(ftemps);
    return 0;
}
