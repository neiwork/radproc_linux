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

double m,RS,alphapar,betapar,gammapar,rmin,rmax,step;
//double mdot;
double f;
int mgrid;
double r,radius,eTemp,eDensity,iDensity,magField;

int main()
{
    double *x;
    int check=0;
  
    x=dvector(1,N);
    constants();
    r=rmin;
    FILE *ftemps, *fmdots, *fqs,*fne;
    ftemps=fopen("temperatures.txt","w");
	fne=fopen("ne.txt","w");
	fmdots=fopen("mdots.txt","w");
	fqs=fopen("qrates.txt","w");
	
    double inittempi=1.0e12;
	double inittempe=1.0e9;
	x[1]=(BOLTZMANN*inittempi)/(PROTONMASS*CLIGHT2);
    x[2]=(BOLTZMANN*inittempe)/(ELECTRONMASS*CLIGHT2);
    x[3]=1.0e-3;
    for (int i=1;i<=mgrid;i++) {
        r=r*step;
        radius=r*RS;
		double rad=r;
		broydn(x,N,&check,vecfunc);
		double mdot=x[3];
		//double f=x[3];
        double qie=qiefunc(x[1],x[2],mdot);
		//double qie=qiefunc2(x[1],x[2],f);
		double qSy,nucrit,qCbr,qCsy,qBr;
		qsync(&qSy,&nucrit,x[2]);
		qBr=qbremss(x[2]);
		qC(qBr,qSy,nucrit,x[2],&qCbr,&qCsy,mdot);
        double qemi=qem(x[2],mdot);
		double qplus=qp(mdot);
		//qBr=qbremss2(x[2]);
		//qsync2(&qSy,&nucrit,x[2]);
		//qC2(qBr,qSy,nucrit,x[2],&qCbr,&qCsy,f);
		//double qemi=qem2(x[2],f);
        //double qplus=qp2(f);
        double realtempi=PROTONMASS*CLIGHT2*x[1]/BOLTZMANN;
        double realtempe=ELECTRONMASS*CLIGHT2*x[2]/BOLTZMANN;
        fprintf(ftemps,"%f %f %f %f %f %f\n",log10(r),log10(realtempi),log10(realtempe),
			log10(ne(f)),log10(ni(f)),log10(magf(f)));
		fprintf(fne,"%f\n",ne(f));
		fprintf(fmdots,"%f %f %f\n",log10(r),log10(mdot),log10(mdotcrit()));
		fprintf(fqs,"%f %f %f %f\n",log10(r),log10(qplus),log10(qie),log10(qemi));
    }
    fclose(ftemps);
	fclose(fmdots);
	fclose(fqs);
    return 0;
}