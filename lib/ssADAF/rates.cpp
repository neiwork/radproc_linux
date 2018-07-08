#include "ssfunctions.h"
#include <fmath/constants.h>
#include <math.h>
extern "C" {
#include <nrMath/bessel.h>
#include <nrMath/nr.h>
}
#define TOL 1.0e-3

/*#define pi 3.14159265359
#define boltzmann 1.38e-16
#define PLANCK 6.6260755e-27
#define electronMass 9.109e-28
#define PROTONMASS 1.67262e-24
#define cLight2 8.9875518e20
#define stefanboltzmann 5.6704e-5*/


//Eq 3.4-3.8
double qbremss(double tempe)
{
    extern double eDensity,iDensity;
    double fei,fee,qei,qee;
  
    if (tempe <= 1.0) {
        fei=1.0159*sqrt(tempe)*(1.0+1.781*pow(tempe,1.34));
        fee=2.56e-22*pow(tempe,1.5)*(1.0+1.1*tempe+tempe*tempe-1.25*pow(tempe,2.5));
    } else {
        fei=1.4324*tempe*(log(1.123*tempe+0.48)+1.5);
        fee=3.4e-22*tempe*(log(1.1232*tempe)+1.28);
    }

    qei=1.48e-22*eDensity*iDensity*fei;  //Eq 3.5 ne² ? //deberia ser eDensity*eDensity
    qee=eDensity*eDensity*fee;

    return qei+qee;
}

//Eq. 3.14 
double auxf1(double x)
{
    extern double radius,eDensity,magField,eTemp;
    return exp(1.8899*pow(x,1.0/3.0))-2.49e-10*4.0*pi*eDensity*radius/magField/
        (eTemp*eTemp*eTemp*bessk(2,1.0/eTemp)*pow(x,5.0/3.0))*(0.5316+0.4*pow(x,0.25)+sqrt(x));
}																

void qsync(double *qSy, double *nucrit, double tempe)
{
    //int zbrac(double (*func)(double),double *x1, double *x2);
    //double zbrent(double (*func)(double), double x1, double x2, double tol);
    extern double radius,magField;
  
    double x1=1.0e-5;
    double x2=1.0e3;
    zbrac(auxf1,&x1,&x2);
    double xcrit=zbrent(auxf1,x1,x2,TOL);
    double nu0=2.8e6*magField;                                                        // [Hz]
    double auxnucrit=1.5*nu0*tempe*tempe*xcrit;  //Eq. 3.15
    *nucrit=auxnucrit;
    //printf("nu0 = %f, nuc = %f, r = %f\n",nu0,nucrit,r);
    double auxqSy=2.0*pi*electronMass/3.0 *tempe*(*nucrit)*(*nucrit)*(*nucrit)/radius;  //Eq 3.18
    *qSy=auxqSy;
}

//Eq. 3.23; 3.24
void qC(double qbr, double qsy, double nuc, double tempe, double *qCbr, double *qCsy, double mdot) 
{
    double prob,a,eta1,eta2,eta3;
    double xc=planck*nuc/(electronMass*cLight2);
    prob=1.0-exp(-taues(mdot));
    a=1.0+4.0*tempe*(1.0+4.0*tempe);
    eta1=prob*(a-1.0)/(1.0-prob*a);
    eta3=-1.0-log(prob)/log(a);
    eta2=pow(3.0,-eta3)*eta1;
    double auxqCbr=qbr*tempe*eta1*((1.0-xc/tempe)-3.0/(eta3+1.0) * (pow(1.0/3.0,eta3+1.0)-
						 pow(xc/(3.0*tempe),eta3+1.0))); //el xc lo pone straub y dice que corrije el typo de narayan & Yi 1995
    *qCbr=auxqCbr;
    double auxqCsy=qsy*(eta1-eta2*pow(xc/tempe,eta3));
    *qCsy=auxqCsy;
}

//Eq 3.1 o 3.3
double qiefunc(double tempi, double tempe, double mdot)  
{
    extern double r;
    double rad=r;
	double sumtemps=tempi+tempe;
    double xe=1.0/tempe;
    double xi=1.0/tempi;
    double xm=xe+xi;
    double besselk1,besselk0,besselk2e,besselk2i,result;
    if (xm > 400.0) {
        if (xi > 200.0) {
			if (xe > 200.0) {
				result=3.326668e-22*ne(mdot)*ni(mdot)*(1836.15*tempi-tempe)*sqrt(2.0*xm/(pi*xe*xi))*
					((2.0*sumtemps*sumtemps+1.0)/sumtemps+2.0);
			} else {
				besselk2e=bessk(2,xe);
				//result=3.328e-22*ne(f)*ni(f)*(1836.15*tempi-tempe)*sqrt(aux3/aux2)
				//    * ((2.0*aux1*aux1+1.0)/aux1 + 2.0) * exp(-(aux2+aux3))/besselk2e;
				result=3.326668e-22*ne(mdot)*ni(mdot)*(1836.15*tempi-tempe)*sqrt(xi/xm)*
					((2.0*sumtemps*sumtemps+1.0)/sumtemps+2.0)*exp(-xe)/besselk2e;
			}
        } else {
            besselk2i=bessk(2,xi);
            //result=3.328e-22*ne(f)*ni(f)*(1836.15*tempi-tempe)*sqrt(aux4/aux2)
            //       * ((2.0*aux1*aux1+1.0)/aux1 + 2.0) * exp(-(aux2+aux4))/besselk2i;
            result=3.326668e-22*ne(mdot)*ni(mdot)*(1836.15*tempi-tempe)*sqrt(xe/xm)*
                   ((2.0*sumtemps*sumtemps+1.0)/sumtemps+2.0)*exp(-xi)/besselk2i;
        }
    } else {
        besselk1=bessk1(xm);
        besselk0=bessk0(xm);
        besselk2i=bessk(2,xi);
        besselk2e=bessk(2,xe);
        //result=3.328e-22*ne(f)*ni(f)*(1836.15*tempi-tempe)*
        //    ((2.0*aux1*aux1+1.0)/aux1 * besselk1/besselk2i + 2.0*besselk0/besselk2i) / besselk2e;
        result=3.326668e-22*ne(mdot)*ni(mdot)*(1836.15*tempi-tempe)*
            ((2.0*sumtemps*sumtemps+1.0)/sumtemps * besselk1/besselk2i + 2.0*besselk0/besselk2i)/
			besselk2e;
    }
    
    return result;
}

//Eq. 3.33
double qem(double tempe, double mdot) 
{
    extern double eDensity,iDensity,magField,eTemp;
    double qBr,qSy,qCbr,qCsy,nucrit;
    eDensity=ne(mdot);
    iDensity=ni(mdot);
    magField=magf(mdot);
    eTemp=tempe;
    qsync(&qSy,&nucrit,tempe);
    qBr=qbremss(tempe);
    qC(qBr,qSy,nucrit,tempe,&qCbr,&qCsy,mdot);
    double qemi=qBr+qSy+qCbr+qCsy;
    double aux=4.0*stefanBoltzmann*pow(electronMass*cLight2*tempe/boltzmann,4)/height();
    double tauabs=qemi/aux;
    double tau=tauabs+taues(mdot);
    double result=aux*pow(1.5*tau+sqrt(3.0)+aux*pow(qemi,-1.0),-1.0);
    
    return result;
}