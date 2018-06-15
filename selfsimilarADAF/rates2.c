#include <math.h>
#include "ssfunctions2.h"
#include "rates2.h"
#include "nr.h"
#define TOL 1.0e-3
#define PI 3.14159265359
#define BOLTZMANN 1.38e-16
#define PLANCK 6.6260755e-27
#define ELECTRONMASS 9.109e-28
#define PROTONMASS 1.67262e-24
#define CLIGHT2 8.9875518e20
#define SIGMASB 5.6704e-5

double qbremss(double temp)
{
    extern double eDensity,iDensity;
    double fei,fee,qei,qee;
  
    double tempe=(BOLTZMANN*temp)/(ELECTRONMASS*CLIGHT2);
    if (tempe <= 1.0) {
        fei=1.0159*sqrt(tempe)*(1.0+1.781*pow(tempe,1.34));
        fee=pow(tempe,1.5)*(1.0+1.1*tempe+tempe*tempe-1.25*pow(tempe,2.5));
    } else {
        fei=1.4324*tempe*(log(1.123*tempe+0.48)+1.5);
        fee=tempe*(log(1.1232*tempe)+1.28);
    }
    qei=1.48e-22*eDensity*iDensity*fei;
    qee=3.4e-22*eDensity*eDensity*fee;
    return qei+qee;
}

double auxf1(double x)
{
    extern double radius,eDensity,magField,eTemp;
    return exp(1.8899*pow(x,1.0/3.0))-2.49e-10*4.0*PI*eDensity*radius/magField/
        (eTemp*eTemp*eTemp*bessk(2,1.0/eTemp)*pow(x,5.0/3.0))*(0.5316+0.4*pow(x,0.25)+sqrt(x));
}

#include "nr.h"
void qsync(double *qSy, double *nucrit, double tempe)
{
    int zbrac(double (*func)(double),double *x1, double *x2);
    double zbrent(double (*func)(double), double x1, double x2, double tol);
    extern double radius,magField;
  
    double normtempe=(BOLTZMANN*tempe)/(ELECTRONMASS*CLIGHT2);
    double x1=1.0e-5;
    double x2=1.0e3;
    zbrac(auxf1,&x1,&x2);
    double xcrit=zbrent(auxf1,x1,x2,TOL);
    double nu0=2.8e6*magField;                                                        // [Hz]
    double auxnucrit=1.5*nu0*normtempe*normtempe*xcrit;
    *nucrit=auxnucrit;
    //printf("nu0 = %f, nuc = %f, r = %f\n",nu0,nucrit,r);
    double auxqSy=2.0*PI/3.0/CLIGHT2*BOLTZMANN*tempe*(*nucrit)*(*nucrit)*(*nucrit)/radius;
    *qSy=auxqSy;
}

void qC(double qbr, double qsy, double nuc, double temp, double *qCbr, double *qCsy, double f)
{
    double prob,a,eta1,eta2,eta3;
    double xc=PLANCK*nuc/(ELECTRONMASS*CLIGHT2);
    double tempe=(BOLTZMANN*temp)/(ELECTRONMASS*CLIGHT2);
    prob=1.0-exp(-taues(f));
    a=1.0+4.0*tempe*(1.0+4.0*tempe);
    eta1=prob*(a-1.0)/(1.0-prob*a);
    eta3=-1.0-log(prob)/log(a);
    eta2=pow(3.0,-eta3)*eta1;
    double auxqCbr=qbr*tempe*eta1*((1.0-xc/tempe)-3.0/(eta3+1.0) * (pow(1.0/3.0,eta3+1.0)-
						 pow(xc/(3.0*tempe),eta3+1.0)));
    *qCbr=auxqCbr;
    double auxqCsy=qsy*(eta1-eta2*pow(xc/tempe,eta3));
    *qCsy=auxqCsy;
}

double qiefunc2(double tempions, double tempelectrons, double f)
{
    extern double r;
    double rad=r;
    double tempi=(BOLTZMANN*tempions)/(PROTONMASS*CLIGHT2);
    double tempe=(BOLTZMANN*tempelectrons)/(ELECTRONMASS*CLIGHT2);
    double aux1=(tempi+tempe);
    double aux3=1.0/tempi;
    double aux4=1.0/tempe;
    double aux2=aux3+aux4;
    double besselk1,besselk0,besselk2e,besselk2i,result;
    if (aux2 > 200.0) {
        if (aux3 > 100.0) {
            besselk2e=bessk(2,aux4);
            //result=3.328e-22*ne(f)*ni(f)*(1836.15*tempi-tempe)*sqrt(aux3/aux2)
            //    * ((2.0*aux1*aux1+1.0)/aux1 + 2.0) * exp(-(aux2+aux3))/besselk2e;
            result=5.61e-32*ne(f)*ni(f)*(tempions-tempelectrons)*sqrt(aux3/aux2)
                * ((2.0*aux1*aux1+1.0)/aux1 + 2.0) * exp(-(aux2-aux3))/besselk2e;
        } else {
            besselk2i=bessk(2,aux3);
            //result=3.328e-22*ne(f)*ni(f)*(1836.15*tempi-tempe)*sqrt(aux4/aux2)
            //       * ((2.0*aux1*aux1+1.0)/aux1 + 2.0) * exp(-(aux2+aux4))/besselk2i;
            result=5.61e-32*ne(f)*ni(f)*(tempions-tempelectrons)*sqrt(aux4/aux2)
                   * ((2.0*aux1*aux1+1.0)/aux1 + 2.0) * exp(-(aux2-aux4))/besselk2i;
        }
    } else {
        besselk1=bessk1(aux2);
        besselk0=bessk0(aux2);
        besselk2i=bessk(2,aux3);
        besselk2e=bessk(2,aux4);
        //result=3.328e-22*ne(f)*ni(f)*(1836.15*tempi-tempe)*
        //    ((2.0*aux1*aux1+1.0)/aux1 * besselk1/besselk2i + 2.0*besselk0/besselk2i) / besselk2e;
        result=5.61e-32*ne(f)*ni(f)*(tempions-tempelectrons)*
            ((2.0*aux1*aux1+1.0)/aux1 * besselk1/besselk2i + 2.0*besselk0/besselk2i) / besselk2e;
    }
    
    return result;
}

double qem2(double tempe, double f)
{
    extern double eDensity,iDensity,magField,eTemp;
    double qBr,qSy,qCbr,qCsy,nucrit;
    eDensity=ne(f);
    iDensity=ni(f);
    magField=magf(f);
    eTemp=(BOLTZMANN*tempe)/(ELECTRONMASS*CLIGHT2);
    qsync(&qSy,&nucrit,tempe);
    qBr=qbremss(tempe);
    qC(qBr,qSy,nucrit,tempe,&qCbr,&qCsy,f);
    double qemi=qBr+qSy+qCbr+qCsy;
    double aux=4.0*SIGMASB*pow(tempe,4)/height(f);
    double tauabs=qemi/aux;
    double tau=tauabs+taues(f);
    double result=aux*pow(1.5*tau+sqrt(3.0)+aux*pow(qemi,-1.0),-1.0);
    
    return qemi;
}