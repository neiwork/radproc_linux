#include <math.h>
#include <fmath/laguerre.h>
#include <fmath/mathFunctions.h>
#include <fparameters/nrutil.h>
#include <fmath/physics.h>
#include "coppiProbDist.h"

double rate(double frec, double gamma)
{
	double beta=sqrt(1.0-1.0/(gamma*gamma));
	double cte=3.0*cLight*thomson/(32.0*(gamma*frec)*(gamma*frec)*beta);
	double min=2.0*gamma*frec*(1.0-beta);
	double max=2.0*gamma*frec*(1.0+beta);
	int nX=100;
	double sum=0.0;
	double x_int=pow(max/min,1.0/nX);
	double x = min;
	for (int i=0;i<nX;++i)   
	{
		double dx = x*(x_int-1.0);
		double f=(1.0-4.0/x-8.0/(x*x))*log(1.0+x)+0.5+8.0/x-0.5/((1.0+x)*(1.0+x));
		sum = sum + f*dx;
		x = x*x_int;
	}
	return cte*sum;
}

double frecprom(double frecprim, double gamma)  //<w>
{
	double beta=sqrt(1.0-1.0/(gamma*gamma));
	double min = -1.0;
	double max = 1.0;
	int nMu = 100;

	double sum = 0.0;
	double mu_int = (max-min)/nMu;
	double mu = min;
	for (int i=0; i<nMu;++i)
	{
		double x = gamma*frecprim*(1.0-beta*mu);
		double a = gamma*gamma*beta*frecprim*(mu-beta);
		double factor1=3.0*cLight*thomson*x/(8.0*gamma*frecprim*rate(frecprim,gamma));
		double factor2=2.0*(a+gamma*(x-2.0)-(2.0*gamma+a)/x-a*(6.0/(x*x)+3.0/(x*x*x)))/(1.0+2.0*x);
		double factor3=log(1.0+2.0*x)*(gamma-a+3.0*a*(1.0/x+1.0/(x*x)))/(x*x);
		double factor4=(1.0 - 1.0/((1.0+2.0*x)*(1.0+2.0*x)))*(a*(3.0/(x*x)+1.0/(x*x*x))
                                    +(a+gamma)/x +2.0*gamma)/(2.0*x);
		double factor5=(1.0-1.0/P3(1.0+2.0*x))/3.0 *(gamma+a*(1.0/x+1.0/(x*x)))/3.0-2.0*a/(x*x*x);
		double f=factor1*(factor2+factor3+factor4+factor5);
		sum+=f*mu_int/2.0;
		mu+=mu_int;
	}
    return sum;
}

double sqrdfrecprom(double frecprim, double gamma)  //<w^2>
{
	double beta = sqrt(1.0-1.0/(gamma*gamma));
    double gamma2=gamma*gamma;
    double frecprim2=frecprim*frecprim;
	if(frecprim*gamma < 1.0)
	{
		return 14.0*frecprim2*(gamma2*gamma2)*(1.0-176.0*gamma*frecprim/35.0)/5.0; 
	}
	else
	{
		double a1 = 64.0*beta*(gamma2*frecprim2);
		double a2 = (6.0*gamma*frecprim*(1.0+beta)*(2.0*gamma2+1.0)+6.0*gamma2+3.0)*
                                log(2.0*gamma*frecprim*(1.0+beta)+1.0);
		double a3 = (6.0*gamma*frecprim*(1.0-beta)*(2.0*gamma2+1.0)+6.0*gamma2+3.0)*
                                log(2.0*gamma*frecprim*(1.0-beta)+1.0);
		double a4 = 9.0*(1.0-beta*beta)*(1.0-beta*beta)*gamma*gamma2*frecprim/32.0- 
                            (58.0*gamma2+1.0)/(64.0*gamma*frecprim)+7.0*(2.0*gamma2*(1.0-beta*beta)+1.0)/32.0;

		return (a2 - a3 + a4)/rate(frecprim,gamma)/a1;
	}
}

double probability(double frec, double frecprim, double gamma)
{
    double frecave=frecprom(frecprim,gamma);
    double varfrec=sqrdfrecprom(frecprim,gamma)-frecave*frecave;
    double d=new_min(sqrt(3.0*varfrec),frecave,varfrec);
    double x=d-abs(frec- frecave);
    return (x >= 0.0 ? 0.5/d : 0.0);
}

double k(double frec,double frecprim,double normtemp)
{
    int N=6;
    double alf=0.0;
    void gaulag(double *x, double *w, int n, double alf);
    double *abscissas,*weights,*gamma,*c;
    abscissas=dvector(1,N);
    weights=dvector(1,N);
    gamma=dvector(1,N);
    c=dvector(1,N);
    double sum1=0.0,sum2=0.0;
    gaulag(abscissas, weights,N,alf);
    for (int i=1;i<=N;i++) {
        gamma[i]=abscissas[i]*normtemp+1.0;
        c[i]=weights[i]*gamma[i]*sqrt(gamma[i]*gamma[i]-1.0);
        sum1+=c[i]*probability(frec,frecprim,gamma[i])*rate(frecprim,gamma[i]);
        sum2+=c[i];
    }
    return sum1/sum2;
}