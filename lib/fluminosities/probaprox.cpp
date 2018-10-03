#include "probaprox.h"
#include <math.h>
#include <fmath/physics.h>
#include <fmath/mathFunctions.h>

double rate1(double om, double gamma)
{
	double beta=sqrt(1.0-1.0/(gamma*gamma));
    
    if (beta < 1.0e-2 && om < 1.0e-2) {
        return cLight*thomson*(1.0-2.0*gamma*om);
    } else {
        double cte=3.0*cLight*thomson/(32.0*(gamma*om)*(gamma*om)*beta);
        double min=2.0*gamma*om*(1.0-beta);
        double max=2.0*gamma*om*(1.0+beta);
        int nX=100;
        double sum=0.0;
        double x_int=pow(max/min,1.0/nX);
        double x = min;
        for (int i=0;i<nX;++i)   
        {
            double dx = x*(x_int-1.0);
            double f=(1.0-4.0/x-8.0/(x*x))*log(1.0+x)+0.5+8.0/x-0.5/((1.0+x)*(1.0+x));
            sum += f*dx;
            x = x*x_int;
        }
        double result = cte*sum;
        return result;
    }
}

/*double rate(double om, double gamma)
{
    double x = gamma*om;
    if (x > 0.02 && x < 900) {
        double cte = 3.0*cLight/(8.0*x);
        double result = (1.0-2.0/x-2.0/(x*x))*log(1.0+2.0*x)+0.5+4.0/x-0.5/P2(1.0+2.0*x);
        return cte*result;
    } else {
        return 0.0;
    }
}*/

double omprom(double omprim, double gamma)  //<w>
{
	double beta=sqrt(1.0-1.0/(gamma*gamma));
	double result;
    if (beta < 1.0e-2 && omprim < 1.0e-3) {
        result= (1.0 + 4.0/3.0 * beta*beta - omprim)*omprim;
    } else {
        double min = -0.99;
        double max = 0.99;
        int nMu = 10;
        double sum = 0.0;
        double dmu = (max-min)/nMu;
        double mu = min;
        double r=rate1(omprim,gamma);
        for (int i=0; i<nMu;++i)
        {
            double x = gamma*omprim*(1.0-beta*mu);
            double a = gamma*gamma*beta*omprim*(mu-beta);
            double factor1=3.0*cLight*thomson*x/(8.0*gamma*omprim*r);
            double factor2=2.0/(1.0+2.0*x)*(a+gamma*(x-2.0)-(2.0*gamma+a)/x-a*(6.0/(x*x)+3.0/(x*x*x)));
            double factor3=log(1.0+2.0*x)*(gamma-a+3.0*a*(1.0/x+1.0/(x*x)))/(x*x);
            double factor4=(1.0 - 1.0/((1.0+2.0*x)*(1.0+2.0*x)))*(a*(3.0/(x*x)+1.0/(x*x*x))
                                    +(a+gamma)/x +2.0*gamma)/(2.0*x);
            double factor5=(1.0-1.0/P3(1.0+2.0*x))/3.0 *(gamma+a*(1.0/x+1.0/(x*x)))-2.0*a/(x*x*x);
            double f=factor1*(factor2+factor3+factor4+factor5);
            sum += f*dmu;
            mu += dmu;
        }
        result= sum/2.0;
    }
    return result;
}

double sqrdauxom(double omprim, double gamma, double beta, double alpha, double mu)
{
    double aux1=P2(gamma*(1.0-beta*mu)+gamma*beta*alpha*(mu-beta))+(0.5*beta*beta*(1.0-alpha*alpha)*
                    (1.0-mu*mu));
    double aux2=P2(1.0+gamma*(1.0-beta*mu)*(1.0-alpha)*omprim);
    double result = gamma*gamma*omprim*omprim * aux1/aux2;
    return result;
}

double pmualpha(double omprim, double gamma, double beta, double alpha, double mu)
{
    double x = gamma*omprim*(1.0-beta*mu);
    double xprim = x/(1.0+(1.0-alpha)*x);
    double y=xprim/x;
    double dom_da = 3.0/8.0 * thomson * y*y * (1.0/y + y - 1.0 + alpha*alpha);
    double result=cLight*(1.0-beta*mu)*dom_da/rate1(omprim,gamma);
    return result;
}

double sqrdomprom(double omprim, double gamma)  //<w^2>
{
    double beta = sqrt(1.0-1.0/(gamma*gamma));
    int nAlpha=10;
    int nMu=10;
    double max=0.99;
    double min=-0.99;
    double dmu=(max-min)/nMu;
    double dalpha=(max-min)/nAlpha;
    
    double mu = min;
    double sum=0.0;
    for (int i=0;i<nMu;i++) {
        double alpha = min;
        for (int j=0;j<nAlpha;j++) {
            double auxom2 = sqrdauxom(omprim,gamma,beta,alpha,mu);
            double pmua = pmualpha(omprim,gamma,beta,alpha,mu);
            sum += dmu*dalpha*auxom2*pmua;
            alpha += dalpha;
        }
        mu += dmu;
    }
    return sum/2.0;
}

double probaprox(double om, double omprim, double gamma,double ommin,double ommax)
{
    double omave=omprom(omprim,gamma);
    double varom=sqrdomprom(omprim,gamma)-omave*omave;
    double A=sqrt(3.0*varom);
    double d=new_min(A,omave-ommin,ommax-omave);
    double x=d-abs(om- omave);
    double result = (x >= 0.0 ? 0.5/d : 0.0);
    return result;
}