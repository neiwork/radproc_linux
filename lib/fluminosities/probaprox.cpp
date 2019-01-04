#include "probaprox.h"
#include <math.h>
#include <fmath/physics.h>
#include <fmath/mathFunctions.h>
#include <fmath/integrators.h>
#include <gsl/gsl_math.h>

#define zero 1.0e-7

double f2(double x, void *params)
{
	return (1.0-4.0/x-8.0/(x*x))*log(1.0+x)+0.5+8.0/x-0.5/((1.0+x)*(1.0+x));
}

double rateExact2(double om, double gamma)
{
	gsl_function gsl_f;
	gsl_f.function = &f2;
	double error;
	int status;
	double beta = sqrt(1.0-1.0/(gamma*gamma));
    double y = gamma*om;
    double cte = 3.0*cLight*thomson/(32.0*y*y*beta);
    double min = 2.0*y*(1.0-beta)*(1.0+zero);
    double max = 2.0*y*(1.0+beta)*(1.0-zero);
	void *params;
	size_t numevals;
    double result = cte * integrator_qags(&gsl_f,min,max,0,1.0e-2,100,&error,&status);
	/*if (status)
		return 0.0;*/
	return (result > 0.0 ? result : 0.0);
};

double rate2(double om, double gamma)
{
	double f(double);
	double beta=sqrt(1.0-1.0/(gamma*gamma));
	double y=gamma*om;
    double cte = 3.0*cLight*thomson/(8.0*y);
	
    if (beta < 1.0e-3 && om < 1.0e-3) {
        return cLight*thomson*(1.0-2.0*y);
    } else {
		if (om < 1.0/(2.0*gamma*(1.0+beta))) {
			return cLight*thomson*(1.0-2.0*y/3.0 * (3.0 +beta*beta));
		} else {
			if (2.0 < gamma < 30.0 && 0.01 < om < 30.0) {
				return cte*( (1.0-2.0/y-2.0/(y*y))*log(1.0+2.0*y)+0.5+4.0/y-0.5/P2(1.0+2.0*y) );
			} else {
				if (om > 1.0e2/(4.0*gamma) && gamma > 1.0e2) {
					return cte*log(4.0*y);
				} else {
					return rateExact2(om,gamma);
				}
			}
		}
	}
}

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
        double r=rate2(omprim,gamma);
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

double sqrdauxom(double omprim, double gamma, double alpha, double mu)
{
	double beta = sqrt(1.0-1.0/(gamma*gamma));
	double aux1 = gamma*(1.0-beta*mu)+gamma*beta*alpha*(mu-beta);
	aux1 *= aux1;
    double aux2 = 0.5*beta*beta*(1.0-alpha*alpha)*(1.0-mu*mu);
    double aux3=1.0+gamma*(1.0-beta*mu)*(1.0-alpha)*omprim;
	aux3 *= aux3;
    return gamma*gamma*omprim*omprim * (aux1+aux2)/aux3;
}

double om2pmualpha(double alpha, void *params)
{
	struct three_d_params *p = (struct three_d_params *) params;
	double mu = p->p1;
	double omprim = p->p2;
	double gamma = p->p3;
	
	double om2 = sqrdauxom(omprim,gamma,alpha,mu);
	double beta = sqrt(1.0-1.0/(gamma*gamma));
    double x = gamma*omprim*(1.0-beta*mu);
    double xprim = x/(1.0+(1.0-alpha)*x);
    double y=xprim/x;
    double dom_da = 3.0/8.0 * thomson * y*y * (1.0/y + y - 1.0 + alpha*alpha);
    double result=cLight*(1.0-beta*mu)*dom_da/rate2(omprim,gamma);
    return result*om2;
}

double intpmualpha(double mu, void *params)
{
	double error;
	int status;
	
	struct two_d_params *p = (struct two_d_params *) params;
	double omprim = p->p1;
	double gamma = p->p2;
	
	struct three_d_params om2pmualpha_params = {mu,omprim,gamma};
	gsl_function gsl_om2pmualpha;
		gsl_om2pmualpha.function = &om2pmualpha;
		gsl_om2pmualpha.params = &om2pmualpha_params;
		
	return integrator_qags(&gsl_om2pmualpha,-2.0,1.0,0,1.0e-2,100,&error,&status);
}

double sqrdomprom(double omprim, double gamma)  //<w^2>
{
    double error;
	int status;
    
	struct two_d_params intpmualpha_params = {omprim,gamma};
	gsl_function gsl_intpmualpha;
		gsl_intpmualpha.function = &intpmualpha;
		gsl_intpmualpha.params = &intpmualpha_params;
	
	double sum = integrator_qags(&gsl_intpmualpha,-2.0,1.0,0,1.0e-2,100,&error,&status);
    return sum/2.0;
}

double probaprox(double om, double omprim, double gamma,double ommin,double ommax)
{
    //double omave=omprom(omprim,gamma);
	double omave = (4.0/3.0)*gamma*gamma*omprim;
    double varom=sqrdomprom(omprim,gamma)-omave*omave;
    double A=sqrt(3.0*varom);
    double d=new_min(A,omave-ommin,ommax-omave);
    double x=d-abs(om- omave);
    double result = (x >= 0.0 ? 0.5/d : 0.0);
    return result;
}