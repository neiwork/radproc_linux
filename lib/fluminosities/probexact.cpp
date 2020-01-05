#include "probexact.h"
#include <math.h>
#include <fmath/physics.h>
#include <fmath/RungeKutta.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>
#include <nrMath/integrators.h>
#include <nrMath/nrutil.h>
#include <nrMath/laguerre.h>

#include <iostream>
using namespace std;
#define zero 1.0e-7

double zedaminus(double rho, double eps, double g)
{
    double beta = sqrt(1.0-1.0/(g*g));
    double d = 1.0+eps-rho*eps;
    double aux = d*d-1.0/(g*g);
    return rho*(d-sqrt(aux));
}

double zedaplus(double rho, double eps, double g)
{
    double beta = sqrt(1.0-1.0/(g*g));
    double d = 1.0+eps-rho*eps;
    double aux = d*d-1.0/(g*g);
    return rho*(d+sqrt(aux));
}

double extrinf(double rho, void *params)
{
  struct two_d_params *p = (struct two_d_params *) params;
  double eps = p->p1;
  double g = p->p2;
  double beta = sqrt(1.0-1.0/(g*g));
  return zedaplus(rho,eps,g) - (1.0-beta);
}

double extrsup(double rho, void *params)
{
  struct two_d_params *p = (struct two_d_params *) params;
  double eps = p->p1;
  double g = p->p2;
  double beta = sqrt(1.0-1.0/(g*g));
  double zeda = zedaminus(rho,eps,g);
  
  return zedaminus(rho,eps,g) - (1.0+beta);
}

double f(double x, void *params)
{
	return (1.0-4.0/x-8.0/(x*x))*log(1.0+x)+0.5+8.0/x-0.5/((1.0+x)*(1.0+x));
}

double rateExact(double om, double gamma)
{
	gsl_function gsl_f;
	gsl_f.function = &f;
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

double rate(double om, double gamma)
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
					return rateExact(om,gamma);
				}
			}
		}
	}
}

double rateThermal(double normtemp, double energy)
{
	double alf=0.0;
	int nG=10;
    double *abscissas,*weights;
    abscissas=dvector(1,nG);
    weights=dvector(1,nG);
    gaulag(abscissas,weights,nG,alf);
	
	double om = energy/(electronMass*cLight2);
	double sum1 = 0.0;
	double sum2 = 0.0;
	for (int jG=1;jG<=nG;jG++) {
		double gamma = abscissas[jG]*normtemp+1.0;
        double cWeights = weights[jG]*gamma*sqrt(gamma*gamma-1.0);
		sum1 += (cWeights*rate(om,gamma));
		sum2 += cWeights;
	}
	return sum1/sum2;
}

double dPdz(double z,double om, double omprim, double gamma)
{
	double g2 = gamma*gamma;
    double beta = sqrt(1.0-1.0/g2);
	double zeda = 1.0-beta*z;
    double eps = omprim/gamma;
    double rho = om/omprim;
    double k = gamma/om;
    double k2 = k*k;
    double aux = rho*(beta*beta+eps*eps+2.0*beta*eps*z);
    double y0 = (eps+beta*z)*(rho+eps*rho-zeda)/aux;
    double aux1 = rho*rho*beta*beta;
    double aux2 = 2.0*rho*eps*(1.0-rho)*zeda;
    double aux3 = (rho-zeda)*(rho-zeda);
    double sqrt1 = aux1+aux2-aux3;
	if (sqrt1 < 0.0)
	{
		double sqrt1 = 0.0;
	}
    double delta = beta*sqrt(1.0-z*z)*sqrt(sqrt1) / aux;
    double a = zeda-(1.0-y0)/k;
    double a2 = a*a;
    double b = delta/k;
    double b2 = b*b;
    
    double num = 3.0*thomson*cLight*om;
    double r = rate(omprim,gamma);
	if (r <= 0.0) 
		return 0.0;
	
    double a2b2 = a2-b2;
    double minussqrta2b2 = 1.0/sqrt(a2b2);
	double msqrta2b2_3 = minussqrta2b2*minussqrta2b2*minussqrta2b2;
    
    double den = 16.0*g2*g2*omprim*omprim*r*sqrt(beta*beta+eps*eps+2.0*beta*eps*z)*zeda;
    
    double factor = 2.0*y0*k-a*k2+minussqrta2b2*(1.0+y0*y0-2.0*a*y0*k+a2*k2)+1.0/(k2*zeda) *
    (k2 + k2*minussqrta2b2*a*(2.0*b-a)/(a-b)+msqrta2b2_3*(a*(1.0-y0)*(1.0-y0)+2.0*k*
    b2*(1-y0)-b*a2*k2));
        
    double result = num*factor/den;
    return (result > 0.0 ? result : 0.0);
}

double dPdzR(double z,double om, double omprim, double gamma)
{
	double g2 = gamma*gamma;
    double beta = sqrt(1.0-1.0/g2);
	double zeda = 1.0-beta*z;
    double eps = omprim/gamma;
    double rho = om/omprim;
    double k = gamma/om;
    double k2 = k*k;
    double aux = rho*(beta*beta+eps*eps+2.0*beta*eps*z);
    double y0 = (eps+beta*z)*(rho+eps*rho-zeda)/aux;
    double aux1 = rho*rho*beta*beta;
    double aux2 = 2.0*rho*eps*(1.0-rho)*zeda;
    double aux3 = (rho-zeda)*(rho-zeda);
    double sqrt1 = aux1+aux2-aux3;
	if (sqrt1 < 0.0)
	{
		double sqrt1 = 0.0;
	}
    double delta = beta*sqrt(1.0-z*z)*sqrt(sqrt1) / aux;
    double a = zeda-(1.0-y0)/k;
    double a2 = a*a;
    double b = delta/k;
    double b2 = b*b;
    
    double num = 3.0*om;
	
    double a2b2 = a2-b2;
    double minussqrta2b2 = 1.0/sqrt(a2b2);
	double msqrta2b2_3 = minussqrta2b2*minussqrta2b2*minussqrta2b2;
    
    double den = 16.0*g2*g2*omprim*omprim*sqrt(beta*beta+eps*eps+2.0*beta*eps*z)*zeda;
    
    double factor = 2.0*y0*k-a*k2+minussqrta2b2*(1.0+y0*y0-2.0*a*y0*k+a2*k2)+1.0/(k2*zeda) *
    (k2 + k2*minussqrta2b2*a*(2.0*b-a)/(a-b)+msqrta2b2_3*(a*(1.0-y0)*(1.0-y0)+2.0*k*
    b2*(1-y0)-b*a2*k2));
        
    double result = num*factor/den;
    return (result > 0.0 ? result : 0.0);
}

double gsldPdz(double z, void *params)
{
	struct three_d_params *p = (struct three_d_params *) params;
	double om = p->p1;
	double omprim = p->p2;
	double gamma = p->p3;
	return dPdz(z,om,omprim,gamma);
}

double gsldPdzR(double z, void *params)
{
	struct three_d_params *p = (struct three_d_params *) params;
	double om = p->p1;
	double omprim = p->p2;
	double gamma = p->p3;
	return dPdzR(z,om,omprim,gamma);
}

double probexactNew(double om,double omprim,double gamma) 
{
	gsl_function gsl_dPdz;
	struct three_d_params dPdz_params = {om,omprim,gamma};
	gsl_dPdz.function = &gsldPdz;
	gsl_dPdz.params = &dPdz_params;
	
	double rho = om/omprim;
    double beta = sqrt(1.0-1.0/(gamma*gamma));
    double eps= omprim/gamma;
    double d=1.0+eps-eps*rho;
    double aux1=sqrt(d*d-1.0/(gamma*gamma));
    double zminus= (1.0 - rho*(d-aux1))/beta;
    double zplus=(1.0-rho*(d+aux1))/beta;
    double zmax = (1.0 > zminus ? zminus : 1.0);
	zmax *= (1.0-zero);
    double zmin= (-1.0 < zplus ? zplus : -1.0);
	zmin *= (1.0+zero);
	double error;
	int status;
	double result = 0.0;
	if (zmax > zmin)
		result = integrator_qags(&gsl_dPdz,zmin,zmax,0,1.0e-2,100,&error,&status);
	/*if (status)
		return 0.0;*/
	return (result > 0.0 ? result : 0.0);
}

double probexactR(double om,double omprim,double gamma) 
{
	gsl_function gsl_dPdzR;
	struct three_d_params dPdzR_params = {om,omprim,gamma};
	gsl_dPdzR.function = &gsldPdzR;
	gsl_dPdzR.params = &dPdzR_params;
	
	double rho = om/omprim;
    double beta = sqrt(1.0-1.0/(gamma*gamma));
    double eps= omprim/gamma;
    double d=1.0+eps-eps*rho;
    double aux1=sqrt(d*d-1.0/(gamma*gamma));
    double zminus= (1.0 - rho*(d-aux1))/beta;
    double zplus=(1.0-rho*(d+aux1))/beta;
    double zmax = (1.0 > zminus ? zminus : 1.0);
	zmax *= (1.0-zero);
    double zmin= (-1.0 < zplus ? zplus : -1.0);
	zmin *= (1.0+zero);
	double error;
	int status;
	double result = 0.0;
	if (zmax > zmin)
		result = integrator_qags(&gsl_dPdzR,zmin,zmax,0,1.0e-2,100,&error,&status);
	/*if (status)
		return 0.0;*/
	return (result > 0.0 ? result : 0.0);
}

double gslprobexact(double om, void *params)
{
	struct two_d_params *p = (struct two_d_params *) params;
	double omprim = p->p1;
	double gamma = p->p2;
	return probexactNew(om,omprim,gamma);
}

double probexact(double om,double omprim,double gamma) 
{
	double rho = om/omprim;
    double beta = sqrt(1.0-1.0/(gamma*gamma));
    double eps= omprim/gamma;
    double d=1.0+eps-eps*rho;
    double aux1=sqrt(d*d-1.0/(gamma*gamma));
    double zminus= (1.0 - rho*(d-aux1))/beta;
    double zplus=(1.0-rho*(d+aux1))/beta;
    double zmax = (1.0 > zminus ? zminus : 1.0);
	zmax -= zero;
    double zmin= (-1.0 < zplus ? zplus : -1.0);
	zmin += zero;
	double error;
	
	double result = 0.0;
	if (zmax > zmin)
		result = RungeKuttaSimple(zmin,zmax,[&om,&omprim,&gamma]
				(double z){return dPdz(z,om,omprim,gamma);});
	
	return (result > 0.0 ? result : 0.0);
}