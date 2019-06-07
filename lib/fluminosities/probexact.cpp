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

double gsldPdz(double z, void *params)
{
	struct three_d_params *p = (struct three_d_params *) params;
	double om = p->p1;
	double omprim = p->p2;
	double gamma = p->p3;
	return dPdz(z,om,omprim,gamma);
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

double probInterpolated(Vector comp, double nu, double nuPrim, double gamma)
{
	int nGamma, nNu, nNuPrim;
	nGamma = nNu = nNuPrim = 30;
	double gammaMin = 1.0;
	double gammaMax = 5.0e3;
	double pasoGamma = pow(gammaMax/gammaMin,1.0/nGamma);
	double nuPrimMin = 1.0e8;
	double nuPrimMax = 1.0e20;
	double pasoNuPrim = pow(nuPrimMax/nuPrimMin,1.0/nNuPrim);
	double nuMin = 1.0e10;
	double nuMax = 1.0e22;
	double pasoNu = pow(nuMax/nuMin,1.0/nNu);

	size_t kGamma, kNuPrim, kNu;
	kGamma = kNuPrim = kNu = 0;
	double gammaAux = gammaMin;
	while (gammaAux < gamma) {
		gammaAux *= pasoGamma;
		kGamma++;
	}
	double nuPrimAux = nuPrimMin;
	while (nuPrimAux < nuPrim) {
		nuPrimAux *= pasoNuPrim;
		kNuPrim++;
	}
	double nuAux = nuMin;
	while (nuAux < nu) {
		nuAux *= pasoNu;
		kNu++;
	}
	
	double logGamma1 = log(gammaAux/pasoGamma);
	double logGamma2 = log(gammaAux);
	double logGamma = log(gamma);
	double logNuPrim1 = log(nuPrimAux/pasoNuPrim);
	double logNuPrim2 = log(nuPrimAux);
	double logNuPrim = log(nuPrim);
	double logNu1 = log(nuAux/pasoNu);
	double logNu2 = log(nuAux);
	double logNu = log(nu);
	
	double prob111 = comp[((kGamma-1)*nNuPrim+(kNuPrim-1))*nNu+(kNu-1)];
	double prob112 = comp[((kGamma-1)*nNuPrim+(kNuPrim-1))*nNu+kNu];
	double prob121 = comp[((kGamma-1)*nNuPrim+kNuPrim)*nNu+(kNu-1)];
	double prob122 = comp[((kGamma-1)*nNuPrim+kNuPrim)*nNu+kNu];
	double prob211 = comp[(kGamma*nNuPrim+(kNuPrim-1))*nNu+(kNu-1)];
	double prob212 = comp[(kGamma*nNuPrim+(kNuPrim-1))*nNu+kNu];
	double prob221 = comp[(kGamma*nNuPrim+kNuPrim)*nNu+(kNu-1)];
	double prob222 = comp[(kGamma*nNuPrim+kNuPrim)*nNu+kNu];
	
	double prob11 = (prob112-prob111)/(logNu2-logNu1) * (logNu-logNu1) + prob111;
	double prob12 = (prob122-prob121)/(logNu2-logNu1)*(logNu-logNu1) + prob121;
	double prob1 = (prob12-prob11)/(logNuPrim2-logNuPrim1)*(logNuPrim-logNuPrim1)+prob11;
	
	double prob21 = (prob212-prob211)/(logNu2-logNu1) * (logNu-logNu1) + prob211;
	double prob22 = (prob222-prob221)/(logNu2-logNu1)*(logNu-logNu1) + prob221;
	double prob2 = (prob22-prob21)/(logNuPrim2-logNuPrim1)*(logNuPrim-logNuPrim1)+prob21;
	
	double prob = (prob2-prob1)/(logGamma2-logGamma1)*(logGamma-logGamma1) + prob1;
	return prob;
}

double probInterpolated2(Vector comp, double nu, double nuPrim, double temp)
{
	int nTemp, nNu, nNuPrim;
	nTemp = nNu = nNuPrim = 30;
	double tempMin = 0.1;
	double tempMax = 10.0;
	double pasoTemp = pow(tempMax/tempMin,1.0/nTemp);
	double nuPrimMin = 1.0e8;
	double nuPrimMax = 1.0e20;
	double pasoNuPrim = pow(nuPrimMax/nuPrimMin,1.0/nNuPrim);
	double nuMin = 1.0e10;
	double nuMax = 1.0e22;
	double pasoNu = pow(nuMax/nuMin,1.0/nNu);

	size_t kTemp, kNuPrim, kNu;
	kTemp = kNuPrim = kNu = 0;
	double tempAux = tempMin;
	while (tempAux < temp) {
		tempAux *= pasoTemp;
		kTemp++;
	}
	double nuPrimAux = nuPrimMin;
	while (nuPrimAux < nuPrim) {
		nuPrimAux *= pasoNuPrim;
		kNuPrim++;
	}
	double nuAux = nuMin;
	while (nuAux < nu) {
		nuAux *= pasoNu;
		kNu++;
	}
	
	double logTemp1 = log(tempAux/pasoTemp);
	double logTemp2 = log(tempAux);
	double logTemp = log(temp);
	double logNuPrim1 = log(nuPrimAux/pasoNuPrim);
	double logNuPrim2 = log(nuPrimAux);
	double logNuPrim = log(nuPrim);
	double logNu1 = log(nuAux/pasoNu);
	double logNu2 = log(nuAux);
	double logNu = log(nu);
	
	double prob111 = comp[((kTemp-1)*nNuPrim+(kNuPrim-1))*nNu+(kNu-1)];
	double prob112 = comp[((kTemp-1)*nNuPrim+(kNuPrim-1))*nNu+kNu];
	double prob121 = comp[((kTemp-1)*nNuPrim+kNuPrim)*nNu+(kNu-1)];
	double prob122 = comp[((kTemp-1)*nNuPrim+kNuPrim)*nNu+kNu];
	double prob211 = comp[(kTemp*nNuPrim+(kNuPrim-1))*nNu+(kNu-1)];
	double prob212 = comp[(kTemp*nNuPrim+(kNuPrim-1))*nNu+kNu];
	double prob221 = comp[(kTemp*nNuPrim+kNuPrim)*nNu+(kNu-1)];
	double prob222 = comp[(kTemp*nNuPrim+kNuPrim)*nNu+kNu];
	
	double prob11 = (prob112-prob111)/(logNu2-logNu1) * (logNu-logNu1) + prob111;
	double prob12 = (prob122-prob121)/(logNu2-logNu1)*(logNu-logNu1) + prob121;
	double prob1 = (prob12-prob11)/(logNuPrim2-logNuPrim1)*(logNuPrim-logNuPrim1)+prob11;
	
	double prob21 = (prob212-prob211)/(logNu2-logNu1) * (logNu-logNu1) + prob211;
	double prob22 = (prob222-prob221)/(logNu2-logNu1)*(logNu-logNu1) + prob221;
	double prob2 = (prob22-prob21)/(logNuPrim2-logNuPrim1)*(logNuPrim-logNuPrim1)+prob21;
	
	double prob = (prob2-prob1)/(logTemp2-logTemp1)*(logTemp-logTemp1) + prob1;
	return prob;
}

double probTemp(Vector comp, double nu, double nuPrim, double normtemp)
{
	void gaulag(double *x, double *w, int n, double alf);
	size_t nG = 5;
	double *abscissas,*weights;
	abscissas=dvector(1,nG);
	weights=dvector(1,nG);
	double alf=0.0;
	gaulag(abscissas, weights,nG,alf);
	double sum1=0.0;
	double sum2=0.0;
	double omPrim = planck*nuPrim/(electronMass*cLight2);
	
	for (size_t jG=1;jG<=nG;jG++) {
		double gamma=abscissas[jG]*normtemp+1.0;
		double cWeights = weights[jG]*gamma*sqrt(gamma*gamma-1.0);
		sum1 += probInterpolated(comp,nu,nuPrim,gamma)*cWeights;
		sum2 += cWeights;
	}
	
	return sum1/sum2;
}