#include "thermalCompton.h"
#include "probexact.h"

#include <nrMath/brent.h>
#include <nrMath/integrators.h>
#include <fmath/physics.h>
#include <fmath/RungeKutta.h>
#include <math.h>
#include <iostream>
using namespace std;

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

#include <nrMath/nrutil.h>
#include <nrMath/laguerre.h>


void comptonRedistribution(Vector& p, size_t nG, size_t nE, size_t nOm, 
							size_t jR, double normtemp, Vector energies, Matrix lum)
{
	void gaulag(double *x, double *w, int n, double alf);
	double *abscissas,*weights;
	abscissas=dvector(1,nG);
	weights=dvector(1,nG);
	double alf=0.0;
	gaulag(abscissas, weights,nG,alf);
	
	gsl_function gsl_extrinf;	gsl_extrinf.function = &extrinf;
			
	for (size_t jG=1;jG<=nG;jG++) {
		double gamma=abscissas[jG]*normtemp+1.0;
		double beta=sqrt(1.0-1.0/(gamma*gamma));
		for (size_t jE=0;jE<nE;jE++) {
			if (lum[jE][jR]*(energies[jE]/planck) > 1.0e10) {
				double omPrim = energies[jE]/(electronMass*cLight2);
				double eps = omPrim/gamma;
						
				struct two_d_params extr_params = {eps,gamma};
				gsl_extrinf.params = &extr_params;
				int status1,status2;
				double omMin = omPrim * brent(&gsl_extrinf,0.0,1.0,&status1,&status2);
				double omMaxAbs = omPrim + (gamma-1.0);
				double omMax = omPrim*(1.0+beta)/(1.0-beta+2.0*omPrim/gamma);
				double aux = beta/(1.0+gamma*(1.0+beta));
				omMax = (eps < aux) ? omMax : omMaxAbs;
				
				double probtotexact = RungeKuttaSimple(omMin,omMax,[&](double om)
					{return probexactNew(om,omPrim,gamma);});
				if (probtotexact > 0.0) {
					double var_int_om=pow(omMax/omMin,1.0/nOm);
					double om=omMin;
					for (size_t jOm=0;jOm<nOm;jOm++) {
						double dOm = om*(var_int_om-1.0);
						double prob = probexactNew(om,omPrim,gamma)*dOm;
						p[((jR*nG+(jG-1))*nE+jE)*nOm+jOm] = prob/probtotexact;
						om *= var_int_om;
					}
				}
			}
		}
	}
}

void comptonNew(Matrix& lumOut, double lumIn, Vector energies, 
				double temp, size_t nE, size_t jEprim, size_t jR) 
{
	gsl_function gsl_extrinf;
	gsl_extrinf.function = &extrinf;
	
    double omPrim=energies[jEprim]/(electronMass*cLight2);
    double var_int=pow(energies[nE-1]/energies[0],1.0/nE);
    double dnuprim=energies[jEprim]/planck * (var_int-1.0);
    double normtemp=boltzmann*temp/(electronMass*cLight2);
    size_t nOm=100;
    size_t nG=10;
    double alf=0.0;
    void gaulag(double *, double *, int, double);
    double *abscissas,*weights;
    abscissas=dvector(1,nG);
    weights=dvector(1,nG);
    double sum2=0.0;
    gaulag(abscissas, weights,nG,alf);
    Vector lumOut1(nE,0.0);
	
    for (size_t jG=1;jG<=nG;jG++) {
        double gamma=abscissas[jG]*normtemp+1.0;
        double beta=sqrt(1.0-1.0/(gamma*gamma));
        double cWeight = weights[jG]*gamma*sqrt(gamma*gamma-1.0);
		double eps = omPrim/gamma;
		struct two_d_params extrinf_params = {eps,gamma};
		gsl_extrinf.params = &extrinf_params;
		int status1,status2;
		double omMin = omPrim * brent(&gsl_extrinf,0.0,1.0,&status1,&status2);
		double omMaxAbs = omPrim + (gamma - 1.0);
		double omMax = omPrim*(1.0+beta)/(1.0-beta+2.0*omPrim/gamma);
		double aux = beta/(1.0+gamma*(1.0+beta));
		omMax = (eps < aux) ? omMax : omMaxAbs;
		
		double probtotexact = RungeKuttaSimple(omMin,omMax,[&](double om)
		{return probexactNew(om,omPrim,gamma);});
		double error;
		
		double var_int_om=pow(omMax/omMin,1.0/nOm);
		double om=omMin;
		size_t jE=1;
		if (probtotexact > 0.0) {

		for (size_t jOm=1;jOm<nOm;jOm++) {
			om *= var_int_om;
			double dOm=om*(var_int_om-1.0);
			double prob = probexactNew(om,omPrim,gamma);//probtotexact;
			if (prob < 0.0)
				cout << "ERROR" << endl;
			
			double lim1=sqrt(energies[1]*energies[0]);
			double lim2;
			double energy = om * (electronMass*cLight2);
			
			double lumOutAux = lumIn * prob * (om/omPrim);
			int count=0;
			while (count == 0 && jE < nE-2) {
				lim2 = sqrt(energies[jE]*energies[jE+1]);
				if ( energy < lim2 && energy > lim1) {
					double dnu = (lim2-lim1)/planck;
					lumOut1[jE] += lumOutAux * cWeight * dOm/dnu;
					count++;
				} else {
					lim1=lim2;
					jE++;
				}
			}
		}}
		sum2 += cWeight;
	}
    for (size_t jE=0;jE<nE;jE++) {
        lumOut[jE][jR] += lumOut1[jE] * dnuprim /sum2;
    }
}

void comptonNew2(Matrix& lumOut, double lumIn, Vector energies, size_t nG, size_t nE,
				size_t nOm, double normtemp, Vector probab, size_t jEprim, size_t jR) 
{
	gsl_function gsl_extrinf;
	gsl_extrinf.function = &extrinf;
	
    double omPrim=energies[jEprim]/(electronMass*cLight2);
    double var_int=pow(energies[nE-1]/energies[0],1.0/nE);
    double dnuprim=energies[jEprim]/planck * (var_int-1.0);
    double alf=0.0;
    double *abscissas,*weights;
    abscissas=dvector(1,nG);
    weights=dvector(1,nG);
    double sum2=0.0;
    gaulag(abscissas,weights,nG,alf);
    Vector lumOut1(nE,0.0);
	
    for (size_t jG=1;jG<=nG;jG++) {
        double gamma = abscissas[jG]*normtemp+1.0;
        double beta = sqrt(1.0-1.0/(gamma*gamma));
        double cWeights = weights[jG]*gamma*sqrt(gamma*gamma-1.0);
		double eps = omPrim/gamma;
		
		struct two_d_params extrinf_params = {eps,gamma};
		gsl_extrinf.params = &extrinf_params;
		int status1,status2;
		
		double omMin = omPrim * brent(&gsl_extrinf,0.0,1.0,&status1,&status2);
		double omMaxAbs = omPrim + (gamma - 1.0);
		double omMax = omPrim*(1.0+beta)/(1.0-beta+2.0*omPrim/gamma);
		double aux = beta/(1.0+gamma*(1.0+beta));
		omMax = (eps < aux) ? omMax : omMaxAbs;
		double error;
		
		double var_int_om=pow(omMax/omMin,1.0/nOm);
		double om=omMin;
		
		size_t jE=1;
		for (size_t jOm=0;jOm<nOm;jOm++) {
			double lumOutAux = lumIn * double(probab[((jR*nG+jG)*nE+jE)*nOm+jOm]);//*(om/omPrim);
			double lim1=sqrt(energies[1]*energies[0]);
			double lim2;
			double energy = om * (electronMass*cLight2);
			int count=0;
			while (count == 0 && jE < nE-1) {
				lim2 = sqrt(energies[jE]*energies[jE+1]);
				if ( energy < lim2 && energy > lim1) {
					double dnu = (lim2-lim1)/planck;
					lumOut1[jE] += lumOutAux * cWeights /dnu;
					count++;
				} else {
					lim1=lim2;
					jE++;
				}
			}
			om *= var_int_om;
		}
		sum2 += cWeights;
	}
	
    for (size_t jE=0;jE<nE;jE++) {
        lumOut[jE][jR] += lumOut1[jE] * dnuprim /sum2;
    }
}