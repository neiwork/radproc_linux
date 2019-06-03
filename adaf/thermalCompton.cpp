#include "thermalCompton.h"
#include <fluminosities/probexact.h>

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
#include <fstream>

double probInt(double om, void* params) {
	struct two_d_params *p = (struct two_d_params *) params;
	double omPrim = p->p1;
	double gamma = p->p2;
	
	return probexactNew(om,omPrim,gamma);
}

void comptonMatrix()
{
	FILE *fileBin;	fileBin=fopen("comptonRedMatrix.bin","wb");
	ofstream file1;	file1.open("comptonProbMatrix.dat",ios::out);
	ofstream file2;	file2.open("comptonProbTot.dat",ios::out);
	
	size_t nGamma, nNuPrim, nNu;
	nGamma = nNuPrim = nNu = 30;
	
	float *comptonprobvector;
	comptonprobvector=(float*)calloc(nGamma*nNuPrim*nNu,sizeof(float));
	
	double gammaMin = 1.0;
	double gammaMax = 5.0e3;
	double pasoGamma = pow(gammaMax/gammaMin,1.0/nGamma);
	double nuPrimMin = 1.0e8;
	double nuPrimMax = 1.0e20;
	double pasoNuPrim = pow(nuPrimMax/nuPrimMin,1.0/nNuPrim);
	double nuMin = 1.0e10;
	double nuMax = 1.0e22;
	double pasoNu = pow(nuMax/nuMin,1.0/nNu);
	
	file2 << "gamma\t omPrim\t omMin\t omMax\t probtot" << endl;
	double gamma = gammaMin;
	size_t jCount = 0;
	for (size_t jG=0;jG<nGamma;jG++) {
		gamma *= pasoGamma;
		double beta = sqrt(1.0-1.0/(gamma*gamma));
		double nuPrim = nuPrimMin;
		for (size_t jjNu=0;jjNu<nNuPrim;jjNu++) {
			nuPrim *= pasoNuPrim;
			double omPrim = planck*nuPrim/(electronMass*cLight2);
			double eps = omPrim/gamma;
			double nu = nuMin;
			
			double omMin = omPrim * (1.0-beta)/(1.0+2.0*beta*omPrim/gamma);
			double omMaxAbs = omPrim + (gamma-1.0);
			double omMax = omPrim*(1.0+beta)/(1.0-beta+2.0*omPrim/gamma);
			double aux = beta/(1.0+gamma*(1.0+beta));
			omMax = (eps < aux) ? omMax : omMaxAbs;
			
			struct two_d_params prob_params = {omPrim,gamma};
			gsl_function gsl_prob;
				gsl_prob.function = &probInt;	gsl_prob.params = &prob_params;
			double error;	int status;
			double probtot = RungeKuttaSimple(omMin,omMax,[&](double om)
						{return probexactNew(om,omPrim,gamma);});
			file2 << (float)gamma << "\t" << (float)nuPrim << "\t" <<
					omMin << "\t" << omMax << "\t" << probtot << endl;
			for (size_t jNu=0;jNu<nNu;jNu++) {
				nu *= pasoNu;
				double om = planck*nu / (electronMass*cLight2);
				double prob = 0.0;
				if (om > omMin && om < omMax && probtot > 0.9 && probtot < 1.2) {
					prob = probexactNew(om,omPrim,gamma)*planck/(electronMass*cLight2)/probtot;
				}
				file1 << prob << "\t";
				comptonprobvector[jCount++] = (float) prob;
			}
		}
	}
	fwrite(comptonprobvector,sizeof(float),jCount,fileBin);
	file1.close();
	file2.close();
	fclose(fileBin);
}

void comptonMatrix2()
{
	ofstream file1;	file1.open("comptonProbMatrix2.dat",ios::out);
	
	size_t nTemp, nNuPrim, nNu, nG;
	nTemp = nNuPrim = nNu = 30;
	nG = 5;
	
	double tempMin = 0.1;
	double tempMax = 10.0;
	double pasoTemp = pow(tempMax/tempMin,1.0/nTemp);
	double nuPrimMin = 1.0e8;
	double nuPrimMax = 1.0e20;
	double pasoNuPrim = pow(nuPrimMax/nuPrimMin,1.0/nNuPrim);
	double nuMin = 1.0e10;
	double nuMax = 1.0e22;
	double pasoNu = pow(nuMax/nuMin,1.0/nNu);
	
	void gaulag(double *x, double *w, int n, double alf);
	double *abscissas,*weights;
	abscissas=dvector(1,nG);
	weights=dvector(1,nG);
	double alf=0.0;
	gaulag(abscissas, weights,nG,alf);
	
	double temp = tempMin;
	for (size_t jTemp=0;jTemp<nTemp;jTemp++) {
		temp *= pasoTemp;
			
		double nuPrim = nuPrimMin;
		for (size_t jjNu=0;jjNu<nNuPrim;jjNu++) {
			nuPrim *= pasoNuPrim;
			double omPrim = planck*nuPrim/(electronMass*cLight2);
			
			double nu = nuMin;
			for (size_t jNu=0;jNu<nNu;jNu++) {
				nu *= pasoNu;
				double om = planck*nu / (electronMass*cLight2);
				double sum1 = 0.0;
				double sum2 = 0.0;
				for (size_t jG=1;jG<=nG;jG++) {
					double gamma=abscissas[jG]*temp+1.0;
					double cWeights = weights[jG]*gamma*sqrt(gamma*gamma-1.0);
					double beta = sqrt(1.0-1.0/(gamma*gamma));
					
					double eps = omPrim/gamma;
			
					double omMin = omPrim * (1.0-beta)/(1.0+2.0*beta*omPrim/gamma);
					double omMaxAbs = omPrim + (gamma-1.0);
					double omMax = omPrim*(1.0+beta)/(1.0-beta+2.0*omPrim/gamma);
					double aux = beta/(1.0+gamma*(1.0+beta));
					omMax = (eps < aux) ? omMax : omMaxAbs;
			
					struct two_d_params prob_params = {omPrim,gamma};
					gsl_function gsl_prob;
						gsl_prob.function = &probInt;	gsl_prob.params = &prob_params;
					double error;	int status;
					double probtot = RungeKuttaSimple(omMin,omMax,[&](double om)
							{return probexactNew(om,omPrim,gamma);});
					double prob = 0.0;
					if (om > omMin && om < omMax && probtot > 0.8 && probtot < 1.2) {
						prob = probexactNew(om,omPrim,gamma)*planck/(electronMass*cLight2)
									/probtot;
					}
					sum1 += prob*cWeights;
					sum2 += cWeights;
				}
				double result = sum1/sum2;
				file1 << result << endl;
				cout << jTemp << "\t" << jjNu << "\t" << jNu << endl;
			}
		}
	}
	file1.close();
}

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
		for (size_t jjE=0;jjE<nE;jjE++) {
			if (lum[jjE][jR]*(energies[jjE]/planck) > 1.0e10) {
				double omPrim = energies[jjE]/(electronMass*cLight2);
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
				if (probtotexact > 0.0) {// && probtotexact < 1.2) {
					double var_int_om=pow(omMax/omMin,1.0/nOm);
					double om=omMin;
					for (size_t jOm=0;jOm<nOm;jOm++) {
						double dOm = om*(var_int_om-1.0);
						double prob = probexactNew(om,omPrim,gamma)*dOm;
						p[((jR*nG+(jG-1))*nE+jjE)*nOm+jOm] = prob/probtotexact;
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
			double lumOutAux = lumIn * double(probab[((jR*nG+jG)*nE+jE)*nOm+jOm])*(om/omPrim);
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

double comptonNewNew(Vector comp, Vector lumIn, double normtemp, double nu, Vector energies,
					size_t nE)
{
	double nuPrimMin = 1.0e8;
	double nuPrimMax = 1.0e20;
	double pasoNuPrim = pow(energies[nE-1]/energies[0],1.0/nE);
	double sum = 0.0;
	for (size_t jjNu=1;jjNu<nE;jjNu++) {
		double nuPrim = energies[jjNu]/planck;
		if (nuPrim > nuPrimMin && nuPrim < nuPrimMax) {
			double dnuPrim = nuPrim * (pasoNuPrim-1.0);
			sum += lumIn[jjNu]*(nu/nuPrim)*probTemp(comp,nu,nuPrim,normtemp)*dnuPrim;
		}
	}
	return sum;
}

double comptonNewLocal(Vector comptonVec, Vector lum, double normtemp, size_t jE,
								Vector energies, size_t nE)
{
	double alf=0.0;
	int nG=10;
    double *abscissas,*weights;
    abscissas=dvector(1,nG);
    weights=dvector(1,nG);
    gaulag(abscissas,weights,nG,alf);
	
	double nu = energies[jE]/planck;
	double om = energies[jE]/(electronMass*cLight2);
	double nuPrimMin = 1.0e8;
	double nuPrimMax = 1.0e22;
	double pasoNuPrim = pow(energies[nE-1]/energies[0],1.0/nE);
	
	double sum = 0.0;
	double sum2=0.0;
	for (size_t jG=1;jG<=nG;jG++) {
		double gamma = abscissas[jG]*normtemp+1.0;
        double cWeights = weights[jG]*gamma*sqrt(gamma*gamma-1.0);
		double rOm = rate(om,gamma);
		double sum1 = 0.0;
		for (size_t jjNu=1;jjNu<nE;jjNu++) {
			double nuPrim = energies[jjNu]/planck;
			double omPrim = energies[jjNu]/(electronMass*cLight2);
			double rOmPrim = rate(omPrim,gamma);
			if (nuPrim > nuPrimMin && nuPrim < nuPrimMax && rOmPrim > 0.0) {
				double dnuPrim = nuPrim * (pasoNuPrim-1.0);
				sum1 += lum[jjNu] * (nu/nuPrim) * rOmPrim *
						probInterpolated(comptonVec,nu,nuPrim,gamma) * dnuPrim;
			}
		}
		sum += sum1*cWeights;
		if (rOm > 0.0) sum2 += cWeights*rOm;
	}
	return sum/sum2;
}