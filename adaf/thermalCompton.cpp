#include "thermalCompton.h"
#include "globalVariables.h"
#include <fluminosities/probexact.h>
#include <nrMath/brent.h>
#include <nrMath/integrators.h>
#include <fmath/physics.h>
#include <fmath/RungeKutta.h>
#include <fmath/interpolation.h>
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

double probInterpolated(Vector comp, double nu, double nuPrim, double gamma)
{
	double pasoGamma = pow(gammaMaxCompton/gammaMinCompton,1.0/nGammaCompton);
	double pasoNuPrim = pow(nuPrimMaxCompton/nuPrimMinCompton,1.0/nNuPrimCompton);
	double pasoNu = pow(nuMaxCompton/nuMinCompton,1.0/nNuCompton);

	size_t kGamma, kNuPrim, kNu;
	kGamma = kNuPrim = kNu = 0;
	double gammaAux = gammaMinCompton;
	while (gammaAux < gamma) {
		gammaAux *= pasoGamma;
		kGamma++;
	}
	double nuPrimAux = nuPrimMinCompton;
	while (nuPrimAux < nuPrim) {
		nuPrimAux *= pasoNuPrim;
		kNuPrim++;
	}
	double nuAux = nuMinCompton;
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
	
	double prob111 = comp[((kGamma-1)*nNuPrimCompton+(kNuPrim-1))*nNuCompton+(kNu-1)];
	double prob112 = comp[((kGamma-1)*nNuPrimCompton+(kNuPrim-1))*nNuCompton+kNu];
	double prob121 = comp[((kGamma-1)*nNuPrimCompton+kNuPrim)*nNuCompton+(kNu-1)];
	double prob122 = comp[((kGamma-1)*nNuPrimCompton+kNuPrim)*nNuCompton+kNu];
	double prob211 = comp[(kGamma*nNuPrimCompton+(kNuPrim-1))*nNuCompton+(kNu-1)];
	double prob212 = comp[(kGamma*nNuPrimCompton+(kNuPrim-1))*nNuCompton+kNu];
	double prob221 = comp[(kGamma*nNuPrimCompton+kNuPrim)*nNuCompton+(kNu-1)];
	double prob222 = comp[(kGamma*nNuPrimCompton+kNuPrim)*nNuCompton+kNu];
	
	double prob11 = (prob112-prob111)/(logNu2-logNu1) * (logNu-logNu1) + prob111;
	double prob12 = (prob122-prob121)/(logNu2-logNu1)*(logNu-logNu1) + prob121;
	double prob1 = (prob12-prob11)/(logNuPrim2-logNuPrim1)*(logNuPrim-logNuPrim1)+prob11;
	
	double prob21 = (prob212-prob211)/(logNu2-logNu1) * (logNu-logNu1) + prob211;
	double prob22 = (prob222-prob221)/(logNu2-logNu1)*(logNu-logNu1) + prob221;
	double prob2 = (prob22-prob21)/(logNuPrim2-logNuPrim1)*(logNuPrim-logNuPrim1)+prob21;
	
	double prob = (prob2-prob1)/(logGamma2-logGamma1)*(logGamma-logGamma1) + prob1;
	return prob;
}

double probTemp(Vector comp, double nu, double nuPrim, double normtemp)
{
	void gaulag(double *x, double *w, int n, double alf);
	size_t nG = 6;
	double *abscissas,*weights;
	abscissas=dvector(1,nG);
	weights=dvector(1,nG);
	double alf=0.0;
	gaulag(abscissas, weights,nG,alf);
	double sum1=0.0;
	double sum2=0.0;
	
	for (size_t jG=1;jG<=nG;jG++) {
		double gamma=abscissas[jG]*normtemp+1.0;
		double cWeights = weights[jG]*gamma*sqrt(gamma*gamma-1.0);
		sum1 += probInterpolated(comp,nu,nuPrim,gamma)*cWeights;
		sum2 += cWeights;
	}
	
	return sum1/sum2;
}


double probInterpolated2(Vector tempVec, Vector nuPrimVec, Vector nuVec, Vector probVec,
							double nu, double nuPrim, double temp)
{
	size_t kTemp, kNuPrim, kNu;
	double logTemp = log(temp); double logNuPrim = log(nuPrim); double logNu = log(nu);
	
	/*
	locate(tempVec,nTempCompton,logTemp,kTemp);
	locate(nuPrimVec,nNuPrimCompton,logNuPrim,kNuPrim);
	locate(nuVec,nNuCompton,logNu,kNu);
	
	if (kTemp >= nTempCompton || kTemp < 0) return 0.0;
	if (kNuPrim >= nNuPrimCompton || kNuPrim < 0) return 0.0;
	if (kNu >= nNuCompton || kNu < 0) return 0.0;*/
	
	////////////////////////////////////

	kTemp = kNuPrim = kNu = 0;
	double pasoTemp = pow(tempMaxCompton/tempMinCompton,1.0/nTempCompton);
	double pasoNuPrim = pow(nuPrimMaxCompton/nuPrimMinCompton,1.0/nNuPrimCompton);
	double pasoNu = pow(nuMaxCompton/nuMinCompton,1.0/nNuCompton);
	
	double tempAux = tempMinCompton;
	while (tempAux < temp) {
		tempAux *= pasoTemp;
		kTemp++;
	}
	double nuPrimAux = nuPrimMinCompton;
	while (nuPrimAux < nuPrim) {
		nuPrimAux *= pasoNuPrim;
		kNuPrim++;
	}
	double nuAux = nuMinCompton;
	while (nuAux < nu) {
		nuAux *= pasoNu;
		kNu++;
	}
	
	kTemp--; kNuPrim--; kNu--;
	
	double logTemp1 = log(tempAux/pasoTemp);
	double logTemp2 = log(tempAux);
	double logNuPrim1 = log(nuPrimAux/pasoNuPrim);
	double logNuPrim2 = log(nuPrimAux);
	double logNu1 = log(nuAux/pasoNu);
	double logNu2 = log(nuAux);

	//////////////////////////////////
	/*
	double logTemp1 = tempVec[kTemp];
	double logTemp2 = tempVec[kTemp+1];
	double logNuPrim1 = nuPrimVec[kNuPrim];
	double logNuPrim2 = nuPrimVec[kNuPrim+1];
	double logNu1 = nuVec[kNu];
	double logNu2 = nuVec[kNu+1];*/
	
	double prob111 = probVec[(kTemp*nNuPrimCompton+kNuPrim)*nNuCompton+kNu];
	double prob112 = probVec[(kTemp*nNuPrimCompton+kNuPrim)*nNuCompton+(kNu+1)];
	double prob121 = probVec[(kTemp*nNuPrimCompton+(kNuPrim+1))*nNuCompton+kNu];
	double prob122 = probVec[(kTemp*nNuPrimCompton+(kNuPrim+1))*nNuCompton+(kNu+1)];
	double prob211 = probVec[((kTemp+1)*nNuPrimCompton+kNuPrim)*nNuCompton+kNu];
	double prob212 = probVec[((kTemp+1)*nNuPrimCompton+kNuPrim)*nNuCompton+(kNu+1)];
	double prob221 = probVec[((kTemp+1)*nNuPrimCompton+(kNuPrim+1))*nNuCompton+kNu];
	double prob222 = probVec[((kTemp+1)*nNuPrimCompton+(kNuPrim+1))*nNuCompton+(kNu+1)];
	
	double prob11 = (prob112-prob111)/(logNu2-logNu1) * (logNu-logNu1) + prob111;
	double prob12 = (prob122-prob121)/(logNu2-logNu1)*(logNu-logNu1) + prob121;
	double prob1 = (prob12-prob11)/(logNuPrim2-logNuPrim1)*(logNuPrim-logNuPrim1)+prob11;
	
	double prob21 = (prob212-prob211)/(logNu2-logNu1) * (logNu-logNu1) + prob211;
	double prob22 = (prob222-prob221)/(logNu2-logNu1)*(logNu-logNu1) + prob221;
	double prob2 = (prob22-prob21)/(logNuPrim2-logNuPrim1)*(logNuPrim-logNuPrim1)+prob21;
	
	double prob = (prob2-prob1)/(logTemp2-logTemp1)*(logTemp-logTemp1) + prob1;
	return prob;
}

void comptonMatrix()
{
	ofstream file1;	file1.open("comptonProbMatrix.dat",ios::out);
	ofstream file2;	file2.open("comptonProbTot.dat",ios::out);

	double pasoGamma = pow(gammaMaxCompton/gammaMinCompton,1.0/nGammaCompton);
	double pasoNuPrim = pow(nuPrimMaxCompton/nuPrimMinCompton,1.0/nNuPrimCompton);
	double pasoNu = pow(nuMaxCompton/nuMinCompton,1.0/nNuCompton);
	
	file2 << "gamma\t omPrim\t omMin\t omMax\t probtot" << endl;
	gsl_function gsl_extrinf;
	gsl_extrinf.function = &extrinf;
	double gamma = gammaMinCompton;
	for (size_t jG=0;jG<nGammaCompton;jG++) {
		gamma *= pasoGamma;
		double beta = sqrt(1.0-1.0/(gamma*gamma));
		double nuPrim = nuPrimMinCompton;
		for (size_t jjNu=0;jjNu<nNuPrimCompton;jjNu++) {
			nuPrim *= pasoNuPrim;
			double omPrim = planck*nuPrim/(electronMass*cLight2);
			double eps = omPrim/gamma;
			
			double omMin = omPrim * (1.0-beta)/(1.0+2.0*beta*omPrim/gamma);
			struct two_d_params extrinf_params = {eps,gamma};
			gsl_extrinf.params = &extrinf_params;
			int status1,status2;
		
			omMin = omPrim * brent(&gsl_extrinf,0.0,1.0,&status1,&status2);
			
			
			double omMaxAbs = omPrim + (gamma-1.0);
			double omMax = omPrim*(1.0+beta)/(1.0-beta+2.0*omPrim/gamma);
			double aux = beta/(1.0+gamma*(1.0+beta));
			omMax = (eps < aux) ? omMax : omMaxAbs;

			double probtot = RungeKuttaSimple(omMin,omMax,[&](double om)
						{return probexactNew(om,omPrim,gamma);});
			file2 << (float)gamma << "\t" << (float)omPrim << "\t" <<
					omMin << "\t" << omMax << "\t" << probtot << endl;
			double nu = nuMinCompton;
			for (size_t jNu=0;jNu<nNuCompton;jNu++) {
				nu *= pasoNu;
				double om = planck*nu / (electronMass*cLight2);
				double prob = 0.0;
				if (om > omMin && om < omMax) {// && probtot > 0.9 && probtot < 1.1) {
					prob = probexactNew(om,omPrim,gamma)*planck/(electronMass*cLight2)/probtot;
				}
				file1 << prob << endl;
			}
		}
	}
	file1.close();
	file2.close();
}

void comptonMatrix2()
{
	double pasoTemp = pow(tempMaxCompton/tempMinCompton,1.0/nTempCompton);
	double pasoNuPrim = pow(nuPrimMaxCompton/nuPrimMinCompton,1.0/nNuPrimCompton);
	double pasoNu = pow(nuMaxCompton/nuMinCompton,1.0/nNuCompton);

	void gaulag(double *x, double *w, int n, double alf);	size_t nG = 8;
	double *abscissas,*weights;
	abscissas=dvector(1,nG);
	weights=dvector(1,nG);
	double alf=0.0;
	gaulag(abscissas,weights,nG,alf);
	
	ofstream fileTempVec,fileNuPrimVec,fileNuVec,fileProbs,fileSizes;
	
	Vector temp(nTempCompton,tempMinCompton*pasoTemp);
	Vector nuPrim(nNuPrimCompton,nuPrimMinCompton*pasoNuPrim);
	Vector nu(nNuCompton,nuMinCompton*pasoNu);
	Vector probVec(nTempCompton*nNuPrimCompton*nNuCompton,0.0);
	
	fileTempVec.open("tempComptonVec.dat",ios::out);
	fileTempVec << log(temp[0]) << endl;
	for (size_t jTemp=1;jTemp<nTempCompton;jTemp++) {
		temp[jTemp] = temp[jTemp-1]*pasoTemp;
		fileTempVec << log(temp[jTemp]) << endl;
	}
	fileTempVec.close();
	
	fileNuPrimVec.open("nuPrimComptonVec.dat",ios::out);
	fileNuPrimVec << log(nuPrim[0]) << endl;
	for (size_t jjNu=1;jjNu<nNuPrimCompton;jjNu++) {
		nuPrim[jjNu] = nuPrim[jjNu-1]*pasoNuPrim;
		fileNuPrimVec << log(nuPrim[jjNu]) << endl;
	}
	fileNuPrimVec.close();
	
	fileNuVec.open("nuComptonVec.dat",ios::out);
	fileNuVec << log(nu[0]) << endl;
	for (size_t jNu=1;jNu<nNuCompton;jNu++) {
		nu[jNu] = nu[jNu-1]*pasoNu;
		fileNuVec << log(nu[jNu]) << endl;
	}
	fileNuVec.close();
	
	#pragma omp parallel for
	for (size_t jTemp=0;jTemp<nTempCompton;jTemp++) {
		for (size_t jjNu=0;jjNu<nNuPrimCompton;jjNu++) {
			
			double omPrim = planck*nuPrim[jjNu]/(electronMass*cLight2);
			for (size_t jNu=0;jNu<nNuCompton;jNu++) {
				double om = planck*nu[jNu] / (electronMass*cLight2);
				double sum1 = 0.0;
				double sum2 = 0.0;
				for (size_t jG=1;jG<=nG;jG++) {
					double gamma=abscissas[jG]*temp[jTemp]+1.0;
					double cWeights = weights[jG]*gamma*sqrt(gamma*gamma-1.0);
					double beta = sqrt(1.0-1.0/(gamma*gamma));
					
					double eps = omPrim/gamma;
			
					double omMin = omPrim * (1.0-beta)/(1.0+2.0*beta*omPrim/gamma);
					double omMaxAbs = omPrim + (gamma-1.0);
					double omMax = omPrim*(1.0+beta)/(1.0-beta+2.0*omPrim/gamma);
					double aux = beta/(1.0+gamma*(1.0+beta));
					omMax = (eps < aux) ? omMax : omMaxAbs;

					double probtot = RungeKuttaSimple(omMin,omMax,[&](double om)
							{return probexactNew(om,omPrim,gamma);});
					double prob = 0.0;
					if (om > omMin && om < omMax) {// && probtot > 0.8 && probtot < 1.2) {
						prob = probexactNew(om,omPrim,gamma)*planck/(electronMass*cLight2)
									/probtot;
					}
					sum1 += prob*cWeights;
					sum2 += cWeights;
				}
				probVec[(jTemp*nNuPrimCompton+jjNu)*nNuCompton+jNu] = sum1/sum2;
			}
		}
		cout << "jTemp = " << jTemp << endl;
	}
	fileProbs.open("comptonProbMatrix2.dat",ios::out);
	for (size_t k=0;k<nTempCompton*nNuPrimCompton*nNuCompton;k++)
		fileProbs << probVec[k] << endl;

	fileProbs.close();
}

double comptonNewNew(Vector tempVec, Vector nuPrimVec, Vector nuVec, Vector probVec,
					Vector lumIn, double normtemp, double nu, Vector energies, size_t jE)
{
	double pasoNuPrim = pow(energies[nE-1]/energies[0],1.0/nE);
	double sum = 0.0;
	for (size_t jjNu=0;jjNu<nE;jjNu++) {
		double nuPrim = energies[jjNu]/planck;
		if (nuPrim > nuPrimMinCompton && nuPrim < nuPrimMaxCompton && lumIn[jjNu] > 1.0e8)
		{
			double dNuPrim = nuPrim * (pasoNuPrim-1.0);
			double prob = 0.0;
			if (comptonMethod == 0)
				prob = probInterpolated2(tempVec,nuPrimVec,nuVec,probVec,
													nu,nuPrim,normtemp);
			else if (comptonMethod == 1)
				prob = probTemp(probVec,nu,nuPrim,normtemp);

			sum += lumIn[jjNu]*dNuPrim*prob*(nu/nuPrim);
		}
	}
	return sum;
}

double comptonNewNewNew(Vector lumIn, double normtemp, double nu, Vector energies, size_t jE)
{
	gsl_function gsl_extrinf;
	gsl_extrinf.function = &extrinf;
	double alf=0.0;
	int nG=8;
    double *abscissas,*weights;
    abscissas=dvector(1,nG);
    weights=dvector(1,nG);
    gaulag(abscissas,weights,nG,alf);
	double om = energies[jE]/(electronMass*cLight2);
	double pasoNuPrim = pow(energies[nE-1]/energies[0],1.0/nE);
	double sum = 0.0;
	for (size_t jjE=0;jjE<nE;jjE++) {
		double nuPrim = energies[jjE]/planck;
		double dNuPrim = nuPrim * (pasoNuPrim-1.0);
		double omPrim = energies[jjE]/(electronMass*cLight2);
		
		double sum1 = 0.0;
		double sum2 = 0.0;
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
			
			double prob1 = (omMin < om && om < omMax) ? probexactNew(om,omPrim,gamma) : 0.0;
			sum1 += prob1*cWeights;
			sum2 += cWeights;
		}
		double prob = sum1/sum2 * planck /(electronMass*cLight2);
		sum += lumIn[jjE]*(nu/nuPrim)*dNuPrim*prob;
	}
	return sum;
}

double comptonNewLocal(Vector comptonVec, Vector lum, double normtemp, size_t jE,
								Vector energies, size_t nE)
{
	double alf=0.0;
	int nG=8;
    double *abscissas,*weights;
    abscissas=dvector(1,nG);
    weights=dvector(1,nG);
    gaulag(abscissas,weights,nG,alf);
	
	double nu = energies[jE]/planck;
	double pasoNuPrim = pow(energies[nE-1]/energies[0],1.0/nE);
	
	double sum = 0.0;
	double sum2 = 0.0;
	for (int jG=1;jG<=nG;jG++) {
		double gamma = abscissas[jG]*normtemp+1.0;
        double cWeights = weights[jG]*gamma*sqrt(gamma*gamma-1.0);
		double sum1 = 0.0;
		for (size_t jjNu=1;jjNu<nE;jjNu++) {
			double nuPrim = energies[jjNu]/planck;
			double omPrim = energies[jjNu]/(electronMass*cLight2);
			double rOmPrim = rate(omPrim,gamma);
			if (nuPrim > nuPrimMinCompton && nuPrim < nuPrimMaxCompton && rOmPrim > 0.0) {
				double dNuPrim = nuPrim * (pasoNuPrim-1.0);
				sum1 += lum[jjNu] * rOmPrim * (nu/nuPrim) *
						probInterpolated(comptonVec,nu,nuPrim,gamma) * dNuPrim;
			}
		}
		sum += (sum1*cWeights);
		sum2 += cWeights;
	}
	
	return sum/sum2;
}