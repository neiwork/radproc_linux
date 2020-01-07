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
#include <fparameters/SpaceIterator.h>
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
			file2 << "gamma = " << (float)log10(gamma) << "\t nuPrim = " << (float)log10(nuPrim)
					<< "\t nuMin = " << log10(omMin*electronMass*cLight2/planck)
					<< "\t nuMax = " << log10(omMax*electronMass*cLight2/planck)
					<< "\t prob = " << probtot << endl;
			double nu = nuMinCompton;
			for (size_t jNu=0;jNu<nNuCompton;jNu++) {
				nu *= pasoNu;
				double om = planck*nu / (electronMass*cLight2);
				double prob = 0.0;
				if (om > omMin && om < omMax && probtot < 1.1) {
					prob = probexactNew(om,omPrim,gamma)*planck/(electronMass*cLight2);
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

void comptonNewNewPruebaVector(Vector tempVec, Vector nuPrimVec, Vector nuVec, Vector probVec,
			size_t jTemp, Vector energies, size_t jE, Vector& p, double normtemp)
{
	double nu = energies[jE]/planck;
	for (size_t jjNu=0;jjNu<nE;jjNu++) {
		double nuPrim = energies[jjNu]/planck;
		if (comptonMethod == 0)
			p[(jTemp*nE+jjNu)*nE+jE] = probInterpolated2(tempVec,nuPrimVec,nuVec,probVec,
											nu,nuPrim,normtemp);
		else if (comptonMethod == 1)
			p[(jTemp*nE+jjNu)*nE+jE] = probTemp(probVec,nu,nuPrim,normtemp);
	}
}

/*	FILE *emissFile;
	emissFile=fopen("emissivities2.bin","rb");
	float *emiss;
	emiss=calloc(numE*numR*numTheta,sizeof(float));           // Emissivity vector
	fread(emiss,sizeof(float),numE*numR*numTheta,emissFile);
	fclose(emissFile); */

void cNew(State& st, Vector& p, Vector redshift)
{
	FILE *probsFile, *ntFile, *omFile, *ompFile;
	probsFile = fopen("probs.bin","rb");
	ntFile = fopen("nt.bin","rb");
	omFile = fopen("om.bin","rb");
	ompFile = fopen("omp.bin","rb");
	
	float *probVec, *ntVec, *omVec, *ompVec;
	probVec = (float*)calloc(nTempCompton*nNuCompton*nNuPrimCompton,sizeof(float));  
	ntVec = (float*)calloc(nTempCompton,sizeof(float)); 
	omVec = (float*)calloc(nNuCompton,sizeof(float));
	ompVec = (float*)calloc(nNuPrimCompton,sizeof(float));
	fread(probVec,sizeof(float),nTempCompton*nNuCompton*nNuPrimCompton,probsFile);
	fread(ntVec,sizeof(float),nTempCompton,ntFile);
	fread(omVec,sizeof(float),nNuCompton,omFile);
	fread(ompVec,sizeof(float),nNuPrimCompton,ompFile);
	fclose(probsFile);
	fclose(ntFile);
	fclose(omFile);
	fclose(ompFile);
	
	size_t jR = 0;
	st.photon.ps.iterate([&](const SpaceIterator& iR) {
		double lognt = log10(boltzmann*st.tempElectrons.get(iR)/electronRestEnergy);
		double aux = ntVec[0];
		size_t pos_r = 0;
		while (aux < lognt) {
			pos_r++;
			aux = ntVec[pos_r];
		}
		double lognt1 = ntVec[pos_r-1];
		double lognt2 = ntVec[pos_r];
		size_t jE = 0;
		st.photon.ps.iterate([&](const SpaceIterator& iRE) {
			double logom = log10(iRE.val(DIM_E)/redshift[jR]/electronRestEnergy);
			if (logom > omVec[0] && logom < omVec[nNuCompton-1]) {
				double aux = omVec[0];
				size_t pos_om = 0;
				while (aux < logom) {
					pos_om++;
					aux = omVec[pos_om];
				}
				double logom1 = omVec[pos_om-1];
				double logom2 = omVec[pos_om];
				size_t jjE = 0;
				st.photon.ps.iterate([&](const SpaceIterator& iREE) {
					double logomp = log10(iREE.val(DIM_E)/redshift[jR]/electronRestEnergy);
					if (logomp > ompVec[0] && logomp < ompVec[nNuCompton-1]) {
						double aux = ompVec[0];
						size_t pos_omp = 0;
						while (aux < logomp) {
							pos_omp++;
							aux = ompVec[pos_omp];
						}
						double logomp1 = ompVec[pos_omp-1];
						double logomp2 = ompVec[pos_omp];
						
						double prob111 = probVec[((pos_r-1)*nNuPrimCompton+(pos_omp-1))*nNuCompton+(pos_om-1)];
						double prob112 = probVec[((pos_r-1)*nNuPrimCompton+(pos_omp-1))*nNuCompton+pos_om];
						double prob121 = probVec[((pos_r-1)*nNuPrimCompton+pos_omp)*nNuCompton+(pos_om-1)];
						double prob122 = probVec[((pos_r-1)*nNuPrimCompton+pos_omp)*nNuCompton+pos_om];
						double prob211 = probVec[(pos_r*nNuPrimCompton+(pos_omp-1))*nNuCompton+(pos_om-1)];
						double prob212 = probVec[(pos_r*nNuPrimCompton+(pos_omp-1))*nNuCompton+pos_om];
						double prob221 = probVec[(pos_r*nNuPrimCompton+pos_omp)*nNuCompton+(pos_om-1)];
						double prob222 = probVec[(pos_r*nNuPrimCompton+pos_omp)*nNuCompton+pos_om];
						
						double prob11 = (prob112-prob111)/(logom2-logom1) * (logom-logom1) + prob111;
						double prob12 = (prob122-prob121)/(logom2-logom1)*(logom-logom1) + prob121;
						double prob1 = (prob12-prob11)/(logomp2-logomp1)*(logomp-logomp1)+prob11;
						
						double prob21 = (prob212-prob211)/(logom2-logom1) * (logom-logom1) + prob211;
						double prob22 = (prob222-prob221)/(logom2-logom1) * (logom-logom1) + prob221;
						double prob2 = (prob22-prob21)/(logomp2-logomp1)*(logomp-logomp1)+prob21;
						
						double prob = (prob2-prob1)/(lognt2-lognt1)*(lognt-lognt1) + prob1;
						p[(jR*nE+jjE)*nE+jE] = pow(10.0,prob);
					} else
						p[(jR*nE+jjE)*nE+jE] = 0.0;
					jjE++;
				},{-1,0,0});
			} else
				for (size_t jjE=0;jjE<nE;jjE++)
					p[(jR*nE+jjE)*nE+jE] = 0.0;
			jE++;
		},{-1,0,0});
		jR++;
	},{0,-1,0});
}

void comptonNewNewNewPruebaVector(size_t jR, Vector energies, size_t jE, Vector& p, double normtemp)
{
	gsl_function gsl_extrinf;
	gsl_extrinf.function = &extrinf;
	double alf=0.0;
	int nG=8;
	double pasoOm = pow(energies[nE-1]/energies[0],1.0/(nE+1));
    double *abscissas,*weights;
    abscissas=dvector(1,nG);
    weights=dvector(1,nG);
    gaulag(abscissas,weights,nG,alf);
	double om = energies[jE]/(electronMass*cLight2);
	double om1 = om/sqrt(pasoOm);
	double om2 = om*sqrt(pasoOm);
	for (size_t jjNu=0;jjNu<nE;jjNu++) {
		double omPrim = energies[jjNu]/(electronMass*cLight2);
		
		double sum1 = 0.0;
		double sum2 = 0.0;
		for (int jG=1;jG<=nG;jG++) {
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
			
			//double probtotal = RungeKuttaSimple(omMin,omMax,[&](double om1)
			//		{return probexactNew(om1,omPrim,gamma);});

			//double prob1 = (omMin < om1 && om1 < omMax) ? // && probtotal < 1.2) ?
			//				probexactNew(om1,omPrim,gamma) : 0.0;
			//double prob2 = (omMin < om2 && om2 < omMax) ? // && probtotal < 1.2) ?
			//				probexactNew(om2,omPrim,gamma) : 0.0;
			double prob = (omMin < om && om < omMax) ? // && probtotal < 1.2) ?
							//probexactNew(om,omPrim,gamma) : 0.0;
							probexactR(om,omPrim,gamma) : 0.0;
			sum1 += prob*cWeights;
			//sum1 += (prob2+prob1)/2.0 * cWeights;
			sum2 += cWeights;
		}
		p[(jR*nE+jjNu)*nE+jE] = (sum1/sum2) * (planck/(electronMass*cLight2));
	}
}

double compton(Vector p, Vector lumIn, size_t jTemp, Vector energies, size_t jE)
{
	double pasoNuPrim = pow(energies[nE-1]/energies[0],1.0/nE);
	double sum = 0.0;
	for (size_t jjNu=0;jjNu<nE;jjNu++) {
		double nuPrim = energies[jjNu]/planck;
		if (nuPrim > nuPrimMinCompton && nuPrim < nuPrimMaxCompton)
		{
			double dNuPrim = nuPrim * (sqrt(pasoNuPrim)-1.0/sqrt(pasoNuPrim));
			double prob = p[(jTemp*nE+jjNu)*nE+jE];
			sum += lumIn[jjNu]*dNuPrim*prob*(energies[jE]/energies[jjNu]);
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
		for (int jG=1;jG<=nG;jG++) {
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

			double prob1 = (omMin < om && om < omMax) ? // && probtotal < 1.2) ?
							probexactNew(om,omPrim,gamma) : 0.0;
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

double numberPhotons(Vector nPh, double nu, Vector energies)
{
	double lognu_actual = log(nu);
	double aux = log(energies[0]/planck);
	size_t pos_nu = 0;
	while (aux < lognu_actual) {
		pos_nu++;
		aux = log(energies[pos_nu]/planck);
	}
	if (nPh[pos_nu] == 0 || nPh[pos_nu-1] == 0) return 0.0;
	else {
		double m = log(nPh[pos_nu]/nPh[pos_nu-1])/log(energies[pos_nu]/energies[pos_nu-1]);
		double log_nPh = m*(log(nu*planck/energies[pos_nu-1])) + log(nPh[pos_nu-1]);
		return exp(log_nPh);
	}
}

double comptonNewLocal2(Vector nPh, double normtemp, size_t jE, Vector energies)
{
	
	gsl_function gsl_extrinf;
	gsl_extrinf.function = &extrinf;
	
	double alf=0.0;
	int nG=8;
    double *abscissas,*weights;
    abscissas=dvector(1,nG);
    weights=dvector(1,nG);
    gaulag(abscissas,weights,nG,alf);
	
	double om = energies[jE]/electronRestEnergy;
	double nu = energies[jE]/planck;
	
	double sum = 0.0;
	double sum2 = 0.0;
	
	for (int jG=1;jG<=nG;jG++) {
		double gamma = abscissas[jG]*normtemp+1.0;
		double beta = sqrt(1.0-1.0/(gamma*gamma));
        double cWeights = weights[jG]*gamma*sqrt(gamma*gamma-1.0);
		
		
		double factor = 4.0/3.0 * gamma*gamma;
		double omPrim = om/factor;
		double rOmPrim = (om < omPrim + gamma) ? rate(omPrim,gamma) : 0.0;
		double nuPrim = omPrim * electronRestEnergy / planck;
		if (nuPrim > nuPrimMinCompton && nuPrim < nuPrimMaxCompton && rOmPrim > 0.0) {
			double nPhotons = numberPhotons(nPh,nuPrim,energies);
			nPhotons * (1.0+2.0*pi*numberPhotons(nPh,nu,energies)*cLight2*cLight/(nu*nu));
			sum += nPhotons * rOmPrim * cWeights / factor;
		}
		
		/*
		double pasoOmPrim = pow(energies[nE-1]/energies[0],1.0/nE);
		double sum1 = 0.0;
		for (size_t jjNu=0;jjNu<nE;jjNu++) {
			double nuPrim = energies[jjNu]/planck;
			double omPrim = energies[jjNu]/electronRestEnergy;
			double eps = omPrim/gamma;
			struct two_d_params extrinf_params = {eps,gamma};
			gsl_extrinf.params = &extrinf_params;
			int status1,status2;
		
			double omMin = omPrim * brent(&gsl_extrinf,0.0,1.0,&status1,&status2);
			double omMaxAbs = omPrim + (gamma - 1.0);
			double omMax = omPrim*(1.0+beta)/(1.0-beta+2.0*omPrim/gamma);
			double aux = beta/(1.0+gamma*(1.0+beta));
			omMax = (eps < aux) ? omMax : omMaxAbs;
			double rOmPrim = rate(omPrim,gamma);
			if (nuPrim > nuPrimMinCompton && nuPrim < nuPrimMaxCompton && rOmPrim > 0.0) {
				double dOmPrim = omPrim * (pasoOmPrim-1.0);
				double prob = (om < omMax && om > omMin) ? probexactNew(om,omPrim,gamma) : 0.0;
				sum1 += nPh[jjNu] * rOmPrim * prob * dOmPrim;
			}
		}
		sum += sum1*cWeights;
		*/
		
		sum2 += cWeights;
	}
	
	return sum/sum2;
}