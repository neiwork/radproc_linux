// Standard C++ libraries
#include <iomanip>

// Boost libraries
#include <boost/property_tree/ptree.hpp>
#include <boost/math/tools/roots.hpp>

// Standard user-made libraries
#include <fparameters/SpaceIterator.h>
#include <fparameters/parameters.h>
#include <fluminosities/thermalSync.h>
#include <fluminosities/thermalBremss.h>
#include <fluminosities/luminosityHadronic.h>
#include <fluminosities/luminositySynchrotron.h>
#include <fluminosities/blackBody.h>
#include <fluminosities/reflection.h>
#include <fluminosities/probexact.h>
#include <fluminosities/luminosityNTHadronic.h>
#include <fmath/mathFunctions.h>
#include <fmath/physics.h>
#include <fmath/fbisection.h>
#include <fparameters/Dimension.h>
#include "absorption.h"
#include <fmath/RungeKutta.h>
// Project headers
#include "globalVariables.h"
#include "messages.h"
#include "read.h"
#include "thermalCompton.h"
#include "thermalProcesses.h"
#include "adafFunctions.h"
#include "write.h"

// Namespaces
using namespace std;

void localProcesses(State& st, Matrix& lumOutSy, Matrix& lumOutBr, Matrix& lumOutpp,
						Vector& energies, const int flags[], Matrix& lumOut)
{
	size_t jE=0;
	st.photon.ps.iterate([&](const SpaceIterator& itE) {
		energies[jE++] = itE.val(DIM_E);
	},{-1,0,0});

	#pragma omp parallel for
	for (size_t jE=0;jE<nE;jE++) {
		double frequency=energies[jE]/planck;
		double lumSync1,lumSync2,lumBremss1,lumBremss2;
		lumSync1 = lumSync2 = lumBremss1 = lumBremss2 = 0.0;
		size_t jR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itER) {
			double r = itER.val(DIM_R);
			double thetaH = st.thetaH.get(itER);
			double rB1 = r/sqrt(paso_r);
			double rB2 = r*sqrt(paso_r);
			//double area = 4.0*pi*rB2*rB2*cos(thetaH);
			double area = 2.0*pi*( (height_method == 0) ? rB2*rB2*(2.0*cos(thetaH)+P2(sin(thetaH))) :
									2.0*rB2*height_fun(r)+(rB2*rB2-rB1*rB1) );
			double fluxToLum = area;
			double vol = volume(r);
			double emissToLum = vol*4.0*pi;
			double lumBr = 0.0;
			double lumSy = 0.0;
			double lumRJ = 0.0;

			double temp = st.tempElectrons.get(itER);
			double temp_i = st.tempIons.get(itER);
			double magf = st.magf.get(itER);
			double dens_i = st.denf_i.get(itER);
			double dens_e = st.denf_e.get(itER);
			double jSy = jSync(energies[jE],temp,magf,dens_e);
			if (flags[0]) {
				lumRJ = pi*bb(frequency,temp)*fluxToLum;
				lumSy = jSy*emissToLum;
			}
			if (flags[1])
				lumBr = jBremss(energies[jE],temp,dens_i,dens_e)*emissToLum;
			if (flags[2] && temp_i > 1.0e11 && frequency > 1.0e20 && frequency < 1.0e26) {
				double jpp = luminosityHadronic(energies[jE],dens_i,temp_i);
				lumOutpp[jE][jR] = jpp*vol;
			}
            if (flags[0]) {
				if (jR > 0)
					lumSync2 = min((1.0-scattAA[jR-1][jR])*lumSync1+lumSy,lumRJ);
				else 
					lumSync2 = min(lumSy,lumRJ);
				if (lumSync2 > lumSync1)
					lumOutSy[jE][jR] = lumSync2-lumSync1;
				lumSync1 = lumSync2;
            }
			
			if (flags[1]) {
				if (flags[0] && frequency < 1.0e14) {
					if (lumBr > lumRJ)
						lumBr = 0.0;
				}
				lumOutBr[jE][jR]=lumBr;
			}
			
			lumOut[jE][jR] = lumOutSy[jE][jR]+lumOutBr[jE][jR]+lumOutpp[jE][jR];
			jR++;
		},{jE,-1,0});
	}
}

void localProcesses2(State& st, Matrix& lumOutSy, Matrix& lumOutBr,
						Vector& energies, const int flags[], Matrix& lumOut)
{
	size_t jE=0;
	st.photon.ps.iterate([&](const SpaceIterator& itE) {
		energies[jE++] = itE.val(DIM_E);
	},{-1,0,0});

	#pragma omp parallel for
	for (int jE=0;jE<nE;jE++) {
		size_t jR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itER) {
			double r = itER.val(DIM_R);
			double rB2 = r*sqrt(paso_r);
			double rB1 = rB2/paso_r;
			double temp_e = st.tempElectrons.get(itER);
			double temp_i = st.tempIons.get(itER);
			double magf = st.magf.get(itER);
			double dens_i = st.denf_i.get(itER);
			double dens_e = st.denf_e.get(itER);
			double xSy,xBr;
			double unrEnergy = energies[jE]/redshift_to_inf[jR];
			unrEnergy = energies[jE];
			xSy = xBr = 0.0;
			if (flags[0])
				xSy = jSync(unrEnergy,temp_e,magf,dens_e)*4.0*pi;
			if (flags[1])
				xBr = jBremss(unrEnergy,temp_e,dens_i,dens_e)*4.0*pi;
			
			SpaceCoord psc = {jE,jR,0};
			double b_nu = bb(unrEnergy/planck,temp_e);
			double kappa = (xSy+xBr)/(4.0*pi*b_nu);
			double tau = 0.5*sqrt(pi)*kappa*height_fun(r);
			double fluxSy = (2.0*sqrt(3.0)*tau > 1.0e-9) ? 
							2.0*pi/sqrt(3.0)*b_nu*(1.0-exp(-2.0*sqrt(3.0)*tau)) :
							0.5*sqrt(pi)*xSy*height_fun(r);
			double fluxBr = (2.0*sqrt(3.0)*tau < 1.0e-3) ? 0.5*sqrt(pi)*xBr*height_fun(r) : 0.0;
							
			lumOutSy[jE][jR] = 2.0 * pi*(rB2*rB2-rB1*rB1) * fluxSy;
			lumOutBr[jE][jR] = 2.0 * pi*(rB2*rB2-rB1*rB1) * fluxBr;
			lumOut[jE][jR] = lumOutSy[jE][jR] + lumOutBr[jE][jR];
			jR++;
		},{jE,-1,0});
	}
}

void localProcesses3(State& st, Matrix& lumOutSy, Matrix& lumOutBr, Matrix& lumOutpp,
						Vector& energies, const int flags[], Matrix& lumOut)
{
	size_t jE=0;
	st.photon.ps.iterate([&](const SpaceIterator& itE) {
		energies[jE++] = itE.val(DIM_E);
	},{-1,0,0});
	
	#pragma omp parallel for
	for (int jE=0;jE<nE;jE++) {
		double frequency=energies[jE]/planck;
		size_t jR = 0.0;
		st.photon.ps.iterate([&](const SpaceIterator& itR) {
			double r = itR.val(DIM_R);
			double magf = st.magf.get(itR);
			double temp = st.tempElectrons.get(itR);
			double dens_e = st.denf_e.get(itR);
			double dens_i = st.denf_i.get(itR);
			
			double jv_th = jSync(energies[jE],temp,magf,dens_e)+jBremss(energies[jE],temp,dens_i,dens_e);
			//jv_th *= 4*pi;
			double jv_pl = luminositySynchrotron2(energies[jE],st.ntElectron,itR,magf)/frequency/(4*pi);
			double av_th = jv_th / bb(frequency,temp);
			SpaceCoord psc = {jE,jR,0};
			double av_pl = ssaAbsorptionCoeff(energies[jE],magf,st.ntElectron,psc);
			
			//double height = (r/schwRadius > 15) ? r*costhetaH(r) : height_fun(r);
			double height = r*costhetaH(r);
			double Sv = (jv_th+jv_pl)/(av_th+av_pl);
			//double Sv = jv_pl/(av_th+av_pl);
			double Inup = (frequency < 1.0e14) ? Sv * (1.0-exp(-2.0*height*(av_th+av_pl))) :
							2.0*height*(jv_th+jv_pl);
			//				2.0*height*jv_pl;
			double flux = 2.0*pi*r*r*(paso_r-1.0)*Inup;
			if (r/schwRadius <= 15)
				lumOutpp[jE][jR] = 4*pi*flux* (1.0-scattAA[jR][jR]);
			else
				lumOutBr[jE][jR] = 4*pi*flux* (1.0-scattAA[jR][jR]);
			lumOut[jE][jR] = lumOutpp[jE][jR]+lumOutBr[jE][jR];
			jR++;
		},{jE,-1,0});
	}
}


void reflectedSpectrum(Matrix lumOut, Matrix& lumOutRefl, Vector energies,
						double pasoE)
{
	size_t numInt = 40;
	for (size_t jE=nE/2;jE<nE;jE++) {
		double freq = energies[jE]/planck;
		double x = energies[jE]/(electronMass*cLight2);
		double xMax = 1.0e2;
		for (size_t jRcd=0;jRcd<nRcd;jRcd++) {
			double freq0,pasoFreq,aux1;
			int logicalInt = 1;
			if (x > 0.03 && x < xMax) {
				freq0 = energies[jE]/planck;
				pasoFreq = pow(xMax/x,1.0/(numInt+1));
				aux1 = 0.0;
			} else if (x <= 0.03) {
				freq0 = 0.03*electronMass*cLight2/planck;
				pasoFreq = pow(xMax/0.03,1.0/(numInt+1));
				aux1 = greenDeltaFunc(x);
			} else { logicalInt = 0; aux1 = 0.0; }
			
			double lum = 0.0;
			if (logicalInt == 1) {
				for (size_t jjFreq=0;jjFreq<numInt;jjFreq++) {
					double dfreq0 = freq0*(sqrt(pasoFreq)-pow(pasoFreq,-0.5));
					double lumInc = lumShellInc(freq0,jE,jRcd,lumOut,reachAD,energies,nR,redshift_RIAF_to_CD);
					double green = greenFuncRefl(freq,freq0);
					lum += green*lumInc*dfreq0*(freq/freq0);
					freq0 *= pasoFreq;
				}
			}
			lumOutRefl[jE][jRcd] =
					lum + aux1*lumShellInc(freq,jE,jRcd,lumOut,reachAD,energies,nR,redshift_RIAF_to_CD);
					
		}
	}
}

void coldDiskLuminosity(State& st, Matrix lumOut, Matrix& lumOutRefl, 
						Matrix& lumOutCD, Vector energies)
{
	matrixInit(lumOutCD,nE,nRcd,0.0);
	matrixInit(lumOutRefl,nE,nRcd,0.0);

	#pragma omp parallel for
	for (int jRcd=0;jRcd<nRcd;jRcd++) {
		double rCd = st.photon.ps[DIM_Rcd][jRcd];
		
		double lj1 = rCd/sqrt(paso_rCD);
		double lj2 = rCd*sqrt(paso_rCD);
		double area = 2.0*pi*(lj2*lj2-lj1*lj1);
		double aux1 = auxCD(rCd);
		double aux2 = 0.0;
		
		double pasoE = pow(energies[nE-1]/energies[0],1.0/(nE+1));
		reflectedSpectrum(lumOut,lumOutRefl,energies,pasoE);

		for (size_t jjE=0;jjE<nE;jjE++) {
			double frequency = energies[jjE]/planck;
			double dfreq = frequency*(sqrt(pasoE)-1.0/sqrt(pasoE));
			double lumInc = 0.0;
			for (size_t kR=0;kR<nR;kR++) {
				lumInc += (lumOut[jjE][kR]*reachAD[kR][jRcd] * pow(redshift_RIAF_to_CD[kR][jRcd],2));
			}
			aux2 += ((lumInc - lumOutRefl[jjE][jRcd])*dfreq);
		}
		double aux = aux1+aux2;
		double temp = pow(aux/area/ stefanBoltzmann,0.25);
		//cout << "Radius = " << rCd/schwRadius << " Temperature = " << scientific << temp << endl;

		for(size_t jE=0;jE<nE;jE++) {
			double frecuency=energies[jE]/planck;
			lumOutCD[jE][jRcd] = area * pi * bb(frecuency,1.7*temp)/pow(1.7,4);
		};
	}
}

void localCompton(const State& st, Matrix& lumOut, Matrix& lumOutIC, Vector energies)
{
	Matrix lumOutLocal;
	matrixInitCopy(lumOutLocal,nE,nR,lumOut);
	
	Vector comptonProbVec(nTempCompton*nNuPrimCompton*nNuCompton,0.0);
	if (calculateComptonRedMatrix) comptonMatrix();
	vectorRead("comptonProbMatrix.dat",comptonProbVec,comptonProbVec.size());
	
	// PARA CADA CELDA
	size_t jR=0;
	st.photon.ps.iterate([&](const SpaceIterator& itR) {

		double temp = st.tempElectrons.get(itR);
		double normtemp = boltzmann*temp/(electronMass*cLight2);
		if (normtemp >= tempMinCompton && normtemp < tempMaxCompton) {
			
			double eDens = st.denf_e.get(itR);
			double r = itR.val(DIM_R);
			double vol = volume(r);
			
			double tescap = 0.5*sqrt(pi)*height_fun(r)/cLight;
			Vector nOld(nE,0.0),nNew(nE,0.0),nWithoutCompton(nE,0.0);
			double tauu = eDens*thomson*height_fun(r);
			
			for (size_t jE=0;jE<nE;jE++) {
				nOld[jE] = lumOut[jE][jR]/vol / energies[jE] * tescap;
				nWithoutCompton[jE] = lumOut[jE][jR]/vol / energies[jE] * tescap;
			}
			
			double res = 0.0;
			size_t it=1;
			int cond = 1;
			do {
				cout << "jR = " << jR << "\t Iteration number = " << it << endl;
				cond = 0;
				
				for (size_t jE=0;jE<nE;jE++) {
					double frequency = energies[jE]/planck;
					if (frequency > nuMinCompton && frequency < nuMaxCompton) {
						double frequency = energies[jE]/planck;
						double rateScatt = rateThermal(normtemp,frequency)*eDens;
						double denom = 1.0/tescap + rateScatt;
						double integ = comptonNewLocal2(nOld,normtemp,jE,energies);
						nNew[jE] = (1.0/denom) * integ * eDens;
						if (nNew[jE] > nWithoutCompton[jE]) cond = 1;
					}
				}
				
				for (size_t jE=0;jE<nE;jE++) {
					nOld[jE] = nNew[jE];
					lumOut[jE][jR] += nNew[jE]*energies[jE]*vol / tescap;
					lumOutIC[jE][jR] = lumOut[jE][jR]-lumOutLocal[jE][jR];
				}
				cout << "Condition = " << cond << endl;
				++it;
			} while (cond == 1 && it < 40);
		}
		jR++;
	},{0,-1,0});
	
	show_message(msgEnd,Module_thermalCompton);
}

void thermalCompton2(State& st, Matrix& lumOut, Matrix& lumInICm, Matrix& lumOutIC, Vector energies, 
						int processesFlags[])
{
	show_message(msgStart,Module_thermalCompton);
	
	Vector tempVec(nTempCompton,0.0);
	Vector nuPrimVec(nNuPrimCompton,0.0);
	Vector nuVec(nNuCompton,0.0);
	Vector comptonProbVec(nTempCompton*nNuPrimCompton*nNuCompton,0.0);
	
	if (comptonMethod == 0) {
		if (calculateComptonRedMatrix) comptonMatrix2();
		vectorRead("comptonProbMatrix2.dat",comptonProbVec,comptonProbVec.size());
		vectorRead("tempComptonVec.dat",tempVec,tempVec.size());
	} else {
		if (calculateComptonRedMatrix) comptonMatrix();
		vectorRead("comptonProbMatrix.dat",comptonProbVec,comptonProbVec.size());
	}

	vectorRead("nuPrimComptonVec.dat",nuPrimVec,nuPrimVec.size());
	vectorRead("nuComptonVec.dat",nuVec,nuVec.size());
	
	Matrix lumOutLocal,lumOutCD,lumOutRefl;
	matrixInit(lumOutCD, nE, nRcd, 0.0);
	matrixInit(lumOutRefl, nE, nRcd, 0.0);
	matrixInit(lumOutIC, nE, nR, 0.0);
	Vector lumInIC(nE, 0.0);
	matrixInitCopy(lumOutLocal, nE, nR, lumOut);
	Vector p(nR*nE*nE, 0.0);
	
	Vector a(nR, 1.0);
	if (comptonMethod == 1)
		cNew(st, p, a);
		//cNew(st, p, redshift_to_inf);

	if (processesFlags[3])
		coldDiskLuminosity(st, lumOut, lumOutRefl, lumOutCD, energies);
		
	double res;		// Residual.
	size_t it=1;	// Iterations.
	do {
		cout << "Iteration number = " << it << endl;
		res=0.0;
		
		// To compute the residuals.
		Vector lumOld(nE,0.0);
		for (size_t jE=0;jE<nE;jE++) {
			for (size_t jR=0;jR<nR;jR++) {
				lumOld[jE] += lumOut[jE][jR];
			}
			if (processesFlags[3]) {
				for (size_t jRcd=0;jRcd<nRcd;jRcd++) {
					lumOld[jE] += lumOutCD[jE][jRcd] + lumOutRefl[jE][jRcd];
				}
			}
		}
		
		matrixInit(lumOutIC,nE,nR,0.0);
		
		// For each shell.
		size_t jR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itR) {
			double normtemp = boltzmann*st.tempElectrons.get(itR)/(electronMass*cLight2);
			// Compute the scattered luminosity.
			fill(lumInIC.begin(),lumInIC.end(),0.0);
			for (size_t jjE=0;jjE<nE;jjE++) {
				for (size_t jjR=0;jjR<nR;jjR++) {
					Vector lumVec(nE,0.0);
					//lumInIC[jjE] += scattAA[jjR][jR] * lumOut[jjE][jjR] * pow(redshift[jjR][jR],2);
					for (size_t jjjE=0;jjjE<nE;jjjE++)
						lumVec[jjjE] = lumOut[jjjE][jjR];
						
					double localEnergy = energies[jjE] / redshift[jjR][jR];
					lumInIC[jjE] += scattAA[jjR][jR] * pow(redshift[jjR][jR],2) 
										* lumInterp(lumVec,energies,jjE,nE,localEnergy);
				}
				if (processesFlags[3]) {
					for (size_t jRcd=0;jRcd<nRcd;jRcd++) {
						Vector lumVec(nE,0.0);
						for (size_t jjjE=0;jjjE<nE;jjjE++)
							lumVec[jjjE] = lumOutCD[jjjE][jRcd]+lumOutRefl[jjjE][jRcd];
						
						double localEnergy = energies[jjE] / redshift_CD_to_RIAF[jRcd][jR];
						lumInIC[jjE] += scattDA[jRcd][jR] * pow(redshift_CD_to_RIAF[jRcd][jR],2) *
											lumInterp(lumVec,energies,jjE,nE,localEnergy);
											
					}
				}
				lumInICm[jjE][jR] = lumInIC[jjE];
			}
			if (normtemp > tempMinCompton && normtemp < tempMaxCompton) {
				
				for (size_t jE=0;jE<nE;jE++) {
					double frequency = energies[jE]/planck;
					if (comptonMethod == 3) {
						if (it == 1) 
							comptonNewNewNewPruebaVector(jR,energies,jE,p,normtemp);
						if (frequency > nuMinCompton && frequency < nuMaxCompton)
							lumOutIC[jE][jR] = compton(p,lumInIC,jR,energies,jE);
					} else {
						if (frequency > nuMinCompton && frequency < nuMaxCompton)
							lumOutIC[jE][jR] = compton(p,lumInIC,jR,energies,jE);
					}
				}
			}
			jR++;
		},{0,-1,0});
		
		double pasoNuPrim = pow(energies[nE-1]/energies[0],1.0/(nE-1));
		double nPhNS = 0.0;
		double nPhBS = 0.0;
		double nPhAS = 0.0;
		for (size_t jE=0;jE<nE;jE++) {
			double dNu = energies[jE]/planck * (sqrt(pasoNuPrim)-1.0/sqrt(pasoNuPrim));
			double lumNew = 0.0;
			for (size_t jR=0;jR<nR;jR++) {
				nPhNS += lumOutLocal[jE][jR]/energies[jE] * dNu;
				nPhBS += lumInICm[jE][jR]/energies[jE] * dNu;
				nPhAS += lumOutIC[jE][jR]/energies[jE] * dNu;
				double lumAux = lumOutLocal[jE][jR]+lumOutIC[jE][jR];
				lumNew += lumAux;
				lumOut[jE][jR] = lumAux;
			}
			for (size_t jRcd=0;jRcd<nRcd;jRcd++) {
				lumNew += (lumOutCD[jE][jRcd]+lumOutRefl[jE][jRcd]);
			}
			if (lumNew > 0.0 && lumOld[jE] > 0.0)
				res += abs(log10(lumNew/lumOld[jE]));
		}
		cout << endl;
		cout << "Total photons non-scattered per unit time = " 
						<< nPhNS << " s^-1" << endl;
		cout << "Total photons scattered per unit time, before scattering = " 
						 << nPhBS << " s^-1" << endl;
		cout << "Total photons scattered per unit time, after scattering = " 
						 << nPhAS << " s^-1" << endl;
		cout << endl;
		
		res /= nE;
		cout << "Residuo = " << res << endl;
		++it;
	} while (res > 1.0e-3 && it < 100);

	show_message(msgEnd,Module_thermalCompton);
}

/*
void thermalCompton(State& st, Matrix& lumOut, Matrix& lumInICm, Matrix& lumOutIC, Vector energies, 
						int processesFlags[])
{
	show_message(msgStart,Module_thermalCompton);
	Vector tempVec(nTempCompton,0.0);
	Vector nuPrimVec(nNuPrimCompton,0.0);
	Vector nuVec(nNuCompton,0.0);
	Vector comptonProbVec(nTempCompton*nNuPrimCompton*nNuCompton,0.0);
	
	if (comptonMethod == 0) {
		if (calculateComptonRedMatrix) comptonMatrix2();
		vectorRead("comptonProbMatrix2.dat",comptonProbVec,comptonProbVec.size());
		vectorRead("tempComptonVec.dat",tempVec,tempVec.size());
	} else {
		if (calculateComptonRedMatrix) comptonMatrix();
		vectorRead("comptonProbMatrix.dat",comptonProbVec,comptonProbVec.size());
	}

	vectorRead("nuPrimComptonVec.dat",nuPrimVec,nuPrimVec.size());
	vectorRead("nuComptonVec.dat",nuVec,nuVec.size());
	
	Matrix lumOutCD,lumOutRefl,lumOutIC_local;
	matrixInit(lumOutCD,nE,nRcd,0.0);
	matrixInit(lumOutRefl,nE,nRcd,0.0);
	Vector lumInIC(nE,0.0);
	Vector lumBeforeCompton(nE,0.0);
	matrixInitCopy(lumOutIC_local,nE,nR,lumOut);
	Vector p(nR*nE*nE,0.0);
	
	for (size_t jE=0;jE<nE;jE++)
		for (size_t jR=0;jR<nR;jR++)
			lumBeforeCompton[jE] += lumOut[jE][jR];
	
	if (comptonMethod == 1) cNew(st,p,redshift_to_inf);
	
	int cond = 1;		// Residual.
	size_t it=1;	// Iterations.
	do {
		cout << "Iteration number = " << it << endl;
		cond = 0;
		if (processesFlags[3])
			coldDiskLuminosity(st,lumOut,lumOutRefl,lumOutCD,energies);

		Vector lumCompton(nE,0.0);
		
		// For each shell.
		size_t jR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itR) {
			double normtemp = boltzmann*st.tempElectrons.get(itR)/(electronMass*cLight2);
			// Compute the scattered luminosity.
			fill(lumInIC.begin(),lumInIC.end(),0.0);
			for (size_t jjE=0;jjE<nE;jjE++) {
				
				for (size_t jjR=0;jjR<nR;jjR++)
					lumInIC[jjE] += scattAA[jjR][jR] * lumOutIC_local[jjE][jjR];
				
				if (processesFlags[3]) {
					for (size_t jRcd=0;jRcd<nRcd;jRcd++)
						lumInIC[jjE] += scattDA[jRcd][jR] * 
											(lumOutCD[jjE][jRcd]+lumOutRefl[jjE][jRcd]);
				}
				lumInICm[jjE][jR] = lumInIC[jjE];
			}

			for (size_t jE=0;jE<nE;jE++)
				lumOutIC_local[jE][jR] = 0.0;

			if (normtemp > tempMinCompton && normtemp < tempMaxCompton) {
				double A = 1.0+4.0*normtemp*(1.0+4.0*normtemp);
				Vector unrEnergies_vec(nE);
				for (size_t jE=0;jE<nE;jE++)
					unrEnergies_vec[jE] = energies[jE]/redshift_to_inf[jR];
				
				for (size_t jE=0;jE<nE;jE++) {
					double frequency = energies[jE]/planck;
					double unrFrequency = unrEnergies_vec[jE]/planck;
					if (unrFrequency > nuMinCompton && unrFrequency < nuMaxCompton) {
						lumOutIC_local[jE][jR] = compton(p,lumInIC,jR,unrEnergies_vec,jE);
						lumCompton[jE] += lumOutIC_local[jE][jR];
					}
				}
			}
			jR++;
		},{0,-1,0});
		
		for (size_t jE=0;jE<nE;jE++)
			if (lumCompton[jE] > 0.01*lumBeforeCompton[jE]) cond = 1;
		
		
		double pasoNuPrim = pow(energies[nE-1]/energies[0],1.0/(nE-1));
		double nPhNS = 0.0;
		double nPhBS = 0.0;
		double nPhAS = 0.0;
		for (size_t jE=0;jE<nE;jE++) {
			double dNu = energies[jE]/planck * (sqrt(pasoNuPrim)-1.0/sqrt(pasoNuPrim));
			for (size_t jR=0;jR<nR;jR++) {
				nPhBS += lumInICm[jE][jR]/energies[jE] * dNu;
				nPhAS += lumOutIC_local[jE][jR]/energies[jE] * dNu;
				lumOut[jE][jR] += lumOutIC_local[jE][jR];
				lumOutIC[jE][jR] += lumOutIC_local[jE][jR];
			}
		}
		cout << endl;
		cout << "Total photons scattered per unit time, before scattering = " 
					 << nPhBS << " s^-1" << endl;
		cout << "Total photons scattered per unit time, after scattering = " 
					 << nPhAS << " s^-1" << endl;
		cout << endl;
		
		cout << "Condition = " << cond << endl;
		++it;
	} while (cond);
	show_message(msgEnd,Module_thermalCompton);
}
*/

void writeLuminosities(State& st, Vector energies, Matrix lumOutSy, Matrix lumOutBr,
						Matrix lumOutpp, Matrix lumInICm, Matrix lumOutIC, Matrix lumOut, 
						Matrix lumOutCD, Matrix lumOutRefl, const string& filename)
{
	ofstream file1,file2;
	ofstream fileCell;
	file1.open(filename.c_str(),ios::out);
	file2.open("lumRadius.txt",ios::out);
	fileCell.open("lumCell.txt",ios::out);
	
	double lumSy,lumBr,lumICin,lumIC,lumpp,lumTot,lumCD,lumRefl;
	double lumThermalTot = 0.0;
	double lumThermalTotRIAF = 0.0;
	double pasoF = pow(energies[nE-1]/energies[0],1.0/(nE-1));
	for (size_t jE=0;jE<nE;jE++) {
		double E = energies[jE];
		double frequency = E/planck;
		double energyEV = energies[jE]/EV_TO_ERG;
		lumSy = lumBr = lumICin = lumIC = lumpp = lumTot = lumCD = lumRefl = 0.0;
		for (size_t jR=0;jR<nR;jR++) {
			Vector lumSyVec(nE,0.0), lumBrVec(nE,0.0), lumICinVec(nE,0.0), lumICVec(nE,0.0),
					lumppVec(nE,0.0), lumOutVec(nE,0.0);
			for (size_t jjE=0;jjE<nE;jjE++) {
				lumSyVec[jjE] = lumOutSy[jjE][jR];
				lumBrVec[jjE] = lumOutBr[jjE][jR];
				lumICinVec[jjE] = lumInICm[jjE][jR];
				lumICVec[jjE] = lumOutIC[jjE][jR];
				lumppVec[jjE] = lumOutpp[jjE][jR];
				lumOutVec[jjE] = lumOut[jjE][jR];
			}
			double localEnergy = energies[jE]/redshift_to_inf[jR];
			lumSy += lumInterp(lumSyVec,energies,jE,nE,localEnergy) * pow(redshift_to_inf[jR],3) * escapeAi[jR];
			lumBr += lumInterp(lumBrVec,energies,jE,nE,localEnergy) * pow(redshift_to_inf[jR],3) * escapeAi[jR];
			lumICin += lumInterp(lumICinVec,energies,jE,nE,localEnergy) * pow(redshift_to_inf[jR],3) * escapeAi[jR];
			lumIC += lumInterp(lumICVec,energies,jE,nE,localEnergy) * pow(redshift_to_inf[jR],3) * escapeAi[jR];
			lumpp += lumInterp(lumppVec,energies,jE,nE,localEnergy) * pow(redshift_to_inf[jR],3) * escapeAi[jR];
			lumTot += lumInterp(lumOutVec,energies,jE,nE,localEnergy) * pow(redshift_to_inf[jR],3) * escapeAi[jR];
		}
		for (size_t jRcd=0;jRcd<nRcd;jRcd++) {
			double localEnergy = energies[jE]/redshift_CD_to_inf[jRcd];
			Vector lumCDVec(nE,0.0), lumReflVec(nE,0.0);
			for (size_t jjE=0;jjE<nE;jjE++) {
				lumCDVec[jjE] = lumOutCD[jjE][jRcd];
				lumReflVec[jjE] = lumOutRefl[jjE][jRcd];
			}
			lumCD += lumInterp(lumCDVec,energies,jE,nE,localEnergy) * pow(redshift_CD_to_inf[jRcd],3) * escapeDi[jRcd];
			lumRefl += lumInterp(lumReflVec,energies,jE,nE,localEnergy) * pow(redshift_CD_to_inf[jRcd],3) * escapeDi[jRcd];
		}
		double dfreq = frequency*(pasoF-1.0);
		lumThermalTot += (lumTot + lumCD + lumRefl)*dfreq;
		lumThermalTotRIAF += lumTot*dfreq;
        file1
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frequency
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << energyEV
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumSy*frequency
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumBr*frequency
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumIC*frequency
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumpp*frequency
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumCD*frequency
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumRefl*frequency
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumTot*frequency
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumICin*frequency
			<< endl;
	};
	cout << "Total thermal luminosity (power) = " << lumThermalTot << endl;
	
	double lumThermal = 0.0;
	for (size_t jR=0;jR<nR;jR++) {
		double r = st.denf_i.ps[DIM_R][jR]/schwRadius;
		for (size_t jE=0;jE<nE;jE++) {
			double E = energies[jE];
			double frequency = E/planck;
			double dfreq = frequency*(pasoF-1.0);
			lumThermal += lumOut[jE][jR]*dfreq*escapeAi[jR] * pow(redshift_to_inf[jR],3);
		}
		cout << "percentage of the lum produced inside r = " << r*sqrt(paso_r)
			 << " equal to " << lumThermal/lumThermalTotRIAF * 100 << " %" << endl;
	};
	
	double eVar = pow(energies[nE-1]/energies[0],1.0/nE);
	size_t jR=0;
	st.photon.ps.iterate([&](const SpaceIterator& itR) {
		double r = itR.val(DIM_R)/schwRadius;
		double vol = volume(itR.val(DIM_R));
		lumSy = lumBr = lumICin = lumIC = lumpp = lumTot = 0.0;
		for (size_t jE=0;jE<nE;jE++)  {
			double frequency = energies[jE]/planck;
			double dfrequency = frequency * (eVar-1.0);
			lumSy += lumOutSy[jE][jR]*dfrequency;
			lumBr += lumOutBr[jE][jR]*dfrequency;
			lumICin += lumInICm[jE][jR]*dfrequency;
			lumIC += lumOutIC[jE][jR]*dfrequency;
			lumpp += lumOutpp[jE][jR]*dfrequency;
			lumTot += lumOut[jE][jR]*dfrequency;
			fileCell << log10(r) << "\t" << log10(energies[jE]/planck)
					 << "\t" << safeLog10(lumOut[jE][jR]*frequency) << endl;
		}
		
		file2
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << r
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << vol
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumSy
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumBr
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumIC
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumpp
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumTot
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumICin
			<< endl;
		jR++;
	},{0,-1,0});
	
	file1.close();
	file2.close();
	fileCell.close();
}

void photonDensity(State& st, Vector energies, Matrix lumOut)
{
	ofstream file;
	file.open("photonDensity_gap.dat",ios::out);
	
	size_t nZ = nR;
	double zMin = schwRadius;
	double pasoZ = pow(1e6,1.0/nZ);
	double z = zMin;
	st.photon.ps.iterate([&](const SpaceIterator& iZ) {
		z *= pasoZ;
		double Uph = 0.0;
		size_t jE=0;
		double pasoE = pow(energies.back()/energies.front(),1.0/(energies.size()-1));
		st.photon.ps.iterate([&](const SpaceIterator& iZE) {
			double nPh = 0.0;
			double energy = iZE.val(DIM_E);
			double dE = energy * (pasoE - 1.0);
			size_t jR=0;
			st.photon.ps.iterate([&](const SpaceIterator& iZER) {
				double r = iZER.val(DIM_R);
				double dist2 = z*z+r*r;
				nPh += ( lumOut[jE][jR] / (4*pi*dist2*cLight*energy*planck) );
				jR++;
			},{iZE.coord[DIM_E],-1,0});
			Uph += nPh * energy * dE;
			file << z/schwRadius << "\t" << energy << "\t" << nPh << endl;
			st.photon.injection.set(iZE,nPh);
			jE++;
		},{-1,iZ.coord[DIM_R],0});
		double Ub = P2(magneticField(st.denf_e.ps[DIM_R][0])/pow(z/schwRadius,1))/(8*pi);
		cout << "z = " << z/schwRadius << " Rs,\t Uph = " << Uph << " erg cm^-3,\t Uph/Ub = " << Uph/Ub << endl;
	},{0,-1,0});
	file.close();
}

void photonDensityAux(State& st, Vector energies, Matrix lumOut)
{
	Vector z(50,GlobalConfig.get<double>("zGap")*schwRadius);
	double zMax = z[0]*1000;
	double pasoZ = pow(zMax/z[0],1.0/50);
	for (size_t i=1;i<50;i++) z[i] = z[i-1]*pasoZ;
	
	ofstream file;
	file.open("photonDensity_z.dat",ios::out);

	for (size_t i=1;i<50;i++) {
		double Uph = 0.0;
		size_t jE=0;
		double pasoE = pow(energies.back()/energies.front(),1.0/(energies.size()-1));
		st.photon.ps.iterate([&](const SpaceIterator& itE) {
			double nPh = 0.0;
			double energy = itE.val(DIM_E);
			double dE = energy * (pasoE - 1.0);
			size_t jR=0;
			st.photon.ps.iterate([&](const SpaceIterator& itER) {
				double r = itER.val(DIM_R);
				double dist2 = z[i]*z[i]+r*r;
				nPh += ( lumOut[jE][jR] / (4*pi*dist2*cLight*energy*planck) );
				jR++;
			},{itE.coord[DIM_E],-1,0});
			Uph += nPh * dE;
			jE++;
		},{-1,0,0});
		file << z[i]/schwRadius << "\t" << Uph << endl;
	}
	file.close();
}

void targetField(State& st, Matrix lumOut, Matrix lumCD, Matrix lumRefl, Vector energies)
{
	size_t jE=0;
	st.photon.ps.iterate([&](const SpaceIterator& itE) {
		size_t jR=0;
		double E = itE.val(DIM_E);
		st.photon.ps.iterate([&](const SpaceIterator& itER) {
			double r = itER.val(DIM_R);
			double rB2 = r*sqrt(paso_r);
			double rB1 = rB2/paso_r;
			double vol = volume(r);
			double lumReachingShell = 0.0;
			double height = height_fun(r);
			double tau_es_1 = st.denf_e.get(itER)*thomson*height;
			double tau_es_2 = tau_es_1 * (rB2-rB1)/height;
			double tescape = height/cLight * (1.0 + tau_es_1);
			double tcross = (rB2-rB1)/cLight * (1.0 + tau_es_2);
			
			for (size_t jjR=0;jjR<nR;jjR++) {
				Vector lumVec(nE,0.0);
				for (size_t jjE=0;jjE<nE;jjE++)
					lumVec[jjE] = lumOut[jE][jjR];
				
				double localEnergy = E / redshift[jjR][jR];
				double lumLocal = lumInterp(lumVec,energies,jE,nE,localEnergy);
				lumReachingShell += ( (jjR == jR) ? lumOut[jE][jR]*tescape : 
									reachAA[jjR][jR]*pow(redshift[jjR][jR],2)*lumLocal*tcross );
			}
			for (size_t jjRcd=0;jjRcd<nRcd;jjRcd++) {
				Vector lumVec(nE,0.0);
				for (size_t jjE=0;jjE<nE;jjE++)
					lumVec[jjE] = lumCD[jE][jjRcd] + lumRefl[jE][jjRcd];
				
				double localEnergy = E / redshift_CD_to_RIAF[jjRcd][jR];
				lumReachingShell += reachDA[jjRcd][jR]*pow(redshift_CD_to_RIAF[jjRcd][jR],2) * 
									lumInterp(lumVec,energies,jE,nE,localEnergy)*tcross;
			}
			st.photon.distribution.set(itER,lumReachingShell/(vol*planck*E)); //erg^⁻1 cm^-3 */
			//for (size_t jjR=0;jjR<nR;jjR++)
			//	lumReachingShell += reachAA[jjR][jR]*pow(redshift[jjR][jR],2)*lumOut[jE][jjR];
			//for (size_t jjRcd=0;jjRcd<nRcd;jjRcd++)
			//	lumReachingShell += reachDA[jjRcd][jR]*pow(redshift_CD_to_RIAF[jjRcd][jR],2)*
			//							(lumCD[jE][jjRcd]+lumRefl[jE][jjRcd]);
			//st.photon.distribution.set(itER, lumReachingShell/(4.0*pi*r*height*planck*E*cLight)); //erg^⁻1 cm^-3
			jR++;
		},{itE.coord[DIM_E],-1,0});
		jE++;
	},{-1,0,0});
}

void absorptionLumThermal(State& st, Matrix lumOut, Matrix lumCD, Matrix lumRefl, Matrix& lumOut_gg,
							Vector energies)
{
	int cond = 1;
	Matrix lumOut_gg_aux;
	matrixInitCopy(lumOut_gg,nE,nR,lumOut);
	matrixInitCopy(lumOut_gg_aux,nE,nR,lumOut);
	int it = 0;
	do {
		cond = 0;
		targetField(st,lumOut_gg,lumCD,lumRefl,energies);
		st.photon.ps.iterate([&](const SpaceIterator& iR) {
			double height = height_fun(iR.val(DIM_R));
			for (size_t jE=0;jE<nE;jE++) {
				double E = st.photon.ps[DIM_E][jE];
				double kappa_gg = integSimpsonLog(st.photon.emin(),st.photon.emax(),
						[&E,&iR,&st](double Eph)
						{
							if (Eph*E > P2(electronRestEnergy))
								return st.photon.distribution.interpolate({{0,Eph}},&iR.coord)*
										ggCrossSection2(E,Eph);
							else 
								return 0.0;
						},50);
				double tau_gg = 0.5*sqrt(pi)*kappa_gg*height;
				double factor_gg = (tau_gg > 1.0e-5) ? 
									(1.0-exp(-2*sqrt(3.0)*tau_gg))/(2.0*sqrt(3.0)*tau_gg) : 1.0;
				double lum = lumOut[jE][iR.coord[DIM_R]] * factor_gg;
				lumOut_gg[jE][iR.coord[DIM_R]] = lum;
				if (lum > 0.0 && lumOut[jE][iR.coord[DIM_R]] > 0.0)
					if (abs(safeLog10(lum/lumOut_gg_aux[jE][iR.coord[DIM_R]])) > 0.01) cond = 1;
				lumOut_gg_aux[jE][iR.coord[DIM_R]] = lum;
			}
		},{0,-1,0});
		it++;
		cout << "Iteration number " << it << "." << endl;
	} while (cond && it <= 20);
	cout << "Exit in " << it << " iterarations." << endl;
}
/*
void calculateElectronTemp(State& st, Matrix lumOut, Matrix lumInICm, Vector energies, double& res)
{
	ofstream fileNewTemp_e, fileCooling;
	fileNewTemp_e.open("newTempElectrons.txt");
	fileCooling.open("coolingElectrons.txt");
	double eVar = pow(energies[nE-1]/energies[0],1.0/(nE-1));
	size_t newDim = 10000;
	Vector logTe2(newDim,logTe.back());
	double paso_r_new = pow(st.denf_e.ps[DIM_R][nR-1]/st.denf_e.ps[DIM_R][1],1.0/(newDim-1));
	double r = st.denf_e.ps[DIM_R][nR-1];
	double dlogr = log(paso_r_new);
	size_t kR = 1;
	while (r > st.denf_e.ps[DIM_R][1]) {
		double Qmin = Qmin_func(r,lumOut,lumInICm,energies,st);
		double v = radialVel(r);
		double Te = exp(logTe2[kR-1]);
		double Ti = ionTemp(r);
		double pe = massDensityADAF(r) * boltzmann * Te / (atomicMassUnit*eMeanMolecularWeight);
		double Qie = qie_beta(r,Ti,Te);
		double Qp = Qplus(r,Ti,Te);
		double Qs = r/(pe*v) * (delta*Qp + Qie - Qmin);//*min(P3(Te/electronTemp(r)),1.0));
		double normtemp = boltzmann*Te / electronRestEnergy;
		double a_aux = 3.0-6.0/(4.0+5.0*normtemp) + normtemp * 30.0/P2(4.0+5.0*normtemp);
		double dlogrhodlogr = dlogrho_dlogr(r);
		double func = (Qs+dlogrhodlogr) / a_aux;
		
		while (abs(func*dlogr) > 0.01) {
			paso_r_new = pow(paso_r_new,1.0/10);
			dlogr = log(paso_r_new);
		}
		logTe2[kR] = min(log(min(electronTempOriginal(r),electronTemp(r))),
						logTe2[kR-1] - func * dlogr);
		//logTe2[kR] = logTe2[kR-1] - func *dlogr;
		fileNewTemp_e << log10(r/schwRadius) << "\t" << log10(exp(logTe2[kR])) << endl;
		fileCooling << log10(r/schwRadius) << "\t" << Qp << "\t" << Qie << "\t" << Qmin << endl;

		//cout << kR << endl;
		kR++;
		r /= paso_r_new;
		paso_r_new = pow(paso_r_new,10);
		dlogr = log(paso_r_new);
	}
	fileNewTemp_e.close();
	fileCooling.close();
}*/

void thermalRadiation(State& st, const string& filename)
{
	show_message(msgStart,Module_thermalLuminosities);

	Matrix lumOutSy,lumOutBr,lumOutpp,lumInICm,lumOutIC,lumOut,lumOutCD,lumOutRefl,lumOut_gg;
	Matrix lumSy,lumBr,lumIC,lum;
	matrixInit(lumOutSy,nE,nR,0.0);
	matrixInit(lumOutBr,nE,nR,0.0);
	matrixInit(lumOutpp,nE,nR,0.0);
	matrixInit(lumInICm,nE,nR,0.0);
	matrixInit(lumOutIC,nE,nR,0.0);
	matrixInit(lumOutCD,nE,nR,0.0);
	matrixInit(lumOutRefl,nE,nR,0.0);
	matrixInit(lumOut,nE,nR,0.0);
	matrixInit(lumOut_gg,nE,nR,0.0);

	Vector energies(nE,0.0);
	
	int processesFlags[numProcesses];	readThermalProcesses(processesFlags);
	if (processesFlags[0] || processesFlags[1] || processesFlags[2] || processesFlags[3]) {
		double res = 0.0;
		do {
			if (processesFlags[0] || processesFlags[1] || processesFlags[2])
				if (height_method == 0)
					localProcesses(st,lumOutSy,lumOutBr,lumOutpp,energies,processesFlags,lumOut);
				else {
					localProcesses2(st,lumOutSy,lumOutBr,energies,processesFlags,lumOut);
				}
			if (processesFlags[4]) {
				if (comptonMethod == 2)
					localCompton(st,lumOut,lumOutIC,energies);
				else
					thermalCompton2(st,lumOut,lumInICm,lumOutIC,energies,processesFlags);
			}
			if (processesFlags[3])
				coldDiskLuminosity(st,lumOut,lumOutRefl,lumOutCD,energies);
		} while (res > 1.0e-3);
		//photonDensity(st,energies,lumOut);
		//photonDensityAux(st,energies,lumOut);
		absorptionLumThermal(st,lumOut,lumOutCD,lumOutRefl,lumOut_gg,energies);
		writeLuminosities(st,energies,lumOutSy,lumOutBr,lumOutpp,lumInICm,
							lumOutIC,lumOut_gg,lumOutCD,lumOutRefl,filename);
	}	
	show_message(msgEnd,Module_thermalLuminosities);
}