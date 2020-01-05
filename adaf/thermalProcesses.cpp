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
						Vector& energies, Vector& redshift, const int flags[], Matrix& lumOut)
{
	size_t jR=0;
	st.photon.ps.iterate([&](const SpaceIterator& itR) {
		double r = itR.val(DIM_R);
		double redshift_factor = sqrt(1.0-schwRadius/r)*(1.0-P2(radialVel(r)/cLight));
		redshift_factor = (redshift_factor > 0.0) ? redshift_factor : 1.0;
		redshift[jR] = redshift_factor;
		redshift[jR] = 1.0;
		jR++;
	},{0,-1,0});
	
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

void localProcesses2(State& st, Matrix& lumOutSy, Matrix& lumOutBr, Matrix& lumOutpp,
						Vector& energies, Vector& redshift, const int flags[], Matrix& lumOut)
{
	size_t jE=0;
	st.photon.ps.iterate([&](const SpaceIterator& itE) {
		energies[jE++] = itE.val(DIM_E);
	},{-1,0,0});
	
	size_t jR=0;
	st.photon.ps.iterate([&](const SpaceIterator& itR) {
		double r = itR.val(DIM_R);
		double redshift_factor = sqrt(1.0-schwRadius/r)*(1.0-P2(radialVel(r)/cLight));
		redshift_factor = (redshift_factor > 0.0) ? redshift_factor : 1.0;
		redshift[jR] = redshift_factor;
		redshift[jR] = 1.0;
		jR++;
	},{0,-1,0});

	#pragma omp parallel for
	for (int jE=0;jE<nE;jE++) {
		size_t jR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itER) {
			double r = itER.val(DIM_R);
			double temp = st.tempElectrons.get(itER);
			double magf = st.magf.get(itER);
			double dens_i = st.denf_i.get(itER);
			double dens_e = st.denf_e.get(itER);
			double xSy,xBr;
			double unrEnergy = energies[jE]/redshift[jR];
			xSy = xBr = 0.0;
			if (flags[0])
				xSy = jSync(unrEnergy,temp,magf,dens_e)*4.0*pi;
			if (flags[1])
				xBr = jBremss(unrEnergy,temp,dens_i,dens_e)*4.0*pi;
			
			SpaceCoord psc = {jE,jR,0};
			double b_nu = bb(unrEnergy/planck,temp);
			double kappa = (xSy+xBr)/(4.0*pi*b_nu);
			double tau = 0.5*sqrt(pi)*kappa*height_fun(r);
			double flux = (2.0*sqrt(3.0)*tau > 1.0e-5) ? 
							2.0*pi/sqrt(3.0)*b_nu*(1.0-exp(-2.0*sqrt(3.0)*tau)) :
							0.5*sqrt(pi)*(xSy+xBr)*height_fun(r);
			lumOutSy[jE][jR] = 2.0*pi*flux*r*r*(sqrt(paso_r)-sqrt(1.0/paso_r));
			lumOut[jE][jR] = lumOutSy[jE][jR];
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
					double lumInc = lumShellInc(freq0,jE,jRcd,lumOut,reachAD,energies,nR);
					double green = greenFuncRefl(freq,freq0);
					lum += green*lumInc*dfreq0*(freq/freq0);
					freq0 *= pasoFreq;
				}
			}
			lumOutRefl[jE][jRcd] =
					lum + aux1*lumShellInc(freq,jE,jRcd,lumOut,reachAD,energies,nR);
					
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
		
		double lj = rCd/sqrt(paso_rCD);
		double lj1 = rCd*sqrt(paso_rCD);
		double area = 2.0*pi*(lj1*lj1-lj*lj);
		double aux1 = auxCD(rCd);
		double aux2 = 0.0;
		
		double pasoE = pow(energies[nE-1]/energies[0],1.0/(nE+1));
		reflectedSpectrum(lumOut,lumOutRefl,energies,pasoE);

		for (size_t jjE=0;jjE<nE;jjE++) {
			double frequency = energies[jjE]/planck;
			double dfreq = frequency*(sqrt(pasoE)-1.0/sqrt(pasoE));
			double lumInc = 0.0;
			for (size_t kR=0;kR<nR;kR++) {
				lumInc += (lumOut[jjE][kR]*reachAD[kR][jRcd]);
			}
			aux2 += ((lumInc - lumOutRefl[jjE][jRcd])*dfreq);
		}
		double aux = aux1+aux2;
		double temp = pow(aux/area/ stefanBoltzmann,0.25);
		//cout << "Radius = " << rCd/schwRadius << " Temperature = " << scientific << temp << endl;

		for(size_t jE=0;jE<nE;jE++) {
			double frecuency=energies[jE]/planck;
			lumOutCD[jE][jRcd] = area * bb(frecuency,1.7*temp)/pow(1.7,4);
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
					//if (nNew[jE] > 0.0 && nOld[jE] > 0.0)
					//	res += abs(log10(nNew[jE]/nOld[jE]));
				}
				
				for (size_t jE=0;jE<nE;jE++) {
					nOld[jE] = nNew[jE];
					lumOut[jE][jR] += nNew[jE]*energies[jE]*vol / tescap;
					lumOutIC[jE][jR] = lumOut[jE][jR]-lumOutLocal[jE][jR];
				}
				/*
				for (size_t jE=0;jE<nE;jE++) {
					double frequency=energies[jE]/planck;
					if (frequency > nuMinCompton && frequency < nuMaxCompton) {
						double Rprom = eDens*rateThermal(normtemp,frequency);
						lumNew[jE] = pow(1.0/tescap + Rprom,-1) * 
							(comptonNewLocal(comptonProbVec,lumOld,normtemp,jE,
								energies,nE)*eDens+cLight*lumOutLocal[jE][jR]*(area/vol));
					}
					if (lumNew[jE] > 0.0 && lumOld[jE] > 0.0 && it > 1)
						res += abs(log10(lumNew[jE]/lumOld[jE]));
				}
				for (size_t jE=0;jE<nE;jE++) {
					lumOld[jE] = lumNew[jE];
					lumOut[jE][jR] = lumNew[jE];
					lumOutIC[jE][jR] = lumNew[jE]-lumOutLocal[jE][jR];
				}*/
				cout << "Condition = " << cond << endl;
				++it;
			//} while (res > 1.0e-3 && it < 30);
			} while (cond == 1 && it < 40);
		}
		jR++;
	},{0,-1,0});
	
	show_message(msgEnd,Module_thermalCompton);
}

void thermalCompton(State& st, Matrix& lumOut, Matrix& lumInICm, Matrix& lumOutIC, Vector energies, 
						Vector redshift, int processesFlags[])
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
	matrixInit(lumOutCD,nE,nRcd,0.0);
	matrixInit(lumOutRefl,nE,nRcd,0.0);
	Vector lumInIC(nE,0.0);
	matrixInitCopy(lumOutLocal,nE,nR,lumOut);
	Vector p(nR*nE*nE,0.0);
	
	if (comptonMethod == 1) cNew(st,p,redshift);
	
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
					lumOld[jE] += lumOutCD[jE][jRcd]+lumOutRefl[jE][jRcd];
				}
			}
		}
		
		matrixInit(lumOutIC,nE,nR,0.0);
		
		if (processesFlags[3])
			coldDiskLuminosity(st,lumOut,lumOutRefl,lumOutCD,energies);

		// For each shell.
		size_t jR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itR) {
			double normtemp = boltzmann*st.tempElectrons.get(itR)/(electronMass*cLight2);
			// Compute the scattered luminosity.
			fill(lumInIC.begin(),lumInIC.end(),0.0);
			for (size_t jjE=0;jjE<nE;jjE++) {
				// K-N correction:
				//double Rprom = rateThermal(normtemp,energies[jjE])/(thomson*cLight);
				//double Rprom=1.0;
				//if (Rprom > 0.0 && Rprom < 1.001) {
				for (size_t jjR=0;jjR<nR;jjR++)
					lumInIC[jjE] += scattAA[jjR][jR] * lumOut[jjE][jjR];
				
				if (processesFlags[3]) {
					for (size_t jRcd=0;jRcd<nRcd;jRcd++)
						lumInIC[jjE] += scattDA[jRcd][jR] * 
											(lumOutCD[jjE][jRcd]+lumOutRefl[jjE][jRcd]);
				}
				lumInICm[jjE][jR] = lumInIC[jjE];
			}
			if (normtemp > tempMinCompton && normtemp < tempMaxCompton) {
				Vector unrEnergies_vec(nE);
				for (size_t jE=0;jE<nE;jE++)
					unrEnergies_vec[jE] = energies[jE]/redshift[jR];
				
				for (size_t jE=0;jE<nE;jE++) {
					double frequency = energies[jE]/planck;
					double unrFrequency = unrEnergies_vec[jE]/planck;
					if (comptonMethod == 3) {
						if (it == 1) comptonNewNewNewPruebaVector(jR,unrEnergies_vec,jE,p,normtemp);
						if (unrFrequency > nuMinCompton && unrFrequency < nuMaxCompton)
							lumOutIC[jE][jR] = comptonNewNewPrueba(p,lumInIC,jR,unrEnergies_vec,jE);
							//lumOutIC[jE][jR] = comptonNewNewNew(lumInIC,normtemp,frequency,
							//						energies,jE);
					} else {
						//if (it==1) comptonNewNewPruebaVector(tempVec,nuPrimVec,nuVec,comptonProbVec,jR,
						//				energies,jE,p,normtemp);
					
						//if (frequency > nuMinCompton && frequency < nuMaxCompton)
						//	lumOutIC[jE][jR] = comptonNewNew(tempVec,nuPrimVec,nuVec,
						//								comptonProbVec,lumInIC,normtemp,
						//								frequency,energies,jE);
						if (unrFrequency > nuMinCompton && unrFrequency < nuMaxCompton)
							lumOutIC[jE][jR] = comptonNewNewPrueba(p,lumInIC,jR,unrEnergies_vec,jE);
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
		cout << "Total photons scattered per unit time, after scattering = " 
					 << nPhBS << " s^-1" << endl;
		cout << "Total photons scattered per unit time, after scattering = " 
					 << nPhAS << " s^-1" << endl;
		cout << endl;
		
		res /= nE;
		cout << "Residuo = " << res << endl;
		++it;
	} while (res > 1.0e-3 && it < 100);
	//} while (it < 3);
	show_message(msgEnd,Module_thermalCompton);
}

/*
void gravRedshift(const State& st, Vector energies, Matrix lum, Matrix& lumRed)
{
	matrixInit(lumRed,nE,nR,0.0);

	//DEFINO FRONTERAS PARA LAS CELDAS DE ENERGIA
	double eVar = pow(energies[nE-1]/energies[0],1.0/(nE-1));
	Vector benergies(nE+1,0.0);
	benergies[0] = energies[0]/sqrt(eVar);
	for (size_t jE=0;jE<nE;jE++){
		benergies[jE+1] = energies[jE]*sqrt(eVar);
	}
	
	size_t jR=0;
	st.photon.ps.iterate([&](const SpaceIterator& itR) {
		double r = itR.val(DIM_R);
		double redshift = sqrt(1.0-schwRadius/r)*(1.0-P2(radialVel(r)/cLight));
		redshift = (redshift > 0.0) ? redshift : 1.0;
		size_t items[nE];
		for (size_t jE=0;jE<nE;jE++) {
			double reddendEnergy = energies[jE]/redshift;
			size_t count=0;
			size_t jjE=0;
			while(jjE < nE && count == 0) {
				if (reddendEnergy > benergies[jjE] && renergy < benergies[jjE+1]) {
					items[jE] = jjE;
					count++;
				}
				jjE++;
			}
			if (count == 0) items[jE]=nE;
		}
		for (size_t jE=0;jE<nE;jE++) {
			if (items[jE] < nE) {
				double energy_new = energies[items[jE]];
				lumRed[items[jE]][jR] += lum[jE][jR] * (gfactor*reddendEnergy/energy_new);
			}
		}
		jR++;
	},{0,-1,0});
}*/

/*void binEmissivities(const State& st, Vector esc, Vector energies, Matrix lumOut)
{
	ofstream filePar, fileLum;
	FILE *fileEmissBin;
	
	fileLum.open("emissivities.txt", ios::out);
	filePar.open("parSpace.txt", ios::out);
	fileEmissBin=fopen("emissivities.bin","wb");

	float *emissivities;
	emissivities=(float*)calloc(nE*nR*nTheta,sizeof(float));
	
	size_t jj=0;
	double eVar = pow(energies[nE-1]/energies[0],1.0/nE);
	for (size_t jE=0;jE<nE;jE++) {
		double sum=0.0;
        double denergy=energies[jE]*(eVar-1.0);
        size_t jR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itR) {
			double r = itR.val(DIM_R);
			double rB1=r-dr/2.0;
			double rB2=r+dr/2.0;
			double vol=(rB2*rB2*rB2-rB1*rB1*rB1)*volFactor;
			double emissToLum=vol*4.0*pi;
			double sumFactor=0.0;
			st.photon.ps.iterate([&](const SpaceIterator& itTheta) {
				sumFactor += st.tempElectrons.get(itTheta);
			},{0,itR.coord[DIM_R],-1});
			st.photon.ps.iterate([&](const SpaceIterator& itTheta) {
				double lum = lumOut[jE][jR]*esc[jR]*st.tempElectrons.get(itTheta)/(sumFactor*2.0);
				fileLum << jE << "\t" << jR << "\t" 
						  << lum << endl;
				emissivities[jj++] = (float)(lum/emissToLum);
			},{0,itR.coord[DIM_R],-1});
            jR++;
		},{0,-1,0});
	}
	fwrite(emissivities,sizeof(float),jj,fileEmissBin);
	
	float deltaLogE = log10(energies[1]/energies[0]);
	filePar << nE << "\t" << nR << "\t" << nTheta << endl;
	filePar << cuspRadius << "\t" << dr << "\t" << minPolarAngle << "\t"
			<< dtheta << "\t" << log10(energies[0]) << "\t" << deltaLogE << endl;
	filePar << gravRadius << "\t" << blackHoleSpinPar << "\t" 
			<< specificAngMom << "\t" << torusCenterRadius << endl;
	
	fclose(fileEmissBin);
	fileLum.close();
	filePar.close();
}*/

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
	double pasoF = pow(energies[nE-1]/energies[0],1.0/(nE-1));
	for (size_t jE=0;jE<nE;jE++) {
		double frequency = energies[jE]/planck;
		double energyEV = energies[jE]/EV_TO_ERG;
		lumSy = lumBr = lumICin = lumIC = lumpp = lumTot = lumCD = lumRefl = 0.0;
		for (size_t jR=0;jR<nR;jR++) {
			lumSy += lumOutSy[jE][jR] * escapeAi[jR];
			lumBr += lumOutBr[jE][jR] * escapeAi[jR];
			lumICin += lumInICm[jE][jR] * escapeAi[jR];
			lumIC += lumOutIC[jE][jR] * escapeAi[jR];
			lumpp += lumOutpp[jE][jR] * escapeAi[jR];
			lumTot += lumOut[jE][jR] * escapeAi[jR];
		}
		for (size_t jRcd=0;jRcd<nRcd;jRcd++) {
			lumCD += lumOutCD[jE][jRcd] * escapeDi[jRcd];
			lumRefl += lumOutRefl[jE][jRcd] * escapeDi[jRcd];
		}
		double dfreq = frequency*(pasoF-1.0);
		lumThermalTot += (lumTot + lumCD + lumRefl)*dfreq;
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
	double z = GlobalConfig.get<double>("zGap")*schwRadius;
	
	ofstream file;
	file.open("photonDensity1.txt",ios::out);

	size_t jE=0;
	st.photon.ps.iterate([&](const SpaceIterator& itE) {
		double nPh = 0.0;
		double energy = itE.val(DIM_E);
		size_t jR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itER) {
			double r = itER.val(DIM_R);
			double dist2 = z*z+r*r;
			nPh += ( lumOut[jE][jR] / (4*pi*dist2*cLight*energy*planck) );
			jR++;
		},{itE.coord[DIM_E],-1,0});
		file << safeLog10(energy) << "\t" << safeLog10(nPh) << endl;
		jE++;
	},{-1,0,0});
	file.close();
}

void targetField(State& st, Matrix lumOut, Matrix lumCD, Matrix lumRefl)
{
	size_t jE=0;
	st.photon.ps.iterate([&](const SpaceIterator& itE) {
		size_t jR=0;
		double E = itE.val(DIM_E);
		st.photon.ps.iterate([&](const SpaceIterator& itER) {
			double r = itER.val(DIM_R);
			double vol = volume(r);
			double lumReachingShell = 0.0;
			double tcross = r*(sqrt(paso_r)-1.0/sqrt(paso_r))/cLight;
			for (size_t jjR=0;jjR<nR;jjR++)
				lumReachingShell += reachAA[jjR][jR]*lumOut[jE][jjR];
			for (size_t jjRcd=0;jjRcd<nRcd;jjRcd++)
				lumReachingShell += reachDA[jjRcd][jR]*(lumCD[jE][jjRcd]+lumRefl[jE][jjRcd]);
			st.photon.distribution.set(itER,lumReachingShell*tcross/(vol*planck*E)); //erg^⁻1 cm^-3
			jR++;
		},{itE.coord[DIM_E],-1,0});
		jE++;
	},{-1,0,0});
}

/*
void writeBlob(State& st, Vector energies, Matrix lum, double t)
{
	double IRfreq = 1.4e14;
	double frequency1 = energies[0]/planck;
	double frequency2 = frequency1;
	size_t jE = 0;
	while (frequency2 < IRfreq) {
		frequency1 = frequency2;
		frequency2 = energies[++jE]/planck;
	}
	frequency1 = sqrt(frequency1*frequency2);
	double IRlum = 0.0;
	for (size_t jR=0;jR<nR;jR++)
		IRlum += sqrt(lum[jE-1][jR]*lum[jE][jR]) * escapeAi[jR];
	
	IRlum *= frequency1;
	ofstream file;
	file.open("lightCurveBlob.txt",ios::app);
	file << t << "\t" << IRlum << endl;
	file.close();
}*/

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
}

void thermalProcesses(State& st, const string& filename)
{
	show_message(msgStart,Module_thermalLuminosities);

	Matrix lumOutSy,lumOutBr,lumOutpp,lumInICm,lumOutIC,lumOut,lumOutCD,lumOutRefl;
	matrixInit(lumOutSy,nE,nR,0.0);
	matrixInit(lumOutBr,nE,nR,0.0);
	matrixInit(lumOutpp,nE,nR,0.0);
	matrixInit(lumInICm,nE,nR,0.0);
	matrixInit(lumOutIC,nE,nR,0.0);
	matrixInit(lumOutCD,nE,nR,0.0);
	matrixInit(lumOutRefl,nE,nR,0.0);
	matrixInit(lumOut,nE,nR,0.0);

	Vector energies(nE,0.0);
	Vector redshift(nR,0.0);
	
	int processesFlags[numProcesses];	readThermalProcesses(processesFlags);
	if (processesFlags[0] || processesFlags[1] || processesFlags[2] || processesFlags[3]) {
		double res = 0.0;
		do {
			if (processesFlags[0] || processesFlags[1] || processesFlags[2])
				if (height_method == 0)
					localProcesses(st,lumOutSy,lumOutBr,lumOutpp,energies,redshift,processesFlags,lumOut);
				else
					localProcesses2(st,lumOutSy,lumOutBr,lumOutpp,energies,redshift,processesFlags,lumOut);
			if (processesFlags[4]) {
				if (comptonMethod == 2)
					localCompton(st,lumOut,lumOutIC,energies);
				else
					thermalCompton(st,lumOut,lumInICm,lumOutIC,energies,redshift,processesFlags);
			}
			if (processesFlags[3])
				coldDiskLuminosity(st,lumOut,lumOutRefl,lumOutCD,energies);
		} while (res > 1.0e-3);
		photonDensity(st,energies,lumOut);
		writeLuminosities(st,energies,lumOutSy,lumOutBr,lumOutpp,lumInICm,
							lumOutIC,lumOut,lumOutCD,lumOutRefl,filename);
		targetField(st,lumOut,lumOutCD,lumOutRefl);
		if (calculateNewTemp == 1)
			calculateElectronTemp(st,lumOut,lumInICm,energies,res);
	}	
	show_message(msgEnd,Module_thermalLuminosities);
}