#include "secondariesProcesses.h"
#include "thermalCompton.h"
#include "absorption.h"
#include "modelParameters.h"
#include "write.h"
#include "messages.h"
#include "globalVariables.h"
#include "adafFunctions.h"
#include "NTtimescales.h"

#include <fmath/RungeKutta.h>

#include <fluminosities/thermalSync.h>
#include <fluminosities/blackBody.h>
#include <fluminosities/luminositySynchrotron.h>
#include <fluminosities/luminosityPhotoHadronic.h>
#include <fluminosities/luminosityNTHadronic.h>
#include <fluminosities/opticalDepthSSA.h>
#include <fluminosities/luminosityIC.h>
#include <fluminosities/luminosityPairAnnihilation.h>
#include "NTinjection.h"
#include "NTdistribution.h"
#include <fparameters/Dimension.h>
#include <fparameters/SpaceIterator.h>

void secondariesRadiationProcesses(State& st, const std::string& filename)
{
	show_message(msgStart, Module_secondariesLuminosities);

	std::ofstream file;
	file.open(filename.c_str(),std::ios::out);
	file << "log(E/eV)"
		 << '\t' << "eSyn"
		 << "\t" << "muSy"
		 << "\t" << "piSy"
		 << "\t" << "eIC"
		 << "\t" << "piPP"
		 << "\t" << "piPG"
		 << "\t" << "totAbs"
		 << std::endl;
	
	double Ephmin = st.photon.emin();
	double Ephmax = st.photon.emax();
	
	Particle &ntPh = st.ntPhoton;
	size_t nEnt = ntPh.ps[DIM_E].size();
	Vector ntEnergies(nEnt,0.0);
	for (size_t jE=0;jE<nEnt;jE++)
		ntEnergies[jE] = ntPh.ps[DIM_E][jE];
	
	Matrix eSy;					matrixInit(eSy,nEnt,nR,0.0);
	Matrix eIC;					matrixInit(eIC,nEnt,nR,0.0);
	Matrix muSy;				matrixInit(muSy,nEnt,nR,0.0);
	Matrix piSy;				matrixInit(piSy,nEnt,nR,0.0);
	Matrix eeAnn;				matrixInit(eeAnn,nEnt,nR,0.0);
	Matrix piPP;				matrixInit(piPP,nEnt,nR,0.0);
	Matrix piPG;				matrixInit(piPG,nEnt,nR,0.0);
	Matrix totAbs;				matrixInit(totAbs,nEnt,nR,0.0);
	
	#pragma omp parallel for
	for (int E_ix=0;E_ix<nEnt;E_ix++) {
		double E = ntEnergies[E_ix];
		size_t jR = 0;
		st.ntPhoton.ps.iterate([&](const SpaceIterator& iR) {
			double r = iR.val(DIM_R);
			double rB1 = r/sqrt(paso_r);
			double rB2 = r*sqrt(paso_r);
			double vol = volume(r);
			double magf = st.magf.get(iR);
			SpaceCoord psc = {E_ix,iR.coord[DIM_R],0};
			
			double height = height_fun(r);
			double taugg = st.tau_gg.get(psc);
			double kappagg = 2.0/sqrt(pi)*taugg/height;
			
			double eICLocal,eSyLocal,piPPLocal,piPGLocal,piSyLocal,muSyLocal,ePairAnnLocal;
			eICLocal = eSyLocal = piPPLocal = piPGLocal = piSyLocal = muSyLocal = 0.0;
			eICLocal = luminosityIC_2(E,st.ntPair,iR.coord,st.photon.distribution,Ephmin,Ephmax)/E;
			eICLocal += luminosityIC_2(E,st.ntPair,iR.coord,st.ntPhoton.distribution,
										st.ntPhoton.emin(),st.ntPhoton.emax())/E;
			//ePairAnnLocal += luminosityPairAnnihilation(Elocal,st.ntPair,iR)/Elocal;
			eSyLocal = luminositySynchrotronExactSec(E,st.ntPair,iR,magf)/E;
			muSyLocal = luminositySynchrotronExact(E,st.ntMuon,iR,magf)/E;
			piSyLocal = luminositySynchrotronExact(E,st.ntChargedPion,iR,magf)/E;
			
			if (E/EV_TO_ERG > 1.0e7) {	
				piPGLocal = luminosityPhotoHadronic(E,st.ntChargedPion,st.photon.distribution,iR,Ephmin,Ephmax)/E;
				piPPLocal = luminosityNTHadronic(E,st.ntChargedPion,st.denf_i.get(iR),iR)/E;
			}
			
			double alpha_th = jSync(E,st.tempElectrons.get(iR),magf,st.denf_e.get(iR)) / 
									bb(E/planck,st.tempElectrons.get(iR));
			double alpha_ssa = ssaAbsorptionCoeff(E,magf,st.ntElectron,psc);
			double alpha_tot = alpha_th + alpha_ssa;
			double tau = 0.5*sqrt(pi)*alpha_tot*height;
			
			double eSyLocalFlux = (tau > 1.0e-10) ?
								0.5/sqrt(3.0) * (eSyLocal/alpha_tot) * 
								(1.0-exp(-2.0*sqrt(3.0)*tau)) :
								eSyLocal * 0.5 * sqrt(pi) * height;  // flux
			eSy[E_ix][jR]  = eSyLocalFlux * 2.0 * pi*(rB2*rB2-rB1*rB1);  // luminosity
				

			double factor = pi/sqrt(3.0)*(rB2*rB2-rB1*rB1) / vol;
			double eTot = eSyLocal+eICLocal+muSyLocal+piSyLocal+piPPLocal+piPGLocal+ePairAnnLocal;

			//eSy[E_ix][jR] = eSyLocal*vol;
			muSy[E_ix][jR] = muSyLocal*vol;
			piSy[E_ix][jR] = piSyLocal*vol;
			eIC[E_ix][jR] = eICLocal*vol;
			piPP[E_ix][jR] = piPPLocal*vol;
			piPG[E_ix][jR] = piPGLocal*vol;
			eeAnn[E_ix][jR] = ePairAnnLocal*vol;
			
			double totLocal = (taugg+tau > 1.0e-10) ? factor*vol*eTot/(alpha_tot+kappagg) * 
								(1.0-exp(-2.0*sqrt(3.0)*(taugg+tau))) : 
								eTot*vol;
			//totLocal = eTot*vol*exp(-taugg);
			
			double totLocalNA = eTot;
			totAbs[E_ix][jR] = totLocal;
			
			double tau_es = st.denf_e.get(iR)*thomson*height;
			double tescape = height/cLight * (1.0+tau_es);
			SpaceCoord iE = {E_ix,iR.coord[DIM_R],0};
			double rate_gg = kappagg * cLight;
			double NTphot = st.ntPhoton.injection.get(iE);
			
			st.ntPhoton.distribution.set(iE,NTphot + eTot/E * pow(rate_gg+pow(tescape,-1),-1));
			jR++;
		},{E_ix,-1,0});
	}
	
	
	double lum_eSy, lum_muSy, lum_piSy, lum_eIC, lum_piPP, lum_piPG, lum_eeAnn, lum_totAbs;
	for (size_t jE=0;jE<nEnt;jE++) {
		lum_eSy = lum_muSy = lum_piSy = lum_eIC = lum_piPP = lum_piPG = lum_eeAnn = lum_totAbs = 0.0;
		double E = ntEnergies[jE];
		double fmtE = safeLog10(E/EV_TO_ERG);
		for (size_t jR=0;jR<nR;jR++) {
			Vector eSyVec(nEnt,0.0), muSyVec(nEnt,0.0), piSyVec(nEnt,0.0), eICVec(nEnt,0.0), piPPVec(nEnt,0.0),
					piPGVec(nEnt,0.0), eeAnnVec(nEnt,0.0), totAbsVec(nEnt,0.0);
			for (size_t jjE=0;jjE<nEnt;jjE++) {
				eSyVec[jjE] = eSy[jjE][jR];
				muSyVec[jjE] = muSy[jjE][jR];
				piSyVec[jjE] = piSy[jjE][jR];
				eICVec[jjE] = eIC[jjE][jR];
				piPPVec[jjE] = piPP[jjE][jR];
				piPGVec[jjE] = piPG[jjE][jR];
				eeAnnVec[jjE] = eeAnn[jjE][jR];
				totAbsVec[jjE] = totAbs[jjE][jR];
			}
			double localEnergy = E/redshift_to_inf[jR];
			lum_eSy += lumInterp(eSyVec,ntEnergies,jE,nEnt,localEnergy) * pow(redshift_to_inf[jR],3) * escapeAi[jR];
			lum_muSy += lumInterp(muSyVec,ntEnergies,jE,nEnt,localEnergy) * pow(redshift_to_inf[jR],3) * escapeAi[jR];
			lum_muSy += lumInterp(piSyVec,ntEnergies,jE,nEnt,localEnergy) * pow(redshift_to_inf[jR],3) * escapeAi[jR];
			lum_eIC += lumInterp(eICVec,ntEnergies,jE,nEnt,localEnergy) * pow(redshift_to_inf[jR],3) * escapeAi[jR];
			lum_piPP += lumInterp(piPPVec,ntEnergies,jE,nEnt,localEnergy) * pow(redshift_to_inf[jR],3) * escapeAi[jR];
			lum_piPG += lumInterp(piPGVec,ntEnergies,jE,nEnt,localEnergy) * pow(redshift_to_inf[jR],3) * escapeAi[jR];
			lum_eeAnn += lumInterp(eeAnnVec,ntEnergies,jE,nEnt,localEnergy) * pow(redshift_to_inf[jR],3) * escapeAi[jR];
			lum_totAbs += lumInterp(totAbsVec,ntEnergies,jE,nEnt,localEnergy) * pow(redshift_to_inf[jR],3) * escapeAi[jR];
		}
		file << fmtE
			 << '\t' << safeLog10(lum_eSy*E)
			 << '\t' << safeLog10(lum_muSy*E)
			 << '\t' << safeLog10(lum_piSy*E)
			 << "\t" << safeLog10(lum_eIC*E)
			 << "\t" << safeLog10(lum_piPP*E)
			 << "\t" << safeLog10(lum_piPG*E)
			 << '\t' << safeLog10(lum_eeAnn*E)
			 << '\t' << safeLog10(lum_totAbs*E)
			 << std::endl;
	}
	
	file.close();
	show_message(msgEnd,Module_secondariesLuminosities);
}


void secondariesProcesses(State& st)
{
	/*
	show_message(msgStart,Module_secondariesTimescales);
	secondariesTimescales(st.ntChargedPion,st,"pionCoolingTimes.dat");
	secondariesTimescales(st.ntMuon,st,"muonCoolingTimes.dat");
	show_message(msgEnd,Module_secondariesTimescales);
	*/
	
	// PION INJECTION AND TRANSPORT
	injectionChargedPion(st.ntChargedPion,st);
	writeEandRParamSpace("pionInjection",st.ntChargedPion.injection,0,1);
	distributionSecondaries(st.ntChargedPion,st);
	writeEandRParamSpace("pionDistribution",st.ntChargedPion.distribution,0,1);
	
	// MUON INJECTION AND TRANSPORT
	injectionMuon(st.ntMuon,st);
	writeEandRParamSpace("pionInjection",st.ntMuon.injection,0,1);
	distributionSecondaries(st.ntMuon,st);
	writeEandRParamSpace("muonDistribution",st.ntMuon.distribution,0,1);
	
	int cond = 0;
	int it = 0;
	do {
		it++;
		injectionPair(st.ntPair,st,it);
		distributionSecondaries(st.ntPair,st);
		writeEandRParamSpace("secondaryPairDistribution",st.ntPair.distribution,0,1);
		secondariesRadiationProcesses(st,"secondariesLum.dat");
		cout << "Iteration = " << it << "\t Cond = " << cond << endl;
	} while (it<=1);
}