#include "secondariesProcesses.h"
#include "absorption.h"
#include "modelParameters.h"
#include "write.h"
#include "messages.h"
#include "globalVariables.h"
#include "adafFunctions.h"
#include "NTtimescales.h"

#include <fmath/RungeKutta.h>

#include <fluminosities/luminositySynchrotron.h>
#include <fluminosities/luminosityPhotoHadronic.h>
#include <fluminosities/luminosityNTHadronic.h>
#include <fluminosities/opticalDepthSSA.h>
#include <fluminosities/luminosityIC.h>
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
	Matrix eSyLocal;
	matrixInit(eSyLocal,nEnt,nR,0.0);
	
	Vector eSy(nEnt,0.0);
	Vector eIC(nEnt,0.0);
	Vector muSy(nEnt,0.0);
	Vector piSy(nEnt,0.0);
	Vector piPP(nEnt,0.0);
	Vector piPG(nEnt,0.0);
	Vector totAbs(nEnt,0.0);

	#pragma omp parallel for
	for (int E_ix=0;E_ix<nEnt;E_ix++) {
		double E = ntPh.ps[DIM_E][E_ix];
		st.ntPhoton.ps.iterate([&](const SpaceIterator& iR) {
			double r = iR.val(DIM_R);
			double rB1 = r/sqrt(paso_r);
			double rB2 = r*sqrt(paso_r);
			double vol = volume(r);
			double magf = st.magf.get(iR);
			SpaceCoord psc = {E_ix,iR.coord[DIM_R],0};
			
			double height = height_fun(r);
			double tau_gg = st.tau_gg.get(psc);
			double kappa_gg = 2.0/sqrt(pi)*tau_gg/height;
			
			double eICLocal,eSyLocal,piPPLocal,piPGLocal,piSyLocal,muSyLocal;
			eICLocal = eSyLocal = piPPLocal = piPGLocal = piSyLocal = muSyLocal = 0.0;
			eICLocal = luminosityIC_2(E,st.ntPair,iR.coord,st.photon.distribution,Ephmin,Ephmax)/E;
			eICLocal += luminosityIC(E,st.ntPair,iR.coord,st.ntPhoton.distribution,
										st.ntPhoton.emin(),st.ntPhoton.emax())/E;
			eSyLocal = luminositySynchrotronExactSec(E,st.ntPair,iR,magf)/E;
			muSyLocal = luminositySynchrotronExact(E,st.ntMuon,iR,magf)/E;
			piSyLocal = luminositySynchrotronExact(E,st.ntChargedPion,iR,magf)/E;
			
			if (E/EV_TO_ERG > 1.0e7) {
				piPGLocal = luminosityPhotoHadronic(E,st.ntChargedPion,st.photon.distribution,iR,Ephmin,Ephmax)/E;
				piPPLocal = luminosityNTHadronic(E,st.ntChargedPion,st.denf_i.get(iR),iR)/E;
			}

			double factor = pi/sqrt(3.0)*(rB2*rB2-rB1*rB1) / vol;
			double eTot = eSyLocal+eICLocal+muSyLocal+piSyLocal+piPPLocal+piPGLocal;

			eSy[E_ix] += eSyLocal*vol * pow(redshift_to_inf[iR.coord[DIM_R]],3);
			muSy[E_ix] += muSyLocal*vol * pow(redshift_to_inf[iR.coord[DIM_R]],3);
			piSy[E_ix] += piSyLocal*vol * pow(redshift_to_inf[iR.coord[DIM_R]],3);
			eIC[E_ix] += eICLocal*vol * pow(redshift_to_inf[iR.coord[DIM_R]],3);
			piPP[E_ix] += piPPLocal*vol * pow(redshift_to_inf[iR.coord[DIM_R]],3);
			piPG[E_ix] += piPGLocal*vol * pow(redshift_to_inf[iR.coord[DIM_R]],3);
			double totLocal = (tau_gg > 1.0e-10) ? factor*vol*eTot/kappa_gg*(1.0-exp(-2.0*sqrt(3.0)*tau_gg)) : 
								eTot*vol;
			double totLocalNA = eTot * pow(redshift_to_inf[iR.coord[DIM_R]],3);
			totAbs[E_ix] += totLocal  * pow(redshift_to_inf[iR.coord[DIM_R]],3);
			
			double tau_es = st.denf_e.get(iR)*thomson*height;
			double tescape = height/cLight * (1.0+tau_es);
			SpaceCoord iE = {E_ix,iR.coord[DIM_R],0};
			double rate_gg = kappa_gg * cLight;
			double NTphot = st.ntPhoton.injection.get(iE);
			st.ntPhoton.distribution.set(iE,NTphot + eTot/E*pow(rate_gg+pow(tescape,-1),-1));
		},{E_ix,-1,0});
	}
	
	for (size_t jE=0;jE<nEnt;jE++) {
		double E = st.ntPhoton.ps[DIM_E][jE];
		double fmtE = safeLog10(E/EV_TO_ERG);
		file << fmtE
			 << '\t' << safeLog10(eSy[jE]*E)
			 << '\t' << safeLog10(muSy[jE]*E)
			 << '\t' << safeLog10(piSy[jE]*E)
			 << "\t" << safeLog10(eIC[jE]*E)
			 << "\t" << safeLog10(piPP[jE]*E)
			 << "\t" << safeLog10(piPG[jE]*E)
			 << '\t' << safeLog10(totAbs[jE]*E)
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
	
	// NEUTRINO INJECTION
	injectionNeutrino(st.neutrino,st);
	
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