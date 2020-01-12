#include "secondariesProcesses.h"
#include "modelParameters.h"
#include "write.h"
#include "messages.h"
#include "globalVariables.h"
#include "adafFunctions.h"
//#include "absorption.h"

#include <fmath/RungeKutta.h>

#include <fluminosities/luminositySynchrotron.h>
#include <fluminosities/opticalDepthSSA.h>
#include <fluminosities/luminosityIC.h>
#include "NTinjection.h"
#include "NTdistribution.h"
#include <fparameters/Dimension.h>
#include <fparameters/SpaceIterator.h>

void secondariesRadiationProcesses(State& st, const std::string& filename)
{
	show_message(msgStart, Module_luminosities);

	std::ofstream file;
	file.open(filename.c_str(),std::ios::out);
	std::ifstream file2;
	file2.open("opticalDepth.txt", std::ios::in);
	file << "log(E/eV)"
		 << '\t' << "eSyn"
		 << "\t" << "eIC"
		 << "\t" << "eAbs"
		 << std::endl;
		
		
	double Emin = st.photon.emin();
	double Emax = st.photon.emax();
	
	Matrix tau_e;	matrixInit(tau_e,nE,nR,0.0);
	//Matrix tau_p;	matrixInit(tau_p,nE,nR,0.0);
	Matrix tau_gg;	matrixInit(tau_gg,nE,nR,0.0);
	
	for (size_t jR=0;jR<nR;jR++) {
		double logr; // = st.photon.ps[DIM_R][jR] / schwRadius;
		for (size_t jE=0;jE<nE;jE++) {
			double fmtE; // = safeLog10(st.photon.ps[DIM_E][jE]/1.6e-12);
			double tau_p;
			file2 >> fmtE >> logr
						  >> tau_e[jE][jR]
						  >> tau_p //tau_p[jE][jR]
						  >> tau_gg[jE][jR];  
		}
	}
	
	
	Vector eSy(nE,0.0);
	Vector eIC(nE,0.0);
	Vector eAbs(nE,0.0);
	
	Matrix lumNT_ssa;	matrixInit(lumNT_ssa,nE,nR,0.0);
	
	

	#pragma omp parallel for
	for (int E_ix=0;E_ix<nE;E_ix++) {
		double E = st.photon.distribution.ps[DIM_E][E_ix];
		st.photon.ps.iterate([&](const SpaceIterator &i) {
			double r = i.val(DIM_R);
			double rB1 = r/sqrt(paso_r);
			double rB2 = r*sqrt(paso_r);
			double vol = (4.0/3.0)*pi*cos(st.thetaH.get(i))*(rB2*rB2*rB2-rB1*rB1*rB1);
			
			Vector tau(3,0.0);
			double fmtE  = safeLog10(i.val(DIM_E)/1.6e-12);
			tau[1] = tau_e[E_ix][i.coord[DIM_R]];
			tau[0] = tau_gg[E_ix][i.coord[DIM_R]];
			
			double attenuation_gg = exp(-tau[0]);
			double attenuation_ssae = (tau[1] > 1.0e-15) ? (1.0-exp(-tau[1]))/tau[1] : 1.0;
			//double attenuation_ssap = (tau[2] > 1.0e-15) ? (1.0-exp(-tau[2]))/tau[2] : 1.0;
			
			double eSyLocal = luminositySynchrotron(E,st.ntPair,i,st.magf);
			double eICLocal = luminosityIC(E,st.ntPair,i.coord,st.photon.distribution,Emin,Emax);
			

			eSy[E_ix] += eSyLocal*vol;
			eIC[E_ix] += eICLocal*vol;
			eAbs[E_ix] += (eSyLocal+eICLocal)*vol*attenuation_ssae*attenuation_gg;
					
			
			lumNT_ssa[E_ix][i.coord[DIM_R]] = ( (eSyLocal+eICLocal)*attenuation_ssae) * vol;

		},{E_ix,-1,0});
	}
	
	//targetFieldNT(st,lumNT_ssa);
	
	for (size_t jE=0;jE<nE;jE++) {
		double fmtE = safeLog10(st.photon.ps[DIM_E][jE]/1.6e-12);
		file << fmtE
			 << '\t' << safeLog10(eSy[jE])
			 << "\t" << safeLog10(eIC[jE])
			 << '\t' << safeLog10(eAbs[jE])
			 << std::endl;
	}
	
	file.close();
	file2.close();

	show_message(msgEnd,Module_luminosities);
}


void secondariesProcesses(State& st)
{
	injectionChargedPion(st.ntChargedPion,st);
	//distributionFast(st.ntChargedPion,st);
	injectionMuon(st.ntMuon,st);
	//distributionFast(st.ntMuon,st);
	injectionPair(st.ntPair,st);
	distributionFast(st.ntPair,st);
	//secondariesRadiationProcesses(st,"secondariesRadiation");
}