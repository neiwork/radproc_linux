#include "processes.h"
#include "modelParameters.h"
#include "write.h"
#include "messages.h"
#include "globalVariables.h"
#include "adafFunctions.h"
#include "absorption.h"
#include "thermalProcesses.h"

#include <fmath/RungeKutta.h>

#include <fluminosities/luminositySynchrotron.h>
#include <fluminosities/opticalDepthSSA.h>
#include <fluminosities/luminosityIC.h>
#include <fluminosities/luminosityNTHadronic.h>
#include <fluminosities/luminosityPhotoHadronic.h>

#include <fparameters/Dimension.h>
#include <fparameters/parameters.h>
#include <fparameters/SpaceIterator.h>


void targetFieldNT(State& st, Matrix lumNT)
{
	size_t jE=0;
	st.photon.ps.iterate([&](const SpaceIterator& itE) {
		size_t jR=0;
		double E = itE.val(DIM_E);
		st.photon.ps.iterate([&](const SpaceIterator& itER) {
			double r = itER.val(DIM_R);
			double thetaH = st.thetaH.get(itER);
			double area = 4.0*pi*r*r*cos(thetaH);
			double lumReachingShell = 0.0;
			for (size_t jjR=0;jjR<nR;jjR++)
				lumReachingShell += reachAA[jjR][jR]*lumNT[jE][jjR];
			st.photon.injection.set(itER,lumReachingShell/(area*cLight*E*E)); // erg^⁻1 cm^-3
			jR++;
		},{itE.coord[DIM_E],-1,0});
		jE++;
	},{-1,0,0});
}

/*
double opticalDepthSSA(int E_ix, State& st, const Particle& p, double r_current)  
{
	double E = st.photon.ps[DIM_E][E_ix];  //E=Eph
	double mass = p.mass;
	double opacity = 0.0;
	st.photon.ps.iterate([&](const SpaceIterator &i) {
		
		double r = i.val(DIM_R);
		if (r >= r_current){ //ver si la int es de 0 a r_current o de r_current hasta que escape
		
			double rB1=r/sqrt(paso_r);
			double rB2=r*sqrt(paso_r);
			double delta_r = rB2-rB1;
		
			double magneticField = st.magf.get(i);
				
			double integral = intSimple(p.emin(),p.emax(),[&](double x) {
								return fSSA(x,E,p,magneticField,i); });

			double absorptionCoefficient = - P3(planck)*cLight2*integral/(8*pi*P2(E)); 
		
			opacity += absorptionCoefficient*delta_r; // parameters.radius;
		}
		
	},{E_ix,-1,0});
		
		return opacity;
}

double opticalDepthSSA2(int E_ix, State& st, const Particle& p, int iR)
{
	double E = st.photon.ps[DIM_E][E_ix];
	double rMax = st.photon.ps[DIM_R][nR-1];
	double rMin = st.photon.ps[DIM_R][0];
	double pasoprim = pow(rMax/rMin,1.0/(nR*2.0));

	SpaceCoord psc = {E_ix,iR,0};
	double r0 = st.photon.ps[DIM_R][iR];
	double drprim = r0*(pasoprim-1.0);
	double thetaprim = inclination*(pi/180.0);

	size_t nPhi = 5;
	double opticalDepth = 0.0;
	double theta0 = 0.5*pi;
	double phi0 = 0.0;
	double dPhi = 2.0*pi/nPhi;
	for (size_t jPhi=0;jPhi<nPhi;jPhi++) {
		double x0=r0*sin(theta0)*cos(phi0);
		double y0=r0*sin(theta0)*sin(phi0);
		double z0=r0*cos(theta0);
		
		double rprim = drprim;
		double r1 = r0;
		double theta1 = theta0;
		double thetaMinLocal = theta0;
		while (r1 < rMax && r1 > rMin && theta1 > thetaMinLocal && theta1 < pi-thetaMinLocal) {
			double xprim = rprim*sin(thetaprim);
			double zprim = rprim*cos(thetaprim);
			double x1 = x0+xprim;
			double z1 = z0+zprim;
			r1=sqrt(x1*x1+y0*y0+z1*z1);
			
			double cosThetaMinLocal = costhetaH(r1);
			theta1 = atan(sqrt(x1*x1+y0*y0)/abs(z1));
			thetaMinLocal = acos(cosThetaMinLocal);
			
			double magneticField = (r1 < rMax && r1 > rMin && 
							theta1 > thetaMinLocal && theta1 < pi-thetaMinLocal) ? 
							st.magf.interpolate({{1,r1}},&psc) : 0.0;
			double integral = intSimple(p.emin(),p.emax(),[&](double x){
							return fSSA2(x,E,p,magneticField,psc);});
			
			double absorptionCoefficient = -P3(planck)*cLight2*integral/(8*pi*P2(E)); 
			opticalDepth += absorptionCoefficient*drprim;
			
			drprim = r1*(pasoprim-1.0);
			rprim += drprim;
		}
		phi0 += dPhi;
	}
	return opticalDepth/nPhi;
}
*/

void processes(State& st, const std::string& filename)
{
	show_message(msgStart, Module_luminosities);

	std::ofstream file, file2;
	file.open(filename.c_str(),std::ios::out);
	file2.open("opticalDepth.txt", std::ios::out);
	file << "log(E/eV)"
		 << '\t' << "eSyn"
		 << '\t' << "pSyn"
		 << "\t" << "eIC"
		 << "\t" << "pPP"
		 << "\t" << "pPG"
		 << "\t" << "eAbs"
		 << "\t" << "pAbs"
		 << std::endl;
			
	/*file2 << "Log(E/eV)" 
		  << "\t" << "r [M]"
		  << "\t" << "tau_e"
		  << "\t" << "tau_p"
		  << "\t" << "tau_gg"
		  << std::endl;*/
	
	double Emin = st.photon.emin();
	double Emax = st.photon.emax();
	int nE = st.photon.ps[DIM_E].size();
	
	Vector eSy(nE,0.0);
	Vector eIC(nE,0.0);
	Vector pSy(nE,0.0);
	Vector pPP(nE,0.0);
	Vector pPG(nE,0.0);
	Vector eAbs(nE,0.0);
	Vector pAbs(nE,0.0);
	
	Matrix lumNT_ssa;	matrixInit(lumNT_ssa,nE,nR,0.0);
	
	Matrix tau_e;	matrixInit(tau_e,nE,nR,0.0);
	Matrix tau_p;	matrixInit(tau_p,nE,nR,0.0);
	Matrix tau_gg;	matrixInit(tau_gg,nE,nR,0.0);
	
	Vector energies(nE,0.0);
	double lumLepSy = 0.0;
	double lumLepIC = 0.0;
	double lumHad = 0.0;
	double pasoE = pow(st.photon.emax()/st.photon.emin(),1.0/nE);

	#pragma omp parallel for
	for (int E_ix=0;E_ix<nE;E_ix++) {
		double E = st.photon.ps[DIM_E][E_ix];
		energies[E_ix] = E;
		double lumLepSyLocal = 0.0;
		double lumLepICLocal = 0.0;
		double lumHadLocal = 0.0;
		st.photon.ps.iterate([&](const SpaceIterator &i) {
			double r = i.val(DIM_R);
			double rB1 = r/sqrt(paso_r);
			double rB2 = r*sqrt(paso_r);
			double vol = (4.0/3.0)*pi*cos(st.thetaH.get(i))*(rB2*rB2*rB2-rB1*rB1*rB1);
			
			Vector tau(3,0.0);
			double fmtE  = safeLog10(i.val(DIM_E)/1.6e-12);
			if (fmtE < 0.0 || fmtE > 5.0) opticalDepth(tau,E_ix,st,i.coord[DIM_R]);

			double attenuation_gg = exp(-tau[0]);
			double attenuation_ssae = (tau[1] > 1.0e-15) ? (1.0-exp(-tau[1]))/tau[1] : 1.0;
			double attenuation_ssap = (tau[2] > 1.0e-15) ? (1.0-exp(-tau[2]))/tau[2] : 1.0;
			
			double eSyLocal,eICLocal,pSyLocal,pPPLocal,pPGLocal;
			eSyLocal = eICLocal = pSyLocal = pPPLocal = pPGLocal = 0.0;
			
			//////////////////////////////////////////////////////
			/*double r0 = r;
			double dr = r/100;
			double tAux = 0.0;
			while (tAux < tAccBlob) {
				r0 -= dr;
				tAux += dr/(-radialVel(r0));
			}
			double B = (r0 > st.magf.ps[DIM_R][0]) ? 
					sqrt( (1.0-magFieldPar)*8.0*pi*massDensityADAF(r0)*sqrdSoundVel(r0) ) : 0.0;
			eSyLocal = luminositySynchrotron2(E,st.ntElectron,i,B);*/
			//////////////////////////////////////////////////////
			
			eSyLocal = luminositySynchrotron(E,st.ntElectron,i,st.magf);
			//eICLocal = luminosityIC(E,st.ntElectron,i.coord,st.photon.distribution,Emin);
			//pSyLocal = luminositySynchrotron(E,st.ntProton,i,st.magf);
			//pPPLocal = luminosityNTHadronic(E,st.ntProton,st.denf_i.get(i),i);
			//pPGLocal = luminosityPhotoHadronic(E,st.ntProton,st.photon.distribution,i,Emin,Emax);

			eSy[E_ix] += eSyLocal*vol;				// [erg s^-1]
			eIC[E_ix] += eICLocal*vol;
			pSy[E_ix] += pSyLocal*vol;
			pPP[E_ix] += pPPLocal*vol;
			pPG[E_ix] += pPGLocal*vol;
			eAbs[E_ix] += (eSyLocal+eICLocal)*vol*attenuation_ssae*attenuation_gg;
			pAbs[E_ix] += (pSyLocal+pPPLocal+pPGLocal)*vol*attenuation_ssap*attenuation_gg;
			tau_gg[E_ix][i.coord[DIM_R]] = tau[0];
			tau_e[E_ix][i.coord[DIM_R]] = tau[1];
			tau_p[E_ix][i.coord[DIM_R]] = tau[2];
			
			lumNT_ssa[E_ix][i.coord[DIM_R]] = ( (eSyLocal+eICLocal)*attenuation_ssae + 
												(pSyLocal+pPPLocal+pPGLocal)*attenuation_ssap ) * vol;
			lumLepSyLocal += eSyLocal*attenuation_ssae*vol;
			lumLepICLocal += eICLocal*attenuation_ssae*vol;
			lumHadLocal += (pSyLocal + pPPLocal + pPGLocal)*attenuation_ssap*vol;
		},{E_ix,-1,0});
		lumLepSy += lumLepSyLocal*(sqrt(pasoE)-1.0/sqrt(pasoE));
		lumLepIC += lumLepICLocal*(sqrt(pasoE)-1.0/sqrt(pasoE));
		lumHad += lumHadLocal*(sqrt(pasoE)-1.0/sqrt(pasoE));
	}
	
	cout << "Leptonic Sync luminosity = " << lumLepSy << endl;
	cout << "Leptonic IC luminosity = " << lumLepIC << endl;
	cout << "Hadronic luminosity = " << lumHad << endl;
	targetFieldNT(st,lumNT_ssa);
	writeEandRParamSpace("nonThermalPhotonDensity.txt",st.photon.injection,0);
	writeEandRParamSpace("thermalPhotonDensity.txt",st.photon.distribution,0);
	
	for (size_t jR=0;jR<nR;jR++) {
		double r = st.photon.ps[DIM_R][jR] / schwRadius;
		for (size_t jE=0;jE<nE;jE++) {
			double fmtE = safeLog10(st.photon.ps[DIM_E][jE]/1.6e-12);
			file2 << fmtE << "\t" << safeLog10(r)
						  << "\t" << tau_e[jE][jR]
						  << "\t" << tau_p[jE][jR]
						  << "\t" << tau_gg[jE][jR]
						  << std::endl;
		}
	}
	
	for (size_t jE=0;jE<nE;jE++) {
		double fmtE = safeLog10(st.photon.ps[DIM_E][jE]/1.6e-12);
		file << fmtE
			 << '\t' << safeLog10(eSy[jE])
			 << '\t' << safeLog10(pSy[jE])
			 << "\t" << safeLog10(eIC[jE])
			 << "\t" << safeLog10(pPP[jE])
			 << "\t" << safeLog10(pPG[jE])
			 << '\t' << safeLog10(eAbs[jE])
			 << '\t' << safeLog10(pAbs[jE])
			 << std::endl;
	}
	
	if (calculateFlare) {
		//double timeBurst = GlobalConfig.get<double>("nonThermal.flare.timeAfterFlare");
		writeBlob(st,energies,lumNT_ssa,timeAfterFlare);
	}

	file.close();
	file2.close();

	show_message(msgEnd,Module_luminosities);
}

///////////////////
//#   pragma omp parallel for \
	//		private(i, eSyn, eIC) \
	//		shared(st, Qsyn, Qic) \
	//		default(none) \
	//		schedule(static, 1) \
	//		num_threads(2)

//#pragma omp parallel sections
//{
//#pragma omp section
//	{
//		Qsyn.fill([&st](const SpaceIterator &i){
//			return luminositySynchrotron(i.val(DIM_E), st.electron); //estos devuelven erg/s/cm^3, integrar!
//		});
//	}

//#pragma omp section
//	{
//		Qic.fill([&st](const SpaceIterator &i){
//			return luminosityAnisotropicIC(i.val(DIM_E), st.electron, i.val(DIM_R));
//		});
//	}
//}

