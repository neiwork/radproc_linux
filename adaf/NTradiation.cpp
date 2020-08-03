#include "NTradiation.h"
#include "modelParameters.h"
#include "write.h"
#include "messages.h"
#include "globalVariables.h"
#include "adafFunctions.h"
#include "absorption.h"
#include "thermalProcesses.h"
#include "thermalCompton.h"

#include <fmath/RungeKutta.h>

#include <fluminosities/luminositySynchrotron.h>
#include <fluminosities/opticalDepthSSA.h>
#include <fluminosities/luminosityIC.h>
#include <fluminosities/luminosityNTHadronic.h>
#include <fluminosities/luminosityHadronic.h>
#include <fluminosities/luminosityPhotoHadronic.h>
#include <fluminosities/thermalSync.h>
#include <fluminosities/blackBody.h>
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
			double vol = volume(r);
			double lumReachingShell = 0.0;
			double tcross = r*(sqrt(paso_r)-1.0/sqrt(paso_r))/cLight;
			double tescape = height_fun(r)/cLight;
			for (size_t jjR=0;jjR<nR;jjR++)
				lumReachingShell += ( (jjR == jR) ? lumNT[jE][jR]*tescape : 
									reachAA[jjR][jR]*pow(redshift[jjR][jR],2)*lumNT[jE][jjR]*tcross );
			st.photon.injection.set(itER,lumReachingShell/(vol*E*E)); // erg^⁻1 cm^-3
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


void nonThermalRadiation(State& st, const std::string& filename)
{
	show_message(msgStart, Module_luminosities);

	std::ofstream file, fileSy;
	file.open(filename.c_str(),std::ios::out);
	fileSy.open("lumSy.dat", std::ios::out);
	file << "log(E/eV)"
		 << '\t' << "pSyn"
		 << "\t" << "eIC"
		 << "\t" << "pIC"
		 << "\t" << "pPP"
		 << "\t" << "pPG"
		 << "\t" << "totAbs"
		 << std::endl;
	
	fileSy << "log(E/eV)" << '\t' << "eSyn" << std::endl;
	
	double Ephmin = st.photon.emin();
	double Ephmax = st.photon.emax();
	
	Particle &ntPh = st.ntPhoton;
	size_t nEnt = ntPh.ps[DIM_E].size();
	Vector ntEnergies(nEnt,0.0);
	Vector energies(nE,0.0);
	for (size_t jE=0;jE<nEnt;jE++) {
		energies[jE] = st.photon.ps[DIM_E][jE];
		ntEnergies[jE] = st.ntPhoton.ps[DIM_E][jE];
	}
	
	Matrix eSy;			matrixInit(eSy,nE,nR,0.0);
	Matrix eIC;			matrixInit(eIC,nEnt,nR,0.0);
	Matrix pSy;			matrixInit(pSy,nEnt,nR,0.0);
	Matrix pIC;			matrixInit(pIC,nEnt,nR,0.0);
	Matrix pPP;			matrixInit(pPP,nEnt,nR,0.0);
	Matrix pPG;			matrixInit(pPG,nEnt,nR,0.0);
	Matrix totNotAbs;	matrixInit(totNotAbs,nEnt,nR,0.0);
	Matrix totAbs;		matrixInit(totAbs,nEnt,nR,0.0);

	if (calculateNTelectrons) {
		#pragma omp parallel for
		for (int E_ix=0;E_ix<nE;E_ix++) {
			double E = st.photon.ps[DIM_E][E_ix];
			st.photon.ps.iterate([&](const SpaceIterator &iR) {
				double r = iR.val(DIM_R);
				double rB1 = r/sqrt(paso_r);
				double rB2 = rB1*paso_r;
				double vol = volume(r);
				double magf = st.magf.get(iR);
				SpaceCoord psc = {E_ix,iR.coord[DIM_R],0};
				
				// FOR THE OBSERVED LUMINOSITY
				double alpha_th = jSync(E,st.tempElectrons.get(iR),magf,st.denf_e.get(iR)) / 
									bb(E/planck,st.tempElectrons.get(iR));
				double alpha_ssa = ssaAbsorptionCoeff(E,magf,st.ntElectron,psc);
				double alpha_tot = alpha_th + alpha_ssa;
				double height = height_fun(r);
				double tau = 0.5*sqrt(pi)*alpha_tot*height;

				double eSyLocal = luminositySynchrotronExact(E,st.ntElectron,iR,magf)/E;
				
				double factor = 0.5/sqrt(3.0);
				eSyLocal = (tau > 1.0e-10) ?
								0.5/sqrt(3.0) * (eSyLocal/alpha_tot) * 
								(1.0-exp(-2.0*sqrt(3.0)*tau)) :
								eSyLocal * 0.5 * sqrt(pi) * height;  // flux
				eSy[E_ix][iR.coord[DIM_R]]  = eSyLocal * 2.0 * pi*(rB2*rB2-rB1*rB1);  // luminosity
								
				
				double tau_es = st.denf_e.get(iR)*thomson*height;
				double tescape = height/cLight * (1.0+tau_es);
				SpaceCoord iE = {E_ix,iR.coord[DIM_R],0};
				st.photon.distribution.set(iE, st.photon.distribution.get(iE) + 
											eSyLocal / E / height * tescape);
			},{E_ix,-1,0});
		}
	}

	if (calculateNonThermalHE) {
		#pragma omp parallel for
		for (int E_ix=0;E_ix<nEnt;E_ix++) {
			double E = ntPh.ps[DIM_E][E_ix];
			size_t jR = 0;
			st.ntPhoton.ps.iterate([&](const SpaceIterator& iR) {
				double r = iR.val(DIM_R);
				double rB1 = r/sqrt(paso_r);
				double rB2 = r*sqrt(paso_r);
				double vol = volume(r);
				double magf = st.magf.get(iR);
				SpaceCoord psc = {E_ix,iR.coord[DIM_R],0};
				
				double alpha_th = jSync(E,st.tempElectrons.get(iR),magf,st.denf_e.get(iR)) / 
									bb(E/planck,st.tempElectrons.get(iR));
				double kappa_ssa_p = ssaAbsorptionCoeff(E,magf,st.ntProton,psc)+alpha_th;
				double kappagg = 0.0;
				double Ephmin_gg = P2(electronRestEnergy)/E;
				if (Ephmin_gg < Ephmax) {
					kappagg = integSimpsonLog(Ephmin_gg,Ephmax,[&E,&psc,&st](double Eph)
									{
										double nPh = st.photon.distribution.interpolate({{0,Eph}},&psc);
										return ggCrossSection2(E,Eph)*nPh;
									},50);
				}
				double height = height_fun(r);
				double tau_ssa_p = 0.5*sqrt(pi)*kappa_ssa_p*height;
				double taugg = 0.5*sqrt(pi)*kappagg*height;
				st.tau_gg.set(psc,taugg);
				
				double factor = pi/sqrt(3.0)*(rB2*rB2-rB1*rB1) / vol;
				
				double eICLocal, pSyLocal, pICLocal, pPPLocal, pPGLocal;
				eICLocal = pSyLocal = pICLocal = pPPLocal = pPGLocal = 0.0;
				
				if (calculateNTelectrons)
					eICLocal = luminosityIC_2(E,st.ntElectron,iR.coord,st.photon.distribution,Ephmin,Ephmax)/E;
				
				if (calculateNTprotons) {
					pSyLocal = luminositySynchrotronExact(E,st.ntProton,iR,st.magf.get(iR))/E;
					pICLocal = luminosityIC_Th(E,st.ntProton,iR,st.photon.distribution,Ephmin,Ephmax)/E;
				
					if (E/EV_TO_ERG > 1.0e7) {
						pPGLocal = luminosityPhotoHadronic(E,st.ntProton,st.photon.distribution,iR,Ephmin,Ephmax)/E;
						pPPLocal = luminosityNTHadronic(E,st.ntProton,st.denf_i.get(iR),iR)/E;
						if (E/EV_TO_ERG < 1.0e11)
							pPPLocal += 4.0*pi*luminosityHadronic(E,st.denf_i.get(iR),st.tempIons.get(iR))/planck;
					}
					
					pSyLocal = (tau_ssa_p > 1.0e-10) ? 
								factor*(pSyLocal/kappa_ssa_p)*(1.0-exp(-2.0*sqrt(3.0)*tau_ssa_p)) :
								pSyLocal;
				}
				
				double eTot = eICLocal+pSyLocal+pICLocal+pPPLocal+pPGLocal;

				pSy[E_ix][jR] = pSyLocal*vol;
				pIC[E_ix][jR] = pICLocal*vol;
				eIC[E_ix][jR] = eICLocal*vol;
				pPP[E_ix][jR] = pPPLocal*vol;
				pPG[E_ix][jR] = pPGLocal*vol;
				
				double totLocal = (taugg > 1.0e-10) ? 
										factor*vol*eTot/kappagg*(1.0-exp(-2.0*sqrt(3.0)*taugg)) : 
											eTot*vol;
				//totLocal = eTot*vol*exp(-taugg);
				totNotAbs[E_ix][jR] = totLocal;
				totAbs[E_ix][jR] = (eTot*vol - totLocal);
				
				
				double tau_es = st.denf_e.get(iR)*thomson*height;
				double tescape = height/cLight * (1.0+tau_es);
				double rate_gg = kappagg * cLight;
				st.ntPhoton.distribution.set(psc, eTot/E * pow(rate_gg+pow(tescape,-1),-1));
				st.ntPhoton.injection.set(psc, eTot/E * pow(rate_gg+pow(tescape,-1),-1));
				
				jR++;
			},{E_ix,-1,0});
		}
	}
	
	double lum_eSy, lum_pSy, lum_eIC, lum_pIC, lum_pPP, lum_pPG, lum_totNotAbs, lum_totAbs;
	
	if (calculateNonThermalHE) {
		for (size_t jE=0;jE<nEnt;jE++) {
			lum_pSy = lum_eIC = lum_pIC = lum_pPP = lum_pPG = lum_totNotAbs = lum_totAbs = 0.0;
			double E = ntEnergies[jE];
			double fmtE = safeLog10(E/EV_TO_ERG);
			for (size_t jR=0;jR<nR;jR++) {
				Vector pSyVec(nEnt,0.0), eICVec(nEnt,0.0), pICVec(nEnt,0.0), pPPVec(nEnt,0.0),
						pPGVec(nEnt,0.0), totNotAbsVec(nEnt,0.0), totAbsVec(nEnt,0.0);
				for (size_t jjE=0;jjE<nEnt;jjE++) {
					pSyVec[jjE] = pSy[jjE][jR];
					eICVec[jjE] = eIC[jjE][jR];
					pICVec[jjE] = pIC[jjE][jR];
					pPPVec[jjE] = pPP[jjE][jR];
					pPGVec[jjE] = pPG[jjE][jR];
					totNotAbsVec[jjE] = totNotAbs[jjE][jR];
					totAbsVec[jjE] = totAbs[jjE][jR];
				}
				double localEnergy = E/redshift_to_inf[jR];
				lum_pSy += lumInterp(pSyVec,ntEnergies,jE,nEnt,localEnergy) * pow(redshift_to_inf[jR],3) * escapeAi[jR];
				lum_eIC += lumInterp(eICVec,ntEnergies,jE,nEnt,localEnergy) * pow(redshift_to_inf[jR],3) * escapeAi[jR];
				lum_pIC += lumInterp(pICVec,ntEnergies,jE,nEnt,localEnergy) * pow(redshift_to_inf[jR],3) * escapeAi[jR];
				lum_pPP += lumInterp(pPPVec,ntEnergies,jE,nEnt,localEnergy) * pow(redshift_to_inf[jR],3) * escapeAi[jR];
				lum_pPG += lumInterp(pPGVec,ntEnergies,jE,nEnt,localEnergy) * pow(redshift_to_inf[jR],3) * escapeAi[jR];
				lum_totNotAbs += lumInterp(totNotAbsVec,ntEnergies,jE,nEnt,localEnergy) * pow(redshift_to_inf[jR],3) * escapeAi[jR];
				lum_totAbs += lumInterp(totAbsVec,ntEnergies,jE,nEnt,localEnergy) * pow(redshift_to_inf[jR],3) * escapeAi[jR];
			}
			file << fmtE
				 << '\t' << safeLog10(lum_pSy*E)
				 << "\t" << safeLog10(lum_eIC*E)
				 << "\t" << safeLog10(lum_pIC*E)
				 << "\t" << safeLog10(lum_pPP*E)
				 << "\t" << safeLog10(lum_pPG*E)
				 << '\t' << safeLog10(lum_totNotAbs*E)
				 << '\t' << safeLog10(lum_totAbs*E)
				 << std::endl;
		}
	}
	
	ofstream file2;
	file2.open("lumSy_radius.dat");
	for (size_t jE=0;jE<nE;jE++) {
		double E = energies[jE];
		double fmtE = safeLog10(E/EV_TO_ERG);
		lum_eSy = 0.0;
		for (size_t jR=0;jR<nR;jR++) {
			double r = st.denf_e.ps[DIM_R][jR];
			Vector eSyVec(nE,0.0);
			for (size_t jjE=0;jjE<nE;jjE++)
				eSyVec[jjE] = eSy[jjE][jR];
			double localEnergy = E/redshift_to_inf[jR];
			double lumR = lumInterp(eSyVec,energies,jE,nE,localEnergy) * pow(redshift_to_inf[jR],3) * escapeAi[jR];
			lum_eSy += lumR;
			file2 << E/EV_TO_ERG << "\t" << r*sqrt(paso_r)/schwRadius << "\t" << lum_eSy << endl;
		}
		
		
		fileSy << fmtE
			   << '\t' << safeLog10(lum_eSy*E)
			   << std::endl;
	}
	
	file.close();
	fileSy.close();
	file2.close();
	
	ofstream file_r;
	file_r.open("lumNonThermal_r.dat");
	double pasoE = pow(ntEnergies[nEnt-1]/ntEnergies[0],1.0/(nEnt-1));
	double pasoE2 = pow(energies[nEnt-1]/energies[0],1.0/(nE-1));
	double lum_eSy_r, lum_eIC_r, lum_pSy_r, lum_pp_r, lum_pg_r, lum_tot_r;
	for (size_t jR=0;jR<nR;jR++) {
		lum_eSy_r = lum_eIC_r = lum_pSy_r = lum_pp_r = lum_pg_r = lum_tot_r = 0.0;
		double r = st.denf_e.ps[DIM_R][jR];
		for (size_t jE=0;jE<nEnt;jE++) {
			double E = ntEnergies[jE];
			double E_eV = E / EV_TO_ERG;
			double dE = E * (sqrt(pasoE)-1.0/sqrt(pasoE));
			lum_eIC_r += eIC[jE][jR] * dE;
			lum_pSy_r += pSy[jE][jR] * dE;
			lum_pp_r += pPP[jE][jR] * dE;
			lum_pg_r += pPG[jE][jR] * dE;
			lum_tot_r += totNotAbs[jE][jR] * dE;
		}
		for (size_t jE=0;jE<nE;jE++) {
			double E = energies[jE];
			double E_eV = E / EV_TO_ERG;
			double dE = E * (sqrt(pasoE2)-1.0/sqrt(pasoE2));
			lum_eSy_r += eSy[jE][jR] * dE;
		}
		file_r	<< r / schwRadius
				<< '\t' << lum_eSy_r * escapeAi[jR] * pow(redshift_to_inf[jR],4)
				<< '\t' << lum_eIC_r * escapeAi[jR] * pow(redshift_to_inf[jR],4)
				<< '\t' << lum_pSy_r * escapeAi[jR] * pow(redshift_to_inf[jR],4)
				<< '\t' << lum_pp_r * escapeAi[jR] * pow(redshift_to_inf[jR],4)
				<< '\t' << lum_pg_r * escapeAi[jR] * pow(redshift_to_inf[jR],4)
				<< '\t' << lum_tot_r * escapeAi[jR] * pow(redshift_to_inf[jR],4)
				<< std::endl;
	}
	file_r.close();
	
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
