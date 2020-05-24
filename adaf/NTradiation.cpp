#include "NTradiation.h"
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
		 << "\t" << "pPP"
		 << "\t" << "pPG"
		 << "\t" << "totAbs"
		 << std::endl;
	
	fileSy << "log(E/eV)" << '\t' << "eSyn" << std::endl;
	
	double Ephmin = st.photon.emin();
	double Ephmax = st.photon.emax();
	
	Particle &ntPh = st.ntPhoton;
	size_t nEnt = ntPh.ps[DIM_E].size();
	Matrix eSyLocal;
	matrixInit(eSyLocal,nEnt,nR,0.0);
	
	Vector eSy(nEnt,0.0);
	Vector eIC(nEnt,0.0);
	Vector pSy(nEnt,0.0);
	Vector pPP(nEnt,0.0);
	Vector pPG(nEnt,0.0);
	Vector totNotAbs(nEnt,0.0);
	Vector totAbs(nEnt,0.0);

	if (calculateNTelectrons) {
		#pragma omp parallel for
		for (int E_ix=0;E_ix<nE;E_ix++) {
			double E = st.photon.ps[DIM_E][E_ix];
			st.photon.ps.iterate([&](const SpaceIterator &iR) {
				double r = iR.val(DIM_R);
				E = E / redshift_to_inf[iR.coord[DIM_R]];
				double rB1 = r/sqrt(paso_r);
				double rB2 = rB1*paso_r;
				double vol = volume(r);
				double magf = st.magf.get(iR);
				SpaceCoord psc = {E_ix,iR.coord[DIM_R],0};
				
				double alpha_th = jSync(E,st.tempElectrons.get(iR),magf,st.denf_e.get(iR)) / 
									bb(E/planck,st.tempElectrons.get(iR));
				double alpha_ssa = ssaAbsorptionCoeff(E,magf,st.ntElectron,psc);
				double alpha_tot = alpha_th + alpha_ssa;
				double height = height_fun(r);
				double tau = 0.5*sqrt(pi)*alpha_tot*height;

				eSyLocal[E_ix][iR.coord[DIM_R]] = luminositySynchrotronExact(E,st.ntElectron,iR,magf)/E;
				
				double factor = 0.5/sqrt(3.0);
				eSyLocal[E_ix][iR.coord[DIM_R]] = (tau > 1.0e-10) ?
								0.5/sqrt(3.0) * (eSyLocal[E_ix][iR.coord[DIM_R]]/alpha_tot) * 
								(1.0-exp(-2.0*sqrt(3.0)*tau)) :
								eSyLocal[E_ix][iR.coord[DIM_R]] * 0.5 * sqrt(pi) * height;  // flux
				eSy[E_ix] += eSyLocal[E_ix][iR.coord[DIM_R]] * 2.0 * pi*(rB2*rB2-rB1*rB1);  // luminosity
				double tau_es = st.denf_e.get(iR)*thomson*height;
				double tescape = height/cLight * (1.0+tau_es);
				SpaceCoord iE = {E_ix,iR.coord[DIM_R],0};
				st.photon.distribution.set(iE, st.photon.distribution.get(iE) + 
											eSyLocal[E_ix][iR.coord[DIM_R]] / E / height * tescape);
			},{E_ix,-1,0});
		}
	}

	/*
	// Adding nonthermal photons to the thermal ones; specially Sync photons.
	st.photon.ps.iterate([&](const SpaceIterator& i) {
		double e = i.val(DIM_E);
		double nPhNT = st.ntPhoton.distribution.interpolate({{DIM_E,e}},&i.coord);
		st.photon.distribution.set(i,st.photon.distribution.get(i)+nPhNT);
	},{-1,-1,0});
*/
	#pragma omp parallel for
	for (int E_ix=0;E_ix<nEnt;E_ix++) {
		double E = ntPh.ps[DIM_E][E_ix];
		st.ntPhoton.ps.iterate([&](const SpaceIterator& iR) {
			E = E / redshift_to_inf[iR.coord[DIM_R]];
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
			
			double eICLocal,pSyLocal,pPPLocal,pPGLocal;
			eICLocal = pSyLocal = pPPLocal = pPGLocal = 0.0;
			
			if (calculateNTelectrons)
				eICLocal = luminosityIC(E,st.ntElectron,iR.coord,st.photon.distribution,Ephmin,Ephmax)/E;
				
			if (calculateNTprotons) {
				pSyLocal = luminositySynchrotronExact(E,st.ntProton,iR,st.magf.get(iR))/E;
			
				if (E/EV_TO_ERG > 1.0e7) {
					pPGLocal = luminosityPhotoHadronic(E,st.ntProton,st.photon.distribution,iR,Ephmin,Ephmax)/E;
					pPPLocal = luminosityNTHadronic(E,st.ntProton,st.denf_i.get(iR),iR)/E;
					if (E/EV_TO_ERG < 1.0e10){}
						pPPLocal += 4.0*pi*luminosityHadronic(E,st.denf_i.get(iR),st.tempIons.get(iR))/planck;
				}
				
				pSyLocal = (tau_ssa_p > 1.0e-10) ? 
							factor*(pSyLocal/kappa_ssa_p)*(1.0-exp(-2.0*sqrt(3.0)*tau_ssa_p)) :
							pSyLocal;
			}
			
			double eTot = eICLocal+pSyLocal+pPPLocal+pPGLocal;

			pSy[E_ix] += pSyLocal*vol * pow(redshift_to_inf[iR.coord[DIM_R]],3);
			eIC[E_ix] += eICLocal*vol * pow(redshift_to_inf[iR.coord[DIM_R]],3);
			pPP[E_ix] += pPPLocal*vol * pow(redshift_to_inf[iR.coord[DIM_R]],3);
			pPG[E_ix] += pPGLocal*vol * pow(redshift_to_inf[iR.coord[DIM_R]],3);
			double totLocal = (taugg > 1.0e-10) ? 
									factor*vol*eTot/kappagg*(1.0-exp(-2.0*sqrt(3.0)*taugg)) : 
										eTot*vol;
			totNotAbs[E_ix] += totLocal * pow(redshift_to_inf[iR.coord[DIM_R]],3);
			totAbs[E_ix] += (eTot*vol - totLocal) * pow(redshift_to_inf[iR.coord[DIM_R]],3);
			
			double tau_es = st.denf_e.get(iR)*thomson*height;
			double tescape = height/cLight * (1.0+tau_es);
			double rate_gg = kappagg * cLight;
			st.ntPhoton.distribution.set(psc,eTot/E * pow(rate_gg+pow(tescape,-1),-1));
			st.ntPhoton.injection.set(psc,eTot/E * pow(rate_gg+pow(tescape,-1),-1));
		},{E_ix,-1,0});
	}
	
	for (size_t jE=0;jE<nEnt;jE++) {
		double E = st.ntPhoton.ps[DIM_E][jE];
		double fmtE = safeLog10(E/EV_TO_ERG);
		file << fmtE
			 << '\t' << safeLog10(pSy[jE]*E)
			 << "\t" << safeLog10(eIC[jE]*E)
			 << "\t" << safeLog10(pPP[jE]*E)
			 << "\t" << safeLog10(pPG[jE]*E)
			 << '\t' << safeLog10(totNotAbs[jE]*E)
			 << '\t' << safeLog10(totAbs[jE]*E)
			 << std::endl;
	}
	
	for (size_t jE=0;jE<nE;jE++) {
		double E = st.photon.ps[DIM_E][jE];
		double fmtE = safeLog10(E/EV_TO_ERG);
		fileSy << fmtE
			   << '\t' << safeLog10(eSy[jE]*E)
			   << std::endl;
	}
	
	file.close();
	fileSy.close();
	
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

