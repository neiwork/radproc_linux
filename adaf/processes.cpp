#include "processes.h"

#include "modelParameters.h"
#include "write.h"
#include "messages.h"
#include "globalVariables.h"
#include "adafFunctions.h"
#include "ggAbsorption.h"

#include <fmath/RungeKutta.h>

#include <fluminosities/opticalDepthSSA.h>
#include <fluminosities/luminositySynchrotron.h>
#include <fluminosities/luminosityIC.h>
#include <fluminosities/luminosityNTHadronic.h>
#include <fluminosities/luminosityPhotoHadronic.h>

#include <fparameters/Dimension.h>
#include <fparameters/SpaceIterator.h>

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

	SpaceCoord psc = {E_ix,0,0};
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
		do {
			double xprim = rprim*sin(thetaprim);
			double zprim = rprim*cos(thetaprim);
			double x1 = x0+xprim;
			double z1 = z0+zprim;
			r1=sqrt(x1*x1+y0*y0+z1*z1);
			
			double magneticField = (r1 < rMax && r1 > rMin && 
							theta1 > thetaMinLocal && theta1 < pi-thetaMinLocal) ? 
							st.magf.interpolate({{1,r1}},&psc) : 0.0;
			double cosThetaMinLocal = costhetaH(r1);
			theta1 = atan(sqrt(x1*x1+y0*y0)/abs(z1));
			thetaMinLocal = acos(cosThetaMinLocal);
			double integral = intSimple(p.emin(),p.emax(),[&](double x){
							return fSSA2(x,E,p,magneticField,psc);});
			
			double absorptionCoefficient = -P3(planck)*cLight2*integral/(8*pi*P2(E)); 
			opticalDepth += absorptionCoefficient*drprim;
			
			drprim = r1*(pasoprim-1.0);
			rprim += drprim;
		} while (r1 < rMax && r1 > rMin && theta1 > thetaMinLocal && theta1 < pi-thetaMinLocal);
		phi0 += dPhi;
	}
	return opticalDepth/nPhi;
}

/* Takes [emi] =  E^2*[Q(E)] and calculates int(2.0*pi*P2(jetR)*emi dz); 
for [N(E)] = 1/erg, then it just sums over all z and returns erg/s  */

void processes(State& st, const std::string& filename)
{
	show_message(msgStart, Module_luminosities);

	std::ofstream file;
	file.open(filename.c_str(),std::ios::out);
	file << "log(E/eV)"
			<< '\t' << "eSyn"
			<< '\t' << "pSyn"
			<< "\t" << "eIC"
			<< "\t" << "pPP"
			<< "\t" << "pPG"
			<< std::endl;
	
	double Emin = st.photon.emin();
	double Emax = st.photon.emax();
	int nE = st.photon.ps[DIM_E].size();
	
	for (int E_ix=0;E_ix<nE;E_ix++) {

		double E = st.photon.distribution.ps[DIM_E][E_ix];
		double fmtE = log10(E / 1.6e-12);
		double eSynNotAbs,eICNotAbs,pSynNotAbs,pPPNotAbs,pPGNotAbs;
		double eSynAbs,eICAbs,pSynAbs,pPPAbs,pPGAbs;
		eSynNotAbs = eICNotAbs = pSynNotAbs = pPPNotAbs = pPGNotAbs = 0.0;
		eSynAbs = eICAbs = pSynAbs = pPPAbs = pPGAbs = 0.0;

		st.photon.ps.iterate([&](const SpaceIterator &i) {
			double r = i.val(DIM_R);
			double rB1 = r/sqrt(paso_r);
			double rB2 = r*sqrt(paso_r);
			double vol = (4.0/3.0)*pi*cos(st.thetaH.get(i))*(rB2*rB2*rB2-rB1*rB1*rB1);
			
			//double tau_e = opticalDepthSSA(E_ix, st, st.ntElectron, i.val(DIM_R));
			double factorSSA_e = 1.0;
			double factorSSA_p = 1.0;
			double factorGG = 1.0;
			
			if (fmtE < 5.0) { // para que no calcule a todas las energías
				double tau_e = opticalDepthSSA2(E_ix,st,st.ntElectron,i.coord[DIM_R]);
				factorSSA_e = (tau_e > 1.0e-15) ? (1.0-exp(-tau_e))/tau_e : 1.0;
				double tau_p = opticalDepthSSA(E_ix, st, st.ntProton, i.val(DIM_R));
				factorSSA_p = (tau_p > 1.0e-15) ? (1.0-exp(-tau_p))/tau_p : 1.0;
			}
			if (fmtE > 5.0) factorGG = exp(-ggOpticalDepth(E_ix,st,i.coord[DIM_R]));
			//double tau_gg = internalAbs(E_ix,st,r);
			//cout << "E = " << E << "\t r = " << r << "\t factor = " << factorSSA_e << endl;
			double eSyn = luminositySynchrotron(E,st.ntElectron,i,st.magf)*vol;
			double eIC  = luminosityIC(E,st.ntElectron,i.coord,st.photon.distribution,Emin)*vol;//ver unidades del distribution XXX
			double pSyn = luminositySynchrotron(E,st.ntProton,i,st.magf)*vol;
			double pPP  = luminosityNTHadronic(E,st.ntProton,st.denf_i.get(i),i)*vol;
			double pPG  = luminosityPhotoHadronic(E,st.ntProton,st.photon.distribution,i,Emin,Emax)*vol;
			eSynNotAbs += eSyn;
			eICNotAbs += eIC;
			pSynNotAbs += pSyn;
			pPPNotAbs += pPP;
			pPGNotAbs += pPG;
			eSynAbs += eSyn*factorSSA_e*factorGG;
			eICAbs += eIC*factorGG;
			pSynAbs += pSyn*factorSSA_p*factorGG;
			pPPAbs += pPP*factorGG;
			pPGAbs += pPG*factorGG;
		},{E_ix,-1,0});
		
		file << fmtE
			 << '\t' << safeLog10(eSynNotAbs)
			 << '\t' << safeLog10(pSynNotAbs)
			 << "\t" << safeLog10(eICNotAbs)
			 << "\t" << safeLog10(pPPNotAbs)
			 << "\t" << safeLog10(pPGNotAbs)
			 << '\t' << safeLog10(eSynAbs)
			 << '\t' << safeLog10(pSynAbs)
			 << "\t" << safeLog10(eICAbs)
			 << "\t" << safeLog10(pPPAbs)
			 << "\t" << safeLog10(pPGAbs)
			 << std::endl;
	}
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

