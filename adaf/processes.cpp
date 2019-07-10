#include "processes.h"



#include "modelParameters.h"
#include "write.h"
#include "messages.h"
#include "globalVariables.h"

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
		
			double thetaH = st.thetaH.get(i);
			double rB1=r/sqrt(paso_r);
			double rB2=r*sqrt(paso_r);
			double delta_r = rB2-rB1; //(r/sqrt(paso_r))*(paso_r-1.0);
			double vol=(rB2*rB2*rB2-rB1*rB1*rB1)*(4.0/3.0)*pi*cos(thetaH);
		
			double magneticField = st.magf.get(i);
				
			double integral = intSimple(p.emin(), p.emax(), [&](double x) {
								return fSSA(x, E, p, magneticField, i);
								});

			double absorptionCoefficient = - P3(planck)*cLight2*integral/(8*pi*P2(E)); 
		
			opacity += absorptionCoefficient*delta_r/vol; // parameters.radius;
		}
		
	},{E_ix,-1,0});
		
		return opacity;
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
		double eSyn,eIC,pSyn,pPP,pPG;
		eSyn = eIC = pSyn = pPP = pPG = 0.0;

		st.photon.ps.iterate([&](const SpaceIterator &i) {
			
			double tau_e = opticalDepthSSA(E_ix, st, st.ntElectron, i.val(DIM_R));   //E=Eph; i.val(DIM_R) = r_current
			double factorSSA_e = 1.0;
			
			if (tau_e > 1.0e-15){factorSSA_e = (1.0-exp(-tau_e))/tau_e;	} 
			//esto esta para que ene el limite de tau = 0 el factor de 1 y no 0
			
			double tau_p = opticalDepthSSA(E_ix, st, st.ntProton, i.val(DIM_R));   //E=Eph; i.val(DIM_R) = r_current
			double factorSSA_p = 1.0;
			
			if (tau_p > 1.0e-15){factorSSA_p = (1.0-exp(-tau_p))/tau_p;	}
	
			eSyn += factorSSA_e*luminositySynchrotron(E,st.ntElectron,i, st.magf);
			eIC  += 0.0;//luminosityIC(E,st.ntElectron,i.coord,st.photon.distribution,Emin);//ver unidades del distribution XXX
			pSyn += factorSSA_p*luminositySynchrotron(E,st.ntProton,i,st.magf);
			pPP  += 0.0;//luminosityNTHadronic(E,st.ntProton,st.denf_i.get(i),i);
			pPG  += 0.0;//luminosityPhotoHadronic(E,st.ntProton,st.photon.distribution,i,Emin,Emax);
		},{E_ix,-1,0});
		
		file << fmtE
			 << '\t' << safeLog10(eSyn)
			 << '\t' << safeLog10(pSyn)
			 << "\t" << safeLog10(eIC)
			 << "\t" << safeLog10(pPP)
			 << "\t" << safeLog10(pPG)
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

