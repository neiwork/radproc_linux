#include "oneZoneTimeDependent.h"
#include "globalVariables.h"
#include <fparameters/parameters.h>
#include "adafFunctions.h"
#include <fparameters/SpaceIterator.h>
#include <fparameters/Dimension.h>
#include <fluminosities/luminositySynchrotron.h>
#include <fmath/fbisection.h>
#include <fmath/RungeKutta.h>
#include <fmath/physics.h>
#include <flosses/lossesSyn.h>
#include "injection.h"

using namespace std;

double cutOffPL2(double E, double Emin, double Emax)
{
	return pow(E,-pIndex)*exp(-E/Emax)*exp(-Emin/E);
}

void oneZoneDist(Particle& p, State& st) {
	
	p.ps.iterate([&](const SpaceIterator& iR) {
	
		double norm_temp = boltzmann*st.tempElectrons.get(iR)/(p.mass*cLight2);
		double dens = st.denf_e.get(iR);
        
		double aTheta = 3.0 - 6.0/(4.0+5.0*norm_temp); // Gammie & Popham (1998)
		double uth = dens*norm_temp*(p.mass*cLight2)*aTheta;   // erg cm^-3
		double Q0 = etaInj*uth;   // energy injected in the burst in nt particles [erg cm^-3]
		
		
		double rOriginal = iR.val(DIM_R);
		double rNow = rOriginal;
		double dr = 0.1*schwRadius;
		double tAccNow = 0.0;
		while (tAccNow < timeAfterFlare && rNow > exp(logr.front())) {
			tAccNow += dr/abs(radialVel(rNow));
			rNow -= dr;
		}
		double magf = st.magf.get(iR);
		
		double r = iR.val(DIM_R);
		double Emax = eEmax(p,r,magf,-radialVel(r),dens);
		double gammaMin = findGammaMin(st.tempElectrons.get(iR),Emax);
		double Emin = gammaMin*electronMass*cLight2;
		
		double int_E = RungeKuttaSimple(Emin,Emax,[&Emax,&Emin](double E){
			return E*cutOffPL2(E,Emin,Emax);});  //integra E*Q(E)  entre Emin y Emax
		
		double Q0p = Q0/int_E;
		
		p.ps.iterate([&](const SpaceIterator& iRE) {
			double E = iRE.val(DIM_E);
			
			if (rNow < maxRadius && rNow > exp(logr.front()))
				magf = magneticField(rNow);

			double dr = 0.1*schwRadius;
			double dEdt = 0.0;
			double tAccNow = 0.0;
			double rAux = rOriginal;
			while (tAccNow < timeAfterFlare && rAux > exp(logr.front())) {
				dEdt += (lossesSyn(E,magneticField(rAux),p)) * dr / abs(radialVel(rAux));
				tAccNow += dr/abs(radialVel(rNow));
				rAux -= dr;
			}
			dEdt /= timeAfterFlare;
			
			magf = magneticField(rOriginal);
			double tSync = E/lossesSyn(E,magf,p);
			//double tSync = E/dEdt;
			double tAcc = r/(-radialVel(r));
			double tEnd = min(tSync,tAcc);
			double Ne = (timeAfterFlare < tEnd && r < maxRadius && r > minRadius) ? 
							Q0p * cutOffPL2(E,Emin,Emax) * 
							pow(1.0-timeAfterFlare/tSync,pIndex-2.0) : 0.0;
			p.distribution.set(iRE,Ne);
		},{-1,iR.coord[DIM_R],0});
	},{0,-1,0});
}


void multiZoneInjection(Particle& p, State& st) {
	
	p.ps.iterate([&](const SpaceIterator& iR) {
	
		double norm_temp = boltzmann*st.tempElectrons.get(iR)/(p.mass*cLight2);
		double dens = st.denf_e.get(iR);
        
		double aTheta = 3.0 - 6.0/(4.0+5.0*norm_temp); // Gammie & Popham (1998)
		double uth = dens*norm_temp*(p.mass*cLight2)*aTheta;   // erg cm^-3
		double Q0 = etaInj*uth;   // energy injected in the burst in nt particles [erg cm^-3]
		
		double r = iR.val(DIM_R);
		double magf = st.magf.get(iR);
		double Emax = eEmax(p,r,magf,-radialVel(r),dens);
		double gammaMin = findGammaMin(st.tempElectrons.get(iR),Emax);
		double Emin = gammaMin*electronMass*cLight2;
		
		double int_E = RungeKuttaSimple(Emin,Emax,[&Emax,&Emin](double E){
			return E*cutOffPL2(E,Emin,Emax);});  //integra E*Q(E)  entre Emin y Emax
		
		double Q0p = Q0/int_E;
		
		p.ps.iterate([&](const SpaceIterator& iRE) {
			double E = iRE.val(DIM_E);
			double Qe = Q0p * cutOffPL2(E,Emin,Emax);
			p.injection.set(iRE,Qe);
		},{-1,iR.coord[DIM_R],0});
	},{0,-1,0});
}


void multiZoneDist(Particle& p, State& st, double dt) {
	
	size_t nPoints = 100;					// FOR THE CHARACTERISTIC LINE
	double dtChar = dt/(nPoints-1.0);
	p.ps.iterate([&](const SpaceIterator& iR) {
		
		const double r = iR.val(DIM_R);
		
		if (r < maxRadius) {

			Vector Rc(nPoints,r);
			// CHARACTERISTIC LINE ///////////////////////////
			for (size_t jRt=1;jRt<nPoints;jRt++) {
				double deriv = radialVel(Rc[jRt-1]);
				Rc[jRt] = Rc[jRt-1] + deriv*dtChar;
			}
			///////////////////////////////////////////////////

			p.ps.iterate([&](const SpaceIterator& iRE) {
				
				const double loge = log(iRE.val(DIM_E));
			
				if (Rc.back() > minRadius && Rc.back() < maxRadius) {

					Vector logEc(nPoints,loge);
					
					/// CHARACTERISTIC LINE /////////////////////
					for (size_t jEr=1;jEr<nPoints;jEr++) {
						double vr = 1.0;
						double be_e = 0.0;
						double deriv = 0.0;
						if (Rc[jEr-1] > exp(logr.front())) {
							vr = radialVel(Rc[jEr-1]);
							double magf = magneticField(Rc[jEr-1]);
							be_e = -1.575e-3*magf*magf*exp(logEc[jEr-1]);
							deriv = be_e/vr;
						}
						double dRc = Rc[jEr]-Rc[jEr-1];
						logEc[jEr] = logEc[jEr-1] + deriv*dRc;
					}
					///////////////////////////////////////////////////
				
					if (logEc.back() > log(p.emin()) && logEc.back() < log(p.emax())) {
						
						double e = exp(loge);
						double eLast = exp(logEc.back());
						double rLast = Rc.back();
						double vLast = radialVel(rLast);
						double v = radialVel(r);
						double densLast = electronDensity(rLast);
						double magfLast = magneticField(rLast);
						double EmaxLast = eEmax(p,rLast,magfLast,-vLast,densLast);
						double tempLast = electronTemp(rLast);
						double gammaMinLast = findGammaMin(tempLast,EmaxLast);
						double EminLast = gammaMinLast*p.mass*cLight2;
						
						double norm_tempLast = boltzmann*tempLast/(p.mass*cLight2);
						double aThetaLast = 3.0 - 6.0/(4.0+5.0*norm_tempLast); // Gammie & Popham (1998)
						double uthLast = densLast*norm_tempLast*(p.mass*cLight2)*aThetaLast;  // erg cm^-3
						double Q0Last = etaInj*uthLast;   // energy injected in the burst in nt particles [erg cm^-3]
						
						double int_E = RungeKuttaSimple(EminLast,EmaxLast,[&EmaxLast,&EminLast](double E){
									return E*cutOffPL2(E,EminLast,EmaxLast);});  //integra E*Q(E)  entre Emin y Emax
						double Q0p = Q0Last/int_E;
					
						double N = pow(rLast/r,2) * (vLast/v) * 
								pow(eLast/e,2.0) * Q0p*cutOffPL2(eLast,EminLast,EmaxLast);
						p.distribution.set(iRE,N);
					} else
						p.distribution.set(iRE,0.0);
				} else
					p.distribution.set(iRE,0.0);
			},{-1,iR.coord[DIM_R],0});
		} else {
			p.ps.iterate([&](const SpaceIterator& iRE) {
				p.distribution.set(iRE,0.0);
			},{-1,iR.coord[DIM_R],0});
		}
	},{0,-1,0});
}



void flareEmission(State& st, Particle& p, ofstream &file)
{
	double frequency = cLight/2.15e-4;
	double energy = planck*frequency;
	double lum = 0.0;
	p.ps.iterate([&](const SpaceIterator& iR) {
		double rOriginal = iR.val(DIM_R);
		double rNow = rOriginal;
		double dr = 0.1*schwRadius;
		double tAccNow = 0.0;
		while (tAccNow < timeAfterFlare && rNow > exp(logr.front())) {
			tAccNow += dr/abs(radialVel(rNow));
			rNow -= dr;
		}
		if (rNow < maxRadius && rNow > exp(logr.front())) {
			double magf = magneticField(rNow);
			double jv = luminositySynchrotron2(energy,p,iR,magf);
			double Inup = 2.0*height_fun(rNow)*jv;
			double vOrig = radialVel(rOriginal);
			double vNow = radialVel(rNow);
			lum += 2.0*pi*rOriginal*rOriginal*(sqrt(paso_r)-1.0/sqrt(paso_r))*Inup * 
							(vOrig/vNow) * (costhetaH(rOriginal)/costhetaH(rNow));
		}
	},{0,-1,0});
	lum += pow(10.0,34.2);
	file << timeAfterFlare << "\t" << log10(lum) << endl;
}


void flareEmission2(State& st, Particle& p, ofstream &file)
{
	double frequency = cLight/2.15e-4;
	double energy = planck*frequency;
	double lum = 0.0;
	p.ps.iterate([&](const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		if (r < maxRadius) {
			double magf = magneticField(r);
			double jv = luminositySynchrotron2(energy,p,iR,magf);
			double Inup = 2.0*height_fun(r)*jv;
			lum += 2.0*pi*r*r*(sqrt(paso_r)-1.0/sqrt(paso_r))*Inup;
		}
	},{0,-1,0});
	lum += pow(10.0,34.2);
	file << timeAfterFlare << "\t" << log10(lum) << endl;
}