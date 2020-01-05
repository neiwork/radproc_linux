#include "flareTimeDependent.h"
#include "globalVariables.h"
#include <fparameters/parameters.h>
#include "adafFunctions.h"
#include <fparameters/SpaceIterator.h>
#include <fparameters/Dimension.h>
#include <fluminosities/luminositySynchrotron.h>
#include <fluminosities/luminosityIC.h>
#include <fluminosities/thermalSync.h>
#include <fluminosities/blackBody.h>
#include <fmath/fbisection.h>
#include <fmath/RungeKutta.h>
#include <fmath/physics.h>
#include <flosses/lossesSyn.h>
#include "injection.h"
#include "absorption.h"
#include "write.h"

using namespace std;


double eEmaxFlare(Particle& p, double r, double B, double v, double dens)
{
	double accE = GlobalConfig.get<double>("nonThermal.flare.injection.accEfficiency");
	double dr = r*(paso_r-1.0);
    double Emax_adv = accE*dr*cLight*electronCharge*B / v; 
    double Emax_syn = p.mass*cLight2*sqrt(accE*6.0*pi*electronCharge / (thomson*B));
    double Emax_Hillas = electronCharge*B*dr;
	if (p.id == "ntProton") {
		double sigmapp = 34.3e-27;
		double Emax_pp = accE*electronCharge*B/(0.5*dens*sigmapp);
		return min(min(Emax_adv,Emax_pp),Emax_Hillas);
	} else
		return min(min(Emax_syn,Emax_adv),Emax_Hillas);
}

double cutOffPL2(double E, double Emin, double Emax, Particle& p)
{
	Emax = min(Emax,p.emax());
	Emin = max(Emin,p.emin());
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
		while (tAccNow < timeAfterFlare && rNow > exp(logr.front())*schwRadius) {
			tAccNow += dr/abs(radialVel(rNow));
			rNow -= dr;
		}
		double magf = st.magf.get(iR);
		
		double r = iR.val(DIM_R);
		double Emax = eEmax(p,r,magf,-radialVel(r),dens);
		double gammaMin = findGammaMin(st.tempElectrons.get(iR),Emax);
		double Emin = gammaMin*electronMass*cLight2;
		
		double int_E = RungeKuttaSimple(Emin,Emax,[&Emax,&Emin,&p](double E){
			return E*cutOffPL2(E,Emin,Emax,p);});  //integra E*Q(E)  entre Emin y Emax
		
		double Q0p = Q0/int_E;
		
		p.ps.iterate([&](const SpaceIterator& iRE) {
			double E = iRE.val(DIM_E);
			
			if (rNow < maxRadius && rNow > exp(logr.front())*schwRadius)
				magf = magneticField(rNow);

			double dr = 0.1*schwRadius;
			double dEdt = 0.0;
			double tAccNow = 0.0;
			double rAux = rOriginal;
			while (tAccNow < timeAfterFlare && rAux > exp(logr.front())*schwRadius) {
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
							Q0p * cutOffPL2(E,Emin,Emax,p) * 
							pow(1.0-timeAfterFlare/tSync,pIndex-2.0) : 0.0;
			p.distribution.set(iRE,Ne);
		},{-1,iR.coord[DIM_R],0});
	},{0,-1,0});
}


void multiZoneInjection(Particle& p, State& st) {
	
	double sum = 0.0;
	double sumMagF = 0.0;
	p.ps.iterate([&](const SpaceIterator& iR) {
	
		double r = iR.val(DIM_R);
		if (r > minRadius && r < maxRadius) {
			double norm_temp = boltzmann*st.tempElectrons.get(iR)/(p.mass*cLight2);
			double dens = st.denf_e.get(iR);
			double rB2 = r*sqrt(paso_r);
			double rB1 = rB2/paso_r;
			double vol = (4.0/3.0)*pi*(rB2*rB2*rB2-rB1*rB1*rB1)*costhetaH(r);
			
			double aTheta = 3.0 - 6.0/(4.0+5.0*norm_temp); // Gammie & Popham (1998)
			double uth = dens*norm_temp*(p.mass*cLight2)*aTheta;   // erg cm^-3
			double Q0 = etaInj*uth;   // energy injected in the burst in nt particles [erg cm^-3]
			double magf = st.magf.get(iR);
			double Emax = eEmax(p,r,magf,-radialVel(r),dens);
			double gammaMin = findGammaMin(st.tempElectrons.get(iR),Emax);
			double Emin = gammaMin*electronMass*cLight2;
			
			double int_E = RungeKuttaSimple(p.emin(),p.emax(),[&Emax,&Emin,&p](double E){
				return E*cutOffPL2(E,Emin,Emax,p);});  //integra E*Q(E)  entre Emin y Emax
			
			double Q0p = Q0/int_E;
			
			
			//Q0 = Q0 * ( r < 10.0*schwRadius ? 1.0/(10.0*schwRadius)*r : 1.0 );
			p.ps.iterate([&](const SpaceIterator& iRE) {
				double E = iRE.val(DIM_E);
				double Qe = Q0p * cutOffPL2(E,Emin,Emax,p);
				p.injection.set(iRE,Qe);
			},{-1,iR.coord[DIM_R],0});
			sum += Q0*vol;
			sumMagF += magf*magf/(8.0*pi)*vol;
		} else {
			p.ps.iterate([&](const SpaceIterator& iRE) {
				p.injection.set(iRE,0.0);
			},{-1,iR.coord[DIM_R],0});
		}
	},{0,-1,0});
	cout << "Total Bomb Energy = " << sum << endl;
	cout << "Total Magnetic Energy = " << sumMagF << endl;
}


void multiZoneDist(Particle& p, State& st, double dt) {
	
	size_t nPoints = 100;					// FOR THE CHARACTERISTIC LINE
	double dtChar = dt/(nPoints-1.0);
	
	p.ps.iterate([&](const SpaceIterator& iR) {
		
		const double r = iR.val(DIM_R);
		
		if (r < maxRadius && r > schwRadius*2) {
			
			Vector logRc(nPoints,log(r));
			// CHARACTERISTIC LINE ///////////////////////////
			for (size_t jRt=1;jRt<nPoints;jRt++) {
				double rLocal = exp(logRc[jRt-1]);
				double deriv = radialVel(rLocal)/rLocal;
				logRc[jRt] = logRc[jRt-1] + deriv*dtChar;
			}
			///////////////////////////////////////////////////

			p.ps.iterate([&](const SpaceIterator& iRE) {
				
				const double loge = log(iRE.val(DIM_E));
			
				if (logRc.back() > log(minRadius) && logRc.back() < log(maxRadius)) {

					Vector logEc(nPoints,loge);
					
					/// CHARACTERISTIC LINE /////////////////////
					for (size_t jEr=1;jEr<nPoints;jEr++) {
						double rLocal = exp(logRc[jEr-1]);
						double vr_r = 1.0;
						double be_e = 0.0;
						double deriv = 0.0;
						if (rLocal > exp(logr.front())*schwRadius) {
							vr_r = radialVel(rLocal)/rLocal;
							double magf = magneticField(rLocal);
							be_e = -1.575e-3*magf*magf*exp(logEc[jEr-1]);
							deriv = be_e/vr_r;
						}
						double dlogRc = logRc[jEr]-logRc[jEr-1];
						logEc[jEr] = logEc[jEr-1] + deriv*dlogRc;
					}
					///////////////////////////////////////////////////
				
					if (logEc.back() > log(p.emin()) && logEc.back() < log(p.emax())) {
						
						double e = exp(loge);
						double eLast = exp(logEc.back());
						double rLast = exp(logRc.back());
						double vLast = radialVel(rLast);
						double v = radialVel(r);
						double densLast = electronDensity(rLast);
						double magfLast = magneticField(rLast);
						double EmaxLast = eEmaxFlare(p,rLast,magfLast,-vLast,densLast);
						double tempLast = electronTemp(rLast);
						double gammaMinLast = findGammaMin(tempLast,EmaxLast);
						double EminLast = gammaMinLast*p.mass*cLight2;
						
						double norm_tempLast = boltzmann*tempLast/(p.mass*cLight2);
						double aThetaLast = 3.0 - 6.0/(4.0+5.0*norm_tempLast); // Gammie & Popham (1998)
						double uthLast = densLast*norm_tempLast*(p.mass*cLight2)*aThetaLast;  // erg cm^-3
						double Q0Last = etaInj*uthLast;   // energy injected in the burst in nt particles [erg cm^-3]
						
						//Q0Last = Q0Last * ( (r < 10.0*schwRadius) ? 1.0/(10.0*schwRadius)*r : 1.0 );
						
						double int_E = RungeKuttaSimple(EminLast,EmaxLast,[&EmaxLast,&EminLast,&p](double E){
									return E*cutOffPL2(E,EminLast,EmaxLast,p);});  //integra E*Q(E)  entre Emin y Emax
						double Q0p = Q0Last/int_E;
					
						double N = pow(rLast/r,2) * (vLast/v) *
								pow(eLast/e,2.0) * Q0p*cutOffPL2(eLast,EminLast,EmaxLast,p); //(costhetaH(rLast)/costhetaH(r))
						
						// SOLO ENFRIAMIENTO (SIN ACRECIÓN)
						/*
						double Emax = eEmax(p,r,st.magf.get(iRE),-radialVel(r),st.denf_e.get(iRE));
						double gammaMin = findGammaMin(st.tempElectrons.get(iRE),Emax);
						double Emin = gammaMin*electronRestEnergy;
						double int_Eaux = RungeKuttaSimple(Emin,Emax,[&Emax,&Emin,&p](double e){
									return e*cutOffPL2(e,Emin,Emax,p);});
						double magfield = magneticField(r);
						double tEnd = e/lossesSyn(e,magfield,p);
						double normTemp = boltzmann*st.tempElectrons.get(iRE)/electronRestEnergy;
						double aTheta = 3.0 - 6.0/(4.0+5.0*normTemp);
						double uth = st.denf_e.get(iRE)*normTemp*electronRestEnergy*aTheta;
						double Q0aux = etaInj*uth / int_Eaux;
						N = (timeAfterFlare < tEnd && r < maxRadius && r > minRadius) ? 
							Q0aux * cutOffPL2(e,Emin,Emax,p) * 
							pow(1.0-timeAfterFlare/tEnd,pIndex-2.0) : 0.0; */
						
						/*double logeLast = loge;
						double tAux = 0.0;
						while (tAux < timeAfterFlare) {
							eLast = exp(logeLast);
							logeLast += 1.575e-3*magfield*magfield*eLast * (-dtChar);
							tAux -= dtChar;
						}
						eLast = exp((log(eLast)+logeLast)/2.0);
						if (eLast < p.emax() && eLast > p.emin() && tEnd > timeAfterFlare && 
								r < maxRadius && r > minRadius)
							N = pow(eLast/e,2.0) * Q0aux*cutOffPL2(eLast,Emin,Emax,p);
						else
							N = 0.0;*/
						///////////////////////////////////
						
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
		while (tAccNow < timeAfterFlare && rNow > exp(logr.front())*schwRadius) {
			tAccNow += dr/abs(radialVel(rNow));
			rNow -= dr;
		}
		if (rNow < maxRadius && rNow > exp(logr.front())*schwRadius) {
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


void flareEmission2(State& st, Particle& p, ofstream &file_cm, ofstream &file_mm, 
						ofstream &fileNIR, ofstream &fileX)
{
	double frequencycm = cLight/1.3;
	double frequencymm = cLight/1.3e-1;
	double frequencyNIR = cLight/2.15e-4;
	double frequencyX = 1.2e18;
	double energycm = planck*frequencycm;
	double energymm = planck*frequencymm;
	double energyNIR = planck*frequencyNIR;
	double energyX = planck*frequencyX;
	double lum_cm = 0.0;
	double lum_mm = 0.0;
	double lumNIR = 0.0;
	double lumX = 0.0;
	double sumTot = 0.0;
	Vector SED(nE,0.0);
	p.ps.iterate([&](const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		if (r < maxRadius) {
			
			double rB2 = r*sqrt(paso_r);
			double rB1 = rB2/paso_r;
			double vol = (4.0/3.0)*pi*(rB2*rB2*rB2-rB1*rB1*rB1)*costhetaH(r);
			double area = 4.0*pi*r*r*costhetaH(r);
			double temp = st.tempElectrons.get(iR);
			double dens_e = st.denf_e.get(iR);
			double magf = magneticField(r);
			
			/*st.photon.ps.iterate([&](const SpaceIterator& iRE) {
				double E = iRE.val(DIM_E);
				double nPh = vol * luminositySynchrotron2(E,p,iRE,magf) / (area*cLight*E*E);
				st.photon.distribution.set(iRE,nPh);
			},{-1,iR.coord[DIM_R],0});*/

			double jvcm = luminositySynchrotron2(energycm,p,iR,magf)/(4*pi);
			double jvmm = luminositySynchrotron2(energymm,p,iR,magf)/(4*pi);
			double jvNIR = luminositySynchrotron2(energyNIR,p,iR,magf)/(4*pi);
			double jvX = 0.0;
			
			///////////
			//double jvX = (luminositySynchrotron2(energyX,p,iR,magf) + 
			//		luminosityIC(energyX,p,iR,st.photon.distribution,st.photon.emin()))/(4.0*pi);
			
			double jvTh_cm = jSync(energycm,temp,magf,dens_e);
			double jvTh_mm = jSync(energymm,temp,magf,dens_e);
			
			double avTh_cm = jvTh_cm / bb(frequencycm,temp);
			double avTh_mm = jvTh_mm / bb(frequencymm,temp);
			SpaceCoord psc = {0,iR.coord[DIM_R],0};
			double avPl_cm = ssaAbsorptionCoeff(energycm,magf,p,psc);
			double avPl_mm = ssaAbsorptionCoeff(energymm,magf,p,psc);
			
			double height = r*costhetaH(r);
			double Sv_cm = jvcm/(avTh_cm+avPl_cm);
			double Sv_mm = jvmm/(avTh_mm+avPl_mm);
			double Inup_cm = Sv_cm * (1.0-exp(-2.0*height*(avTh_cm+avPl_cm)));
			double Inup_mm = Sv_mm * (1.0-exp(-2.0*height*(avTh_mm+avPl_mm)));
			double flux_cm = 2.0*pi*r*r*(paso_r-1.0)*Inup_cm;
			double flux_mm = 2.0*pi*r*r*(paso_r-1.0)*Inup_mm;
			
			double lumLocalcm = 4.0*pi*flux_cm;
			double lumLocalmm = 4.0*pi*flux_mm;
			double lumLocalNIR = 4.0*pi*jvNIR*vol;
			double lumLocalX = 4.0*pi*jvX*vol;
			//double Inup = 2.0*height_fun(r)*jv;
			//lum += 2.0*pi*r*r*(paso_r-1.0)*Inup;
			lum_cm += lumLocalcm;
			lum_mm += lumLocalmm;
			lumNIR += lumLocalNIR;
			lumX += lumLocalX;
		}
	},{0,-1,0});
	
	lum_cm += pow(10.0,33.23);
	lum_mm += pow(10.0,34.85);
	lumNIR += pow(10.0,34.26);
	lumX += pow(10.0,31.7);
	file_cm << timeAfterFlare << "\t" << log10(lum_cm) << endl;
	file_mm << timeAfterFlare << "\t" << log10(lum_mm) << endl;
	fileNIR << timeAfterFlare << "\t" << log10(lumNIR) << endl;
	fileX << timeAfterFlare << "\t" << log10(lumX) << endl;
}