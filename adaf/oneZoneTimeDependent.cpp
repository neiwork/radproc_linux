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
	static const double primaryIndex = GlobalConfig.get<double>("nonThermal.flare.injection.primaryIndex");
	return pow(E,-primaryIndex)*exp(-E/Emax)*exp(-Emin/E);
}

void oneZoneDist(Particle& p, State& st) {
	
	//double timeAfterFlare = GlobalConfig.get<double>("nonThermal.flare.timeAfterFlare");
	double maxRadius = GlobalConfig.get<double>("nonThermal.flare.maxRadius")*schwRadius;
	double minRadius = GlobalConfig.get<double>("nonThermal.flare.minRadius")*schwRadius;
	static const double etaInj = GlobalConfig.get<double>("nonThermal.flare.injection.energyFraction");
	static const double pIndex = GlobalConfig.get<double>("nonThermal.flare.injection.primaryIndex");
	double Emin = p.emin();
	double Emax = p.emax();
	
	p.ps.iterate([&](const SpaceIterator& iR) {
	
		double norm_temp = boltzmann*st.tempElectrons.get(iR)/(p.mass*cLight2);
		double dens = st.denf_e.get(iR);
        
		double aTheta = 3.0 - 6.0/(4.0+5.0*norm_temp); // Gammie & Popham (1998)
		double uth = dens*norm_temp*(p.mass*cLight2)*aTheta;   // erg cm^-3
		double Q0 = etaInj*uth;   // energy injected in the burst in nt particles [erg cm^-3]
        double magf = st.magf.get(iR);
		double r = iR.val(DIM_R);
		
		Emax = eEmax(p,r,magf,-radialVel(r),dens);
		double gammaMin = findGammaMin(st.tempElectrons.get(iR),Emax);
		Emin = gammaMin*electronMass*cLight2;
		
		double int_E = RungeKuttaSimple(Emin,Emax,[&Emax,&Emin](double E){
			return E*cutOffPL2(E,Emin,Emax);});  //integra E*Q(E)  entre Emin y Emax
		
		double Q0p = Q0/int_E;
		
		p.ps.iterate([&](const SpaceIterator& iRE) {
			double E = iRE.val(DIM_E);
			double tSync = E/lossesSyn(E,magf,p);
			double tAcc = r/(-radialVel(r));
			double tEnd = min(tSync,tAcc);
			double Ne = (timeAfterFlare < tEnd && r < maxRadius && r > minRadius) ? 
							Q0p * cutOffPL2(E,Emin,Emax) * 
							pow(1.0-timeAfterFlare/tSync,pIndex-2.0) : 0.0;
			p.distribution.set(iRE,Ne);
		},{-1,iR.coord[DIM_R],0});
	},{0,-1,0});
}

void flareEmission(State& st, Particle& p, ofstream &file)
{
	double maxRadius = GlobalConfig.get<double>("nonThermal.flare.maxRadius")*schwRadius;
	double minRadius = GlobalConfig.get<double>("nonThermal.flare.minRadius")*schwRadius;
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
			double magf = sqrt((1.0-magFieldPar)*8*pi*massDensityADAF(rNow)*sqrdSoundVel(rNow));
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