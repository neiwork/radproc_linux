#include "injection.h"

#include "messages.h"
#include "globalVariables.h"
#include "adafFunctions.h"
#include <fparameters/parameters.h>
#include <fparameters/SpaceIterator.h>
#include <fparameters/Dimension.h>
#include <fmath/fbisection.h>
#include <fmath/RungeKutta.h>
#include <fmath/physics.h>
#include <boost/property_tree/ptree.hpp>
#include <iostream>

double eEmax(Particle& p, double r, double B, double v, double dens)
{
	double accE = GlobalConfig.get<double>("accEfficiency");
    double Emax_adv = accE*r*cLight*electronCharge*B / v; 
    double Emax_syn = p.mass*cLight2*sqrt(accE*6.0*pi*electronCharge / (thomson*B));
	if (p.id == "ntProton") {
		double sigmapp = 34.3e-27;
		double Emax_pp = accE*electronCharge*B/(0.5*dens*sigmapp);
		return min(Emax_adv,Emax_pp);
	} else
		return min(Emax_syn,Emax_adv);
}

double powerLaw(double E, double Emin, double Emax)
{
	static const double primaryIndex = GlobalConfig.get<double>("primaryIndex");
	double result = pow(E,-primaryIndex)*exp(-E/Emax)*exp(-5*Emin/E);
	return result;
}

void injection(Particle& p, State& st)
{
	static const double etaInj = GlobalConfig.get<double>("etaInj");
	double Emin = p.emin();   //esta es la primera que uso de prueba
	p.ps.iterate([&](const SpaceIterator& i) {
		const double r = i.val(DIM_R);
		const double thetaH = st.thetaH.get(i);
		const double area  = r*r*(2.0*cos(thetaH)+sin(thetaH)*sin(thetaH));

		double Emax = eEmax(p,r,st.magf.get(i),-radialVel(r),st.denf_i.get(i));
		double int_E = RungeKuttaSimple(Emin,Emax,[&Emax,&Emin](double E){
			return E*powerLaw(E, Emin, Emax);});  //integra E*Q(E)  entre Emin y Emax

		double norm_temp, dens;
		if(p.id == "ntElectron"){
			norm_temp = boltzmann*st.tempElectrons.get(i)/(p.mass*cLight2);
			dens = st.denf_e.get(i);
		} else if(p.id == "ntProton"){
			norm_temp = boltzmann*st.tempIons.get(i)/(p.mass*cLight2);
			dens = st.denf_i.get(i);
		}
		double aTheta = 3.0 - 6.0/(4.0+5.0*norm_temp); // Gammie & Popham (1998)
		double uth = dens*norm_temp*(p.mass*cLight2)*aTheta;   // erg cm^-3
		double Q0 = etaInj*uth*(area*cLight);   // power injected in nt particles [erg s^{-1}]
		double Q0p = Q0/int_E;
		p.ps.iterate([&](const SpaceIterator& jE) {
			const double E = jE.val(DIM_E);
			double total = powerLaw(E, Emin, Emax)*Q0p; 
			p.injection.set(jE,total); //en unidades de erg^-1 s^-1
		},{-1,i.coord[DIM_R],0});
	},{0,-1,0});
}
