#include "injection.h"

#include "messages.h"
#include "globalVariables.h"
#include "adafFunctions.h"

#include <finjection/pairInjectionExact.h>
#include <finjection/pairInjection.h>

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
	double accE = GlobalConfig.get<double>("nonThermal.injection.accEfficiency");
    double Emax_adv = accE*r*cLight*electronCharge*B / v; 
    double Emax_syn = p.mass*cLight2*sqrt(accE*6.0*pi*electronCharge / (thomson*B));
    double Emax_Hillas = electronCharge*B*r;
	if (p.id == "ntProton") {
		double sigmapp = 34.3e-27;
		double Emax_pp = accE*electronCharge*B/(0.5*dens*sigmapp);
		return min(min(Emax_adv,Emax_pp),Emax_Hillas);
	} else
		return min(min(Emax_syn,Emax_adv),Emax_Hillas);
}

double cutOffPL(double E, double Emin, double Emax)
{
	static const double primaryIndex = GlobalConfig.get<double>("nonThermal.injection.primaryIndex");
	return pow(E,-primaryIndex)*exp(-E/Emax)*exp(-5*Emin/E);
}

void injection(Particle& p, State& st)
{
	static const double etaInj = GlobalConfig.get<double>("nonThermal.injection.energyFraction");
	double Emin = p.emin();   //esta es la primera que uso de prueba
	
    double sumQ = 0.0;
	p.ps.iterate([&](const SpaceIterator& i) {
		const double r = i.val(DIM_R);
		double rB1 = r/sqrt(paso_r);
		double rB2 = rB1*paso_r;
		const double thetaH = st.thetaH.get(i);
		const double vol = (4.0/3.0)*pi*cos(thetaH)*(rB2*rB2*rB2-rB1*rB1*rB1);

		double Emax = eEmax(p,r,st.magf.get(i),-radialVel(r),st.denf_i.get(i));
		double int_E = RungeKuttaSimple(Emin,p.emax(),[&Emax,&Emin](double E){
			return E*cutOffPL(E,Emin,Emax);});  //integra E*Q(E)  entre Emin y Emax
        
        /*
		double grpersecToSolarMassesperyear = 3600.0*24*265.25/solarMass;
		double Wmr = 1.0e42*4.0*accRateADAF(r)*grpersecToSolarMassesperyear * 
					(st.tempIons.get(i)/st.tempElectrons.get(i)); // [erg s^{-1}]
                    
        sumWmr += accRateADAF(r);
        
		double Q0 = 0.0;
		double delta = 0.1;
		if (p.id == "ntElectron") Q0 = delta * Wmr;
		else if (p.id == "ntProton") Q0 = (1.0-delta)* Wmr;*/
		
        
        double norm_temp, dens;
		if(p.id == "ntElectron") {
			norm_temp = boltzmann*st.tempElectrons.get(i)/(p.mass*cLight2);
			dens = st.denf_e.get(i);
		} else if(p.id == "ntProton") {
			norm_temp = boltzmann*st.tempIons.get(i)/(p.mass*cLight2);
			dens = st.denf_i.get(i);
		}
        
		double aTheta = 3.0 - 6.0/(4.0+5.0*norm_temp); // Gammie & Popham (1998)
		double uth = dens*norm_temp*(p.mass*cLight2)*aTheta;   // erg cm^-3
		double Q0 = etaInj*uth*cLight/r;   // power injected in nt particles [erg s^-1 cm^-3]
        
		double Q0p = Q0/int_E;
		p.ps.iterate([&](const SpaceIterator& jE) {
			const double E = jE.val(DIM_E);
			double total = cutOffPL(E, Emin, Emax)*Q0p;
			p.injection.set(jE,total); //en unidades de erg^-1 s^-1 cm^-3
		},{-1,i.coord[DIM_R],0});
		sumQ += Q0*vol;
	},{0,-1,0});
	cout << "Total power injected in " << p.id << " = " << sumQ << endl; 
}

void injectionPair(Particle& p, State& st)
{
	//static const double etaInj = GlobalConfig.get<double>("nonThermal.injection.energyFraction");
	//double Emin = p.emin();   //esta es la primera que uso de prueba
	
    //double sumQ = 0.0;
	show_message(msgStart, Module_pairInjection);
	p.ps.iterate([&](const SpaceIterator& i) {
		const double E = i.val(DIM_E);

		// en donde pongo st.photon.injection debería ir el paramSpaceValue con la densidad de fotones no termicos
		double result = pairInjectionExact(E, st.photon.injection, st.photon.distribution, i, st.photon.emin(), st.photon.emax());
		double result2 = pairInjection(E, st.photon.injection, st.photon.distribution, i, st.photon.emin(), st.photon.emax());
		
		p.injection.set(i,result);
		p.distribution.set(i,result2);  //en rclTabCtrlealidad son inyecciones ambas, lo hago asi para comparar las dos aproximaciones y ver cual usamos
	},{-1,-1,0});
	show_message(msgEnd, Module_pairInjection);
	//cout << "Total power injected in " << p.id << " = " << sumQ << endl; 
}
