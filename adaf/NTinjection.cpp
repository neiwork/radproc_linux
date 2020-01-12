#include "NTinjection.h"
#include "write.h"
#include "messages.h"
#include "globalVariables.h"
#include "adafFunctions.h"
#include "absorption.h"
#include <finjection/ppPionInj.h>
#include <finjection/pgammaPionInj.h>
#include <finjection/muonInj.h>
#include <flosses/lossesHadronics.h>
#include <flosses/lossesSyn.h>
#include <finjection/pairInjectionExact.h>
#include <finjection/pairgg.h>
#include <finjection/pairMuonDecay.h>
#include <finjection/pairBH.h>
#include <gsl/gsl_sf_bessel.h>

#include <fparameters/parameters.h>
#include <fparameters/SpaceIterator.h>
#include <fparameters/Dimension.h>
#include <fmath/fbisection.h>
#include <fmath/RungeKutta.h>
#include <fmath/physics.h>
#include <boost/property_tree/ptree.hpp>
#include <iostream>
#include <fstream>

double eEmax(Particle& p, double r, double B, double v, double dens)
{
	double accE = GlobalConfig.get<double>("nonThermal.injection.PL.accEfficiency");
	double dr = r*(paso_r-1.0);
	double h = height_fun(r);
    double Emax_adv = accE*dr*cLight*electronCharge*B / v; 
    double Emax_syn = p.mass*cLight2*sqrt(accE*6.0*pi*electronCharge / (thomson*B)) * p.mass/electronMass;
    double Emax_Hillas = electronCharge*B*h;
	double zeda = GlobalConfig.get<double>("nonThermal.injection.SDA.fractionTurbulent");
	double q = GlobalConfig.get<double>("nonThermal.injection.SDA.powerSpectrumIndex");
	double Emax_diff = pow(h*electronCharge*B,(q-3)/(2*q-5)) * pow(p.mass*cLight2,(q-2)/(2*q-5)) / 
							pow(9*zeda*accE,1.0/(2*q-5));
	if (p.id == "ntProton") {
		double sigmapp = 34.3e-27;
		double Emax_pp = accE*electronCharge*B/(0.5*dens*sigmapp);
		return min(min(min(Emax_adv,Emax_pp),min(Emax_Hillas,Emax_syn)),Emax_diff);
	} else
		return min(min(min(Emax_syn,Emax_adv),Emax_Hillas),Emax_diff);
}

double cutOffPL(double E, double Emin, double Emax)
{
	return pow(E,-pIndex)*exp(-E/Emax)*exp(-Emin/E);
}

double auxFun(double g, double gammaMax)
{
	if (pIndex < 2.0) return pow(gammaMax,2.0-pIndex)*pow(g,pIndex+1.0);
	else if (pIndex > 2.0) return g*g*g;
	else return log(gammaMax/g)*g*g;
}
double func(double g, double LHS, double norm_temp, double gammaMax)
{
	return LHS - auxFun(g,gammaMax)*sqrt(g*g-1.0)*exp(-g/norm_temp);
}

double findGammaMin(double temp, double Emax)
{
	double norm_temp = boltzmann*temp / (electronMass*cLight2);
	double aTheta = 3.0 - 6.0/(4.0+5.0*norm_temp); // Gammie & Popham (1998)
	double gammaMax = (pIndex <= 2.0) ? Emax / (electronMass*cLight2) : 1.0;;
	double g1,g2;
	g2 = 100.0;
	double LHS = etaInj*norm_temp*norm_temp*aTheta*gsl_sf_bessel_Kn(2,1.0/norm_temp)*
					(pIndex != 2 ? abs(pIndex-2.0) : 1.0);
	double RHS2 = auxFun(g2,gammaMax)*sqrt(g2*g2-1.0)*exp(-g2/norm_temp);
	while (LHS < RHS2) {
		g2 *= 2.0;
		RHS2 = auxFun(g2,gammaMax)*sqrt(g2*g2-1.0)*exp(-g2/norm_temp);
	}
	g1 = g2/1.1;
	double RHS1 = auxFun(g1,gammaMax)*sqrt(g1*g1-1.0)*exp(-g1/norm_temp);
	while ((LHS-RHS2)*(LHS-RHS1) > 0.0) {
		g1 /= 1.1;
		RHS1 = auxFun(g1,gammaMax)*sqrt(g1*g1-1.0)*exp(-g1/norm_temp);
	}
	
	double g = g1;
	double tol = 1.0e-3;
	while ((g2-g1) >= tol) {
		g = (g1+g2)/2;
		if (func(g,LHS,norm_temp,gammaMax) == 0.0)
			return g;
		else if (func(g,LHS,norm_temp,gammaMax)*func(g1,LHS,norm_temp,gammaMax) < 0.0)
			g2 = g;
		else
			g1 = g;
    }
	return 0.5*(g2+g1);
}

void injection(Particle& p, State& st)
{
	pIndex = GlobalConfig.get<double>("nonThermal.injection.primaryIndex");
    double sumQ = 0.0;
	
	std::ofstream file;
	if (p.id == "ntElectron")
		file.open("gammaMin.txt",ios::out);

	p.ps.iterate([&](const SpaceIterator& iR) {
		if (p.id == "ntProton")
			etaInj = GlobalConfig.get<double>("nonThermal.injection.energyFraction_i");
		else
			etaInj = GlobalConfig.get<double>("nonThermal.injection.energyFraction_e");
			
		const double r = iR.val(DIM_R);
		double rB1 = r/sqrt(paso_r);
		double rB2 = rB1*paso_r;
		const double vol = volume(r);
		
		double gammaMin = 2.0;
		double Emax = eEmax(p,r,st.magf.get(iR),-radialVel(r),st.denf_i.get(iR));

		double Emin = gammaMin*p.mass*cLight2;
		if (p.id == "ntElectron") {
			double temp = st.tempElectrons.get(iR);
			gammaMin = findGammaMin(temp,Emax);
			Emin = gammaMin*electronRestEnergy;
			file << r/schwRadius << "\t" << gammaMin << endl;
		}
		
		//double int_E = RungeKuttaSimple(Emin,Emax,[&Emax,&Emin](double E){
		//	return E*cutOffPL(E,Emin,Emax);});
		double int_E = integSimpson(log(Emin),log(Emax),[Emin,Emax](double loge)
				{
					double e = exp(loge);
					return e*e*cutOffPL(e,Emin,Emax);
				},100);
        
        double norm_temp, dens;
		if(p.id == "ntElectron") {
			norm_temp = boltzmann*st.tempElectrons.get(iR)/(p.mass*cLight2);
			dens = st.denf_e.get(iR);
		} else if(p.id == "ntProton") {
			norm_temp = boltzmann*st.tempIons.get(iR)/(p.mass*cLight2);
			dens = st.denf_i.get(iR);
		}
		
		double aTheta = 3.0 - 6.0/(4.0+5.0*norm_temp); // Gammie & Popham (1998)
		double uth = dens*norm_temp*(p.mass*cLight2)*aTheta;   // erg cm^-3
		double area = 4.0*pi*r*height_fun(r);
		double Q0 = etaInj * uth * area/vol * (-radialVel(r));   // power injected in nt particles [erg cm^-3 s^-1]
		double Q0p = Q0/int_E;
		
		p.ps.iterate([&](const SpaceIterator& iRE) {
			const double E = iRE.val(DIM_E);
			double total = cutOffPL(E,Emin,Emax)*Q0p;
			p.injection.set(iRE,total); //en unidades de erg^-1 s^-1 cm^-3
		},{-1,iR.coord[DIM_R],0});
		sumQ += Q0*vol;
	},{0,-1,0});
	file.close();
	cout << "Total power injected in " << p.id << " = " << sumQ << endl; 
}

void injectionChargedPion(Particle& p, State& st)
{
	show_message(msgStart,Module_pionInjection);
	
	double sumTot = 0.0;
	p.ps.iterate([&](const SpaceIterator& iR) {
		p.ps.iterate([&](const SpaceIterator& iRE) {
			double Epi = iRE.val(DIM_E);
			double injection = ppPionInj(Epi,st.ntProton,st.denf_i.get(iRE),iRE);
			injection += pgammaPionInj(Epi,st.ntProton,st.photon.distribution,iRE,
										st.photon.emin(),st.photon.emax());
			p.injection.set(iRE,injection);
			double tdecay = chargedPionMeanLife*(Epi/(chargedPionMass*cLight2));
			double ratecool = lossesHadronics(Epi,st.denf_i.get(iR),p)/Epi;
			p.distribution.set(iRE,injection*pow(ratecool+pow(tdecay,-1),-1));
		},{-1,iR.coord[DIM_R],0});
		
		sumTot += integSimpson(log(p.emin()),log(p.emax()),[&](double loge)
					{
						double e = exp(loge);
						return e*e*p.injection.interpolate({{DIM_E,e}},&iR.coord);
					},100)*volume(iR.val(DIM_R));
	},{0,-1,0});
	cout << "Total power injected in " << p.id << " = " << sumTot << "\t" << endl; 
	writeEandRParamSpace("pionInjection",p.injection,0,1);
	writeEandRParamSpace("pionDistribution",p.distribution,0,1);
	show_message(msgEnd,Module_pionInjection);
}

void injectionMuon(Particle& p, State& st)
{
    double sumTot = 0.0;
	show_message(msgStart, Module_muonInjection);
	p.ps.iterate([&](const SpaceIterator& iR) {
		p.ps.iterate([&](const SpaceIterator& iRE) {
			const double Emu = iRE.val(DIM_E);
			double injection = muonInjNew(Emu,p,st.ntChargedPion,iRE);
			p.injection.set(iRE,injection);
			double ratecool = lossesSyn(Emu,st.magf.get(iR),p)/Emu;
			double tdecay = muonMeanLife*(Emu/(muonMass*cLight2));
			p.distribution.set(iRE,injection*pow(ratecool+pow(tdecay,-1),-1));
		},{-1,iR.coord[DIM_R],0});
		
		sumTot += integSimpson(log(p.emin()),log(p.emax()),[&](double loge)
					{
						double e = exp(loge);
						return e*e*p.injection.interpolate({{DIM_E,e}},&iR.coord);
					},100)*volume(iR.val(DIM_R));
		
	},{0,-1,0});
	cout << "Total power injected in " << p.id << " = " << sumTot << "\t" << endl;
	writeEandRParamSpace("muonInjection",p.injection,0,1);
	writeEandRParamSpace("muonDistribution",p.distribution,0,1);
	show_message(msgEnd, Module_muonInjection);
}

void injectionPair(Particle& p, State& st)
{
	show_message(msgStart, Module_pairInjection);
	double sumTot = 0.0;
	double sumTotPh = 0.0;
	p.ps.iterate([&](const SpaceIterator& iR) {
		p.ps.iterate([&](const SpaceIterator& iRE) {
			double Eee = iRE.val(DIM_E);
			double injection = 	pairMuonDecayNew(Eee,st.ntMuon,iRE)
							+	pairBH(Eee,st.ntProton,st.photon.distribution,iRE,st.photon.emin(),st.photon.emax())
							+	pairGammaGammaNew(Eee,st.ntPhoton.distribution,st.photon.distribution,
								iRE,st.photon.emin(),st.photon.emax(),st.ntPhoton.emin(),st.ntPhoton.emax());
			p.injection.set(iRE,injection);
		},{-1,iR.coord[DIM_R],0});
		sumTot += integSimpson(log(p.emin()),log(p.emax()),[&](double loge)
					{
						double e = exp(loge);
						return e*e*p.injection.interpolate({{DIM_E,e}},&iR.coord);
					},100)*volume(iR.val(DIM_R));
		sumTotPh += integSimpson(log(st.ntPhoton.emin()),log(st.ntPhoton.emax()),
					[&st,&iR](double logEg)
					{
						double Eg = exp(logEg);
						double ntPh = st.ntPhoton.distribution.interpolate({{DIM_E,Eg}},&iR.coord);
						double rategg = cLight*integSimpson(log(st.photon.emin()),log(st.photon.emax()),
								[Eg,&st,&iR](double logEph)
								{
									double Eph = exp(logEph);
									double nPh = st.photon.distribution.interpolate({{DIM_E,Eph}},&iR.coord);
									return ggCrossSection2(Eg,Eph)*nPh*Eph;
								},30);
						return Eg*Eg*ntPh*rategg;
					},100)*volume(iR.val(DIM_R));
	},{0,-1,0});
	cout << "Total photon power absorbed to create pairs = " << sumTotPh << "\t" << endl;
	cout << "Total power injected in " << p.id << " = " << sumTot << "\t" << endl;
	writeEandRParamSpace("secondaryPairInjection",p.injection,0,1);
	writeEandRParamSpace("secondaryPairDistribution",p.distribution,0,1); 
	show_message(msgEnd, Module_pairInjection);
}
