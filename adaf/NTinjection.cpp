#include "NTinjection.h"
#include "write.h"
#include "messages.h"
#include "globalVariables.h"
#include "adafFunctions.h"
#include "absorption.h"
#include "NTdistribution.h"
#include <finjection/ppPionInj.h>
#include <finjection/pgammaPionInj.h>
#include <finjection/muonInj.h>
#include <finjection/neutrinoInj.h>
#include <flosses/lossesHadronics.h>
#include <flosses/lossesSyn.h>
#include <finjection/pairInjectionExact.h>
#include <finjection/pairgg.h>
#include <finjection/pairMuonDecay.h>
#include <finjection/pairBH.h>
#include <gsl/gsl_sf_bessel.h>
#include "losses.h"
#include <flosses/nonThermalLosses.h>

#include <fparameters/parameters.h>
#include <fparameters/SpaceIterator.h>
#include <fparameters/Dimension.h>
#include <fmath/fbisection.h>
#include <fmath/RungeKutta.h>
#include <fmath/physics.h>
#include <boost/property_tree/ptree.hpp>
#include <iostream>
#include <fstream>

double eEmax(Particle& p, double r, double B, double v, double dens, int jR)
{
	double accE = GlobalConfig.get<double>("nonThermal.injection.PL.accEfficiency");
	double dr = r*(paso_r-1.0);
	double h = height_fun(r);
	double sumAdvTime = 0.0;
	for (int jj=nR-1;jj>=jR;jj--) {
		double rAux = p.ps[DIM_R][jj];
		double Baux = magneticField(rAux);
		double drAux = rAux*(sqrt(paso_r)-1.0/sqrt(paso_r));
		double tCellAux = drAux / abs(radialVel(rAux));
		sumAdvTime += Baux * tCellAux;
	}

    double Emax_adv = accE*cLight*electronCharge*sumAdvTime; 
    double Emax_syn = p.mass*cLight2*sqrt(accE*6.0*pi*electronCharge / (thomson*B)) * p.mass/electronMass;
    double Emax_Hillas = electronCharge*B*h;
	double zeda = GlobalConfig.get<double>("nonThermal.injection.SDA.fractionTurbulent");
	double q = GlobalConfig.get<double>("nonThermal.injection.SDA.powerSpectrumIndex");
	double Emax_diff = pow(h*electronCharge*B,(q-3)/(2*q-5)) * pow(p.mass*cLight2,(q-2)/(2*q-5)) / 
							pow(9*zeda*accE,1.0/(2*q-5));
	if (p.id == "ntProton") {
		double sigmapp = 34.3e-27;
		double Emax_pp = accE*electronCharge*B/(0.5*dens*sigmapp);
		return min(min(Emax_adv,Emax_pp),min(Emax_Hillas,Emax_syn));
	} else
		return min(min(min(Emax_syn,Emax_adv),Emax_Hillas),Emax_diff);
}

double eEmax_numerical(Particle& p, double r, double B, double v, double dens, int jR, State& st,
						SpaceCoord i)
{
	double accE = GlobalConfig.get<double>("nonThermal.injection.PL.accEfficiency");
	double dr = r*(paso_r-1.0);
	double h = height_fun(r);
	double sumAdvTime = 0.0;
	for (int jj=nR-1;jj>=jR;jj--) {
		double rAux = p.ps[DIM_R][jj];
		double Baux = magneticField(rAux);
		double drAux = rAux*(sqrt(paso_r)-1.0/sqrt(paso_r));
		double tCellAux = drAux / abs(radialVel(rAux));
		sumAdvTime += Baux * tCellAux;
	}
    double Emax_adv = accE*cLight*electronCharge*sumAdvTime;
	double Emax_Hillas = electronCharge*B*h;
							
	fun1 rate_cool_acc = [&p,&st,&i,accE,B] (double E)
						{
							double rate_cool = losses(E,p,st,i)/E;
							double rate_acc = accE*cLight*electronCharge*B/E;
							return rate_cool - rate_acc;
						};
	double Emax_cool = 10*p.emin()*cLight2;
	Bisect(rate_cool_acc,2.0*p.emin(),p.emax()/2.0,Emax_cool);
	
	return min(Emax_cool, min(Emax_Hillas,Emax_adv));
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
	
	if (p.id == "ntProton")
		etaInj = GlobalConfig.get<double>("nonThermal.injection.energyFraction_i");
	else if (p.id == "ntElectron")
		etaInj = GlobalConfig.get<double>("nonThermal.injection.energyFraction_e");
	
	std::ofstream file;
	if (p.id == "ntElectron")
		file.open("gammaMin.txt",ios::out);


	// CALCULATION OF THE NORMALIZATION FACTOR
	double sum = 0.0;
	Vector Qfactor(nR,1.0);
	p.ps.iterate([&](const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		double vol = volume(r);
		double norm_temp, dens;
		if (p.id == "ntElectron") {
			norm_temp = boltzmann*st.tempElectrons.get(iR)/(p.mass*cLight2);
			dens = st.denf_e.get(iR);
		} else if (p.id == "ntProton") {
			norm_temp = boltzmann*st.tempIons.get(iR)/(p.mass*cLight2);
			dens = st.denf_i.get(iR);
		}
		double aTheta = 3.0 - 6.0/(4.0+5.0*norm_temp); // Gammie & Popham (1998)
		double uth = aTheta * dens * p.mass*cLight2 * norm_temp;
		double vA = st.magf.get(iR)/sqrt(4.0*pi*massDensityADAF(r));
		double h = height_fun(r);
		Qfactor[iR.coord[DIM_R]] = (accMethod == 0) ? uth * abs(radialVel(r)) / r :
										dens * p.mass*cLight2 * cLight/schwRadius;
		sum += vol * Qfactor[iR.coord[DIM_R]];
	},{0,-1,0});
	Ainjection = etaInj*accRateOut*cLight2 / sum;		// [non-dimensional]
	
////////////////////////////////////////////////////////

	int jR = 0;
	p.ps.iterate([&](const SpaceIterator& iR) {
		
		const double r = iR.val(DIM_R);
		double rB1 = r/sqrt(paso_r);
		double rB2 = rB1*paso_r;
		const double vol = volume(r);
		
		double gammaMin = 6.0;
		double Emax = eEmax_numerical(p,r,st.magf.get(iR),-radialVel(r),st.denf_i.get(iR),jR,st,iR.coord);
		double Emin = gammaMin*p.mass*cLight2;
		if (p.id == "ntElectron") {
			double temp = st.tempElectrons.get(iR);
			gammaMin = findGammaMin(temp,Emax);
			Emin = gammaMin*electronRestEnergy;
			file << r/schwRadius << "\t" << gammaMin << endl;
		}
		
		double int_E = integSimpsonLog(p.emin(),p.emax(),[Emin,Emax](double e)
				{
					return e*cutOffPL(e,Emin,Emax);
				},100);
        
		double Q0p = Ainjection * Qfactor[iR.coord[DIM_R]] / int_E;
		
		p.ps.iterate([&](const SpaceIterator& iRE) {
			const double E = iRE.val(DIM_E);
			double total = cutOffPL(E,Emin,Emax)*Q0p;
			p.injection.set(iRE,total); //en unidades de erg^-1 s^-1 cm^-3
		},{-1,iR.coord[DIM_R],0});
		jR++;
	},{0,-1,0});
	file.close();
	cout << "etaInj Mdot c^2 = " << etaInj*accRateOut*cLight2 << endl;
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
			//p.distribution.set(iRE,injection*pow(ratecool+pow(tdecay,-1),-1));
		},{-1,iR.coord[DIM_R],0});
		
		sumTot += integSimpson(log(p.emin()),log(p.emax()),[&](double loge)
					{
						double e = exp(loge);
						return e*e*p.injection.interpolate({{DIM_E,e}},&iR.coord);
					},100)*volume(iR.val(DIM_R));
	},{0,-1,0});
	cout << "Total power injected in " << p.id << " = " << sumTot << "\t" << endl; 
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
			//p.distribution.set(iRE,injection*pow(ratecool+pow(tdecay,-1),-1));
		},{-1,iR.coord[DIM_R],0});
		
		sumTot += integSimpson(log(p.emin()),log(p.emax()),[&](double loge)
					{
						double e = exp(loge);
						return e*e*p.injection.interpolate({{DIM_E,e}},&iR.coord);
					},100)*volume(iR.val(DIM_R));
		
	},{0,-1,0});
	cout << "Total power injected in " << p.id << " = " << sumTot << "\t" << endl;
	writeEandRParamSpace("muonInjection",p.injection,0,1);
	show_message(msgEnd, Module_muonInjection);
}

void injectionPair(Particle& p, State& st, int it)
{
	show_message(msgStart, Module_pairInjection);
	
	std::ofstream secondariesMuonFile, secondariesBHFile, secondariesGammaGammaFile;
	
	if (it == 1) {
		secondariesMuonFile.open("secondariesInjectionMuon.dat", std::ios::out);
		secondariesBHFile.open("secondariesInjectionBH.dat", std::ios::out);
		secondariesGammaGammaFile.open("secondariesInjectionGammaGamma.dat", std::ios::out);
	}
	
	double sumTot = 0.0;
	double sumTotPh = 0.0;
	p.ps.iterate([&](const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		double vol = volume(r);
		double height = height_fun(iR.val(DIM_R));
		p.ps.iterate([&](const SpaceIterator& iRE) {
			double Eee = iRE.val(DIM_E);
			if (it == 1) {
				double injectionMuon = pairMuonDecayNew(Eee,st.ntMuon,iRE);
				double injectionBH = pairBH(Eee,st.ntProton,st.photon.distribution,iRE,
								st.photon.emin(),st.photon.emax());
				p.injection.set(iRE, injectionMuon + injectionBH);
				secondariesMuonFile << safeLog10(Eee/EV_TO_ERG) << "\t" << iR.coord[DIM_R] << "\t"
								<< r/schwRadius << "\t" << injectionMuon*vol << endl;
				secondariesBHFile << safeLog10(Eee/EV_TO_ERG) << "\t" << iR.coord[DIM_R] << "\t"
								<< r/schwRadius << "\t" << injectionBH*vol << endl;
			}
			double injectionGammaGamma = pairGammaGammaNew(Eee,st.ntPhoton.distribution,st.photon.distribution,
						iRE,st.photon.emin(),st.photon.emax(),st.ntPhoton.emin(),st.ntPhoton.emax());
			//double tcool = Eee/lossesSyn(Eee,st.magf.get(iR),p);
			p.injection.set(iRE,p.injection.get(iRE)+injectionGammaGamma);
			secondariesGammaGammaFile << safeLog10(Eee/EV_TO_ERG) << "\t" << vol << "\t"
								<< r/schwRadius << "\t" << injectionGammaGamma << endl;
		},{-1,iR.coord[DIM_R],0});
		sumTot += integSimpsonLog(p.emin(),p.emax(),[&st,&p,iR](double e)
					{
						//double tcool = e/lossesSyn(e,st.magf.get(iR),p);
						double inj = p.injection.interpolate({{DIM_E,e}},&iR.coord);
						return e*inj;
					},100)*volume(iR.val(DIM_R));
		sumTotPh += integSimpsonLog(st.ntPhoton.emin(),st.ntPhoton.emax(),
					[&st,&iR,height](double Eg)
					{
						double ntPh = st.ntPhoton.distribution.interpolate({{DIM_E,Eg}},&iR.coord);
						double kappa_gg = 2.0/sqrt(pi)*st.tau_gg.interpolate({{DIM_E,Eg}},&iR.coord)/height;
						double rategg = cLight*kappa_gg;
						return Eg*ntPh*rategg;
					},100)*volume(iR.val(DIM_R));
	},{0,-1,0});
	cout << "Total photon power absorbed to create pairs = " << sumTotPh << "\t" << endl;
	cout << "Total power injected in " << p.id << " = " << sumTot << "\t" << endl;
	writeEandRParamSpace("secondaryPairInjection",p.injection,0,1);
	
	if (it == 1) {
		secondariesMuonFile.close();
		secondariesBHFile.close();
	}
	secondariesGammaGammaFile.close();
	
	show_message(msgEnd, Module_pairInjection);
}

void injectionNeutrino(Particle& nu, State& st)
{
    double sumTot = 0.0;
	show_message(msgStart,Module_muonInjection);
	nu.ps.iterate([&](const SpaceIterator& iR) {
		nu.ps.iterate([&](const SpaceIterator& iRE) {
			const double Enu = iRE.val(DIM_E);
			double injection1 = muonNeutrinoInjection(Enu,nu,st.ntMuon,st.ntChargedPion,iRE);
			double injection2 = electronNeutrinoInjection(Enu,nu,st.ntMuon,iRE);
			nu.injection.set(iRE,injection1+injection2);
		},{-1,iR.coord[DIM_R],0});
		sumTot += integSimpsonLog(nu.emin(),nu.emax(),[&](double e)
					{
						return e*nu.injection.interpolate({{DIM_E,e}},&iR.coord);
					},100)*volume(iR.val(DIM_R));
		
	},{0,-1,0});
	cout << "Total power injected in " << nu.id << " = " << sumTot << "\t" << endl;
	writeEandRParamSpace("neutrinoInjection",nu.injection,0,1);
	show_message(msgEnd, Module_neutrinoInjection);
}