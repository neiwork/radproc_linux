#include "NTinjection.h"
#include "write.h"
#include "messages.h"
#include "globalVariables.h"
#include "adafFunctions.h"
#include "absorption.h"
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
#include <fmath/tridiagonal.h>
#include <boost/property_tree/ptree.hpp>
#include <iostream>
#include <fstream>

#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>

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
	double dr = r*(1.0-1.0/paso_r);
	double h = height_fun(r);
	double rho = massDensityADAF(r);
	double vA = B / sqrt(4*pi*rho);
	
	double zeda = GlobalConfig.get<double>("nonThermal.injection.SDA.fractionTurbulent");
	double q = GlobalConfig.get<double>("nonThermal.injection.SDA.powerSpectrumIndex");
	double sumAdvTime = (accMethod == 0) ?
	  integSimpsonLog(r,p.ps[DIM_R][nR-1],[&p,zeda,q] (double rAux)
		{
			double Baux = magneticField(rAux);
			double vR = radialVel(rAux);
			return Baux/abs(vR);
		},50)*accE*electronCharge*cLight
	: integSimpsonLog(r,p.ps[DIM_R][nR-1],[&p,zeda,q] (double rAux)
		{
			double Baux = magneticField(rAux);
			double height_aux = height_fun(rAux);
			double rhoAux = massDensityADAF(rAux);
			double vAlfv = Baux / sqrt(4.0*pi*rhoAux);
			double rL = p.mass*cLight2/(electronCharge*Baux);
			double kMin = 1.0/height_aux;
			double factor = zeda * (cLight*kMin) * P2(vAlfv/cLight) * pow(rL*kMin,q-2.0);
			return factor / abs(radialVel(rAux));
		},50);
	/*
	B = magneticField(r);
	double height = height_fun(r);
	double rL = p.mass*cLight2/(electronCharge*B);
	double kMin = 1.0/height;
	double factor = zeda * (cLight*kMin) * P2(vA/cLight) * pow(rL*kMin,q-2.0);
	sumAdvTime = factor / abs(radialVel(r)) * r; */
	double rL = p.mass*cLight2/(electronCharge*B);
	double gamma_min = 2.0*p.emin() / (p.mass*cLight2);
	double Emax_adv = (accMethod == 0) ? sumAdvTime :
				p.mass*cLight2 * (pow((4.0-q*q)*sumAdvTime + pow(gamma_min,2.0-q),1.0/(2.0-q)));
	double Emax_Hillas = electronCharge*B*h * (vA/cLight);
							
	fun1 rate_cool_acc = [&p,&st,&i,accE,B,h,rho,r] (double E)
						{
							double rate_cool = losses(E,p,st,i)/E;
							double rate_diff = 1.0/diffusionTimeTurbulence(E,h,p,B);
							double rate_acc = (accMethod == 0) ? accelerationRate(E,B) : 
												accelerationRateSDA(E,p,B,h,rho);
							return rate_cool + rate_diff - rate_acc;
						};
	double Emax_cool = 10*p.emin()*cLight2;
	Bisect(rate_cool_acc,2.0*p.emin(),p.emax()/2.0,Emax_cool);
	return min(Emax_cool, min(Emax_Hillas,Emax_adv));
}


double cutOffPL(double E, double Emin, double Emax)
{
	return pow(E,-pIndex)*exp(-E/Emax)*exp(-Emin/E);
}

double SDA_PL(double E, double Emin, double Emax)
{
	return pow(E,0.5)*exp(-pow(27*E/Emax,1.0/3.0))*exp(-Emin/E);
}

double SDA_exact(double E, double Emin, double Emax, Particle& p, double r)
{
	double gmax2 = Emax / (p.mass*cLight2);
	double zeda = GlobalConfig.get<double>("nonThermal.injection.SDA.fractionTurbulent");
	double q = GlobalConfig.get<double>("nonThermal.injection.SDA.powerSpectrumIndex");
	double B = magneticField(r);
	double height = height_fun(r);
	double vA = B / sqrt(4.0*pi*massDensityADAF(r));
	double rL = p.mass*cLight2 / (electronCharge*B);
	double theta = pow(rL/height,4.0-2.0*q) / P2(3*vA/cLight * zeda);
	double g = E / (p.mass*cLight2);
	double ginj = Emin / (p.mass*cLight2);
	double gmin = min(ginj,g);
	double gmax = max(ginj,g);
	double beta = 3.0/(2.0-q);
	double arg1 = sqrt(theta)*pow(gmin,2.0-q)/(2.0-q);
	double arg2 = sqrt(theta)*pow(gmax,2.0-q)/(2.0-q);
	double Ifun = gsl_sf_bessel_Inu((beta-1.0)/2.0,arg1);
	double Kfun = (arg2 < 500.0) ? gsl_sf_bessel_Knu((beta-1.0)/2.0,arg2) : exp(-arg2) * sqrt(pi/(2.0*arg2));
	double gmax1 = pow(theta,-0.5/(2.0-q));
	return sqrt(g/ginj)*pow(g*ginj,(2.0-q)/2.0)*Ifun*Kfun; 
	//* (gmax2 < gmax1 && g > gmax2 ? exp(3) * exp(-pow(27*g/gmax2,1.0/3.0)) : 1.0);
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
		Qfactor[iR.coord[DIM_R]] = (accMethod == 0 ? (r/schwRadius < 200.0 ? uth * st.magf.get(iR) : 0.0) : //	 * abs(radialVel(r)) / r :
										dens * p.mass*cLight2 * cLight/schwRadius * P2(vA/cLight) );
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
		
		double zeda = GlobalConfig.get<double>("nonThermal.injection.SDA.fractionTurbulent");
		double q = GlobalConfig.get<double>("nonThermal.injection.SDA.powerSpectrumIndex");
		double B = magneticField(r);
		double height = height_fun(r);
		double vA = B / sqrt(4.0*pi*massDensityADAF(r));
		double rL = p.mass*cLight2 / (electronCharge*B);
		double theta = pow(rL/height,4.0-2.0*q) / P2(3*vA/cLight * zeda);
		
		double gammaMin = 2.0;
		double Emax = eEmax_numerical(p,r,st.magf.get(iR),-radialVel(r),st.denf_i.get(iR),jR,st,iR.coord);
		double Emin = gammaMin*p.mass*cLight2;
		if (p.id == "ntElectron") {
			double temp = st.tempElectrons.get(iR);
			gammaMin = findGammaMin(temp,Emax);
			Emin = gammaMin*electronRestEnergy;
			file << r/schwRadius << "\t" << gammaMin << endl;
		}
		double int_E = integSimpsonLog(p.emin(),p.emax(),[Emin,Emax,r,&p](double e)
				{
					if (accMethod == 0)
						return e*cutOffPL(e,Emin,Emax);
					else
						return e*SDA_exact(e,Emin,Emax,p,r);
				},100);
        
		double Q0p = Ainjection * Qfactor[iR.coord[DIM_R]] / int_E;
		
		p.ps.iterate([&](const SpaceIterator& iRE) {
			const double E = iRE.val(DIM_E);
			double total = 0.0;
			if (accMethod == 0)
				total = cutOffPL(E,Emin,Emax)*Q0p;
			else
				total = SDA_exact(E,Emin,Emax,p,r)*Q0p;
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
	
	double sumTotp = 0.0;
	double sumTotpp = 0.0;
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
		sumTotp += integSimpsonLog(st.ntProton.emin(),st.ntProton.emax(),[&](double ep)
					{
						return ep*st.ntProton.injection.interpolate({{DIM_E,ep}},&iR.coord);
					},100)*volume(iR.val(DIM_R));
		sumTotpp += integSimpsonLog(st.ntProton.emin(),st.ntProton.emax(),[&](double ep)
					{
						return ep*st.ntProton.distribution.interpolate({{DIM_E,ep}},&iR.coord)
						* lossesHadronics(ep,st.denf_i.get(iR),st.ntProton)/ep ;
					},100)*volume(iR.val(DIM_R));
		sumTot += integSimpsonLog(p.emin(),p.emax(),[&](double e)
					{
						return e*p.injection.interpolate({{DIM_E,e}},&iR.coord);
					},100)*volume(iR.val(DIM_R));
	},{0,-1,0});
	cout << "Total power injected in " << st.ntProton.id << " = " << sumTotp << "\t" << endl; 
	cout << "Total power in pp interactions" << " = " << sumTotpp << "\t" << endl; 
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
	
	//p.ps.iterate([&](const SpaceIterator& iR) {
	double sumTot = 0.0;
	double sumTotPh = 0.0;
	#pragma omp parallel for
	for (size_t jR=0;jR<nR;jR++) {
		SpaceCoord iR = {0,jR,0};
		double r = st.denf_e.ps[DIM_R][jR];
		//double r = iR.val(DIM_R);
		double vol = volume(r);
		double height = height_fun(r);
		p.ps.iterate([&](const SpaceIterator& iRE) {
			double Eee = iRE.val(DIM_E);
			if (it == 1) {
				double injectionMuon = pairMuonDecayNew(Eee,st.ntMuon,iR);
				double injectionBH = pairBH(Eee,st.ntProton,st.photon.distribution,iR,
								st.photon.emin(),st.photon.emax());
				p.injection.set(iRE, injectionMuon + injectionBH);
				secondariesMuonFile << safeLog10(Eee/EV_TO_ERG) << "\t" << iRE.coord[DIM_R] << "\t"
								<< r/schwRadius << "\t" << injectionMuon*vol << endl;
				secondariesBHFile << safeLog10(Eee/EV_TO_ERG) << "\t" << iRE.coord[DIM_R] << "\t"
								<< r/schwRadius << "\t" << injectionBH*vol << endl;
			}
			double injectionGammaGamma = pairGammaGammaNew(Eee,st.ntPhoton.distribution,st.photon.distribution,
						iRE,st.photon.emin(),st.photon.emax(),st.ntPhoton.emin(),st.ntPhoton.emax());
			
			//double tcool = Eee/lossesSyn(Eee,st.magf.get(iR),p);
			p.injection.set(iRE,p.injection.get(iRE)+injectionGammaGamma);
			secondariesGammaGammaFile << safeLog10(Eee/EV_TO_ERG) << "\t" << vol << "\t"
								<< r/schwRadius << "\t" << injectionGammaGamma << endl;
		},{-1,jR,0});
		sumTot += integSimpsonLog(p.emin(),p.emax(),[&st,&p,iR](double e)
					{
						//double tcool = e/lossesSyn(e,st.magf.get(iR),p);
						double inj = p.injection.interpolate({{DIM_E,e}},&iR);
						return e*inj;
					},100)*volume(r);
		sumTotPh += integSimpsonLog(st.ntPhoton.emin(),st.ntPhoton.emax(),
					[&st,&iR,height](double Eg)
					{
						double ntPh = st.ntPhoton.distribution.interpolate({{DIM_E,Eg}},&iR);
						double kappa_gg = 2.0/sqrt(pi)*st.tau_gg.interpolate({{DIM_E,Eg}},&iR)/height;
						double rategg = cLight*kappa_gg;
						return Eg*ntPh*rategg;
					},100)*volume(r);
	};
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


void injectionFokkerPlanckOneZone(Particle& p, State& st)
// This routine solves the transport equation by considering separated one zones at each
// radius in the RIAF. The escape rate is the sum of the accretion rate (1/tacc) and the
// diffusion rate. It considers stochastic acceleration by turbulence at each zone, where 
// the steady transport equation is an ordinary Fokker-Plack diff. equation in the energy.
// It is solved by finite differences (CC70, PP96).
// Problems: It does not treat re-acceleration.
{
	
	// CALCULATION OF THE NORMALIZATION FACTOR
	if (p.id == "ntProton")
		etaInj = GlobalConfig.get<double>("nonThermal.injection.energyFraction_i");
	else if (p.id == "ntElectron")
		etaInj = GlobalConfig.get<double>("nonThermal.injection.energyFraction_e");
	
	double sum = 0.0;
	Vector Qfactor(nR,1.0);
	p.ps.iterate([&](const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		double vol = volume(r);
		double dens = 0.0;
		double temp = 0.0;
		if (p.id == "ntElectron") {
			dens = st.denf_e.get(iR);
			temp = st.tempElectrons.get(iR);
		} else if (p.id == "ntProton") {
			dens = st.denf_i.get(iR);
			temp = st.tempIons.get(iR);
		}
		Qfactor[iR.coord[DIM_R]] = dens * boltzmann * temp * abs(radialVel(r))/r;
		Qfactor[iR.coord[DIM_R]] = dens;
		sum += vol * Qfactor[iR.coord[DIM_R]];
	},{0,-1,0});
	Ainjection = etaInj*accRateOut*cLight2 / sum;		// [non-dimensional]
	
////////////////////////////////////////////////////////

	p.ps.iterate([&](const SpaceIterator& iR) {
		
		double r = iR.val(DIM_R);
		double rB2 = r*sqrt(paso_r);
		double rB1 = rB2/paso_r;
		double vR = radialVel(r);
		double tCell = (rB2-rB1)/abs(vR);
		double vol = volume(r);
		double height = height_fun(r);
		double B = st.magf.get(iR);
		double rho = massDensityADAF(r);
		double tAccretion = accretionTime(r);
		double rateAccretion = 1.0/tAccretion;
		double rateWind = 2.0*s*abs(radialVel(r))/r;
		
		// We define a new mesh of points
		size_t M = 199;
		double g_min = p.emin()/(p.mass*cLight2);
		double g_max = p.emax()/(p.mass*cLight2);
		
		double paso_g = pow(g_max/g_min,1.0/M);
		Vector g_au(M+3,g_min/paso_g);
		Vector delta_g_m_au(M+2,0.0);
		Vector delta_g(M+1,0.0);
		
		double g_inj = 2.0*g_min;
		double E_inj = g_inj*p.mass*cLight2;
		
		double dens = 0.0;
		if (p.id == "ntProton")
			dens = st.denf_i.get(iR);
		else if (p.id == "ntElectron")
			dens = st.denf_e.get(iR);
		
		double Q0 = Ainjection * Qfactor[iR.coord[DIM_R]] / g_inj;
		double sigma = g_inj / 5.0;
		
		for (size_t jg=1;jg<M+3;jg++)	g_au[jg] = g_au[jg-1]*paso_g;
		for (size_t jg=0;jg<M+2;jg++)	delta_g_m_au[jg] = g_au[jg+1]-g_au[jg];
		for (size_t jg=0;jg<M+1;jg++)	delta_g[jg] = 0.5*(g_au[jg+2]-g_au[jg]);
		
		fun1 Bfun = [&p,&st,&iR,B,rho,height,vR,r] (double g) 
						{
							double E = g*p.mass*cLight2;
							double tCool = E/abs(losses(E,p,st,iR));
							double Dg = diffCoeff_g(g,p,height,B,rho);
							return g/tCool - 2.0*Dg/g - 0.5*g*abs(vR)/r;
						};

		fun1 Cfun = [&p,B,height,rho] (double g)
						{
							double Dg = diffCoeff_g(g,p,height,B,rho);
							return Dg;
						};

		fun1 Tfun = [r,&p,height,B,rateAccretion,&rateWind] (double g) 
						{
							double E = g*p.mass*cLight2;
							double rateDiff = 1.0/diffusionTimeTurbulence(E,height,p,B);
							return pow(rateDiff,-1);
						};
		
		fun1 Qfun = [Q0,g_inj,sigma] (double g)
						{
							return Q0*gsl_ran_gaussian_pdf(g-g_inj,sigma);
						};

		Vector d(M+1,0.0);
		Vector a(M+1,0.0), b(M+1,0.0), c(M+1,0.0);
		for (size_t m=0;m<M+1;m++) {
			double Bm_minusHalf = 0.5 * (Bfun(g_au[m]) + Bfun(g_au[m+1]));
			double Bm_plusHalf = 0.5 * (Bfun(g_au[m+1]) + Bfun(g_au[m+2]));
			double Cm_minusHalf = 0.5 * (Cfun(g_au[m]) + Cfun(g_au[m+1]));
			double Cm_plusHalf = 0.5 * (Cfun(g_au[m+1]) + Cfun(g_au[m+2]));
			double Qm = Qfun(g_au[m+1]);
			double Tm = Tfun(g_au[m+1]);
			double wm_minusHalf = Bm_minusHalf/Cm_minusHalf * delta_g_m_au[m];
			double wm_plusHalf = Bm_plusHalf/Cm_plusHalf * delta_g_m_au[m+1];
			double Wm_minusHalf(0.0), Wm_plusHalf(0.0);
			double Wplus_m_minusHalf(0.0), Wminus_m_minusHalf(0.0);
			double Wplus_m_plusHalf(0.0), Wminus_m_plusHalf(0.0);
			
			if ( abs(wm_minusHalf) < 1.0e-3 ) {
				Wm_minusHalf = pow(1.0+gsl_pow_2(wm_minusHalf)/24+gsl_pow_4(wm_minusHalf)/1920,-1);
				Wplus_m_minusHalf = Wm_minusHalf * exp(0.5*wm_minusHalf);
				Wminus_m_minusHalf = Wm_minusHalf * exp(-0.5*wm_minusHalf);
			} else {
				Wm_minusHalf = abs(wm_minusHalf)*exp(-0.5*abs(wm_minusHalf)) /
								(1.0-exp(-abs(wm_minusHalf)));
				if (wm_minusHalf > 0.0) {
					Wplus_m_minusHalf = abs(wm_minusHalf) / (1.0-exp(-abs(wm_minusHalf)));
					Wminus_m_minusHalf = Wm_minusHalf * exp(-0.5*wm_minusHalf);
				} else {
					Wplus_m_minusHalf = Wm_minusHalf * exp(0.5*wm_minusHalf);
					Wminus_m_minusHalf = abs(wm_minusHalf) / (1.0-exp(-abs(wm_minusHalf)));
				}
			}
			if ( abs(wm_plusHalf) < 1.0e-3 ) {
				Wm_plusHalf = pow(1.0+gsl_pow_2(wm_plusHalf)/24+gsl_pow_4(wm_plusHalf)/1920,-1);
				Wplus_m_plusHalf = Wm_plusHalf * exp(0.5*wm_plusHalf);
				Wminus_m_plusHalf = Wm_plusHalf * exp(-0.5*wm_plusHalf);
			} else {
				Wm_plusHalf = abs(wm_plusHalf)*exp(-0.5*abs(wm_plusHalf)) /
								(1.0-exp(-abs(wm_plusHalf)));
				if (wm_plusHalf > 0.0) {
					Wplus_m_plusHalf = abs(wm_plusHalf) / (1.0-exp(-abs(wm_plusHalf)));
					Wminus_m_plusHalf = Wm_plusHalf * exp(-0.5*wm_plusHalf);
				} else {
					Wplus_m_plusHalf = Wm_plusHalf * exp(0.5*wm_plusHalf);
					Wminus_m_plusHalf = abs(wm_plusHalf) / (1.0-exp(-abs(wm_plusHalf)));
				}
			}
			
			if (m >= 1)
				a[m] = - Cm_minusHalf * Wminus_m_minusHalf / ( delta_g[m] * delta_g_m_au[m] );
			if (m <= M-1)
				c[m] = - Cm_plusHalf * Wplus_m_plusHalf / ( delta_g[m] * delta_g_m_au[m+1] );
			b[m] = 1.0/delta_g[m] * ( Cm_minusHalf * Wplus_m_minusHalf / delta_g_m_au[m] +
						Cm_plusHalf * Wminus_m_plusHalf / delta_g_m_au[m+1] ) + 1.0/Tm;
			d[m] = Qm;
		}
		TriDiagSys(a,b,c,d,M);
		size_t m = 1;
		p.ps.iterate([&](const SpaceIterator &iRE) {
			double logg = log10(iRE.val(DIM_E)/(p.mass*cLight2));
			while (log10(g_au[m]) < logg)
				m++;
			double slope = (safeLog10(d[m]/d[m-1]))/log10(g_au[m]/g_au[m-1]);
			double dist = slope * (logg-log10(g_au[m-1])) + safeLog10(d[m-1]);
			dist = pow(10,dist);
			p.injection.set(iRE,dist/(p.mass*cLight2));
		},{-1,iR.coord[DIM_R],0});
	},{0,-1,0});
}