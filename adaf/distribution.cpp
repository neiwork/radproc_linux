#include "distribution.h"
#include "injection.h"
#include "messages.h"
#include "write.h"
#include "adafFunctions.h"
#include "globalVariables.h"
#include "losses.h"
#include <fmath/RungeKutta.h>
#include <flosses/nonThermalLosses.h>
#include <nrMath/brent.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_bessel.h>

#include <fparameters/parameters.h>
#include <iostream>

#include <fparameters/Dimension.h>
#include <fparameters/SpaceIterator.h>

int rangeE(double e, Particle& p)	{ return (e > p.emin() && e < p.emax()); }
int rangeR(double r, Particle& p)	{ return (r > p.ps[DIM_R].first() && r < p.ps[DIM_R].last()); }

double fAux(double r, double pasoR, double E, double magfield)
{
	double cos1 = costhetaH(r);
	double cos2 = costhetaH(r*pasoR);
	double cos0 = costhetaH(r*sqrt(pasoR));
	double dcos = cos2-cos1;
	double tadv = r/abs(radialVel(r));
	double Diffrate = diffusionRate(E,r*costhetaH(r),magfield);
	double q_times = tadv*Diffrate;
	double d2cos = (cos2+cos1-2.0*cos0)/(paso_r-1.0);
	return -((2.0-s-q_times)*(pasoR-1.0) + 4.0/3.0 * dcos) / (cos0 + 1.0/3.0 * dcos/ (paso_r-1.0)) - 
				d2cos/(cos1 + 1.0/3.0 * dcos);
}

double effectiveE(double Ee, double Emax, double t, Particle& p, State& st, const SpaceCoord& i)
{
	int nEeff = 200;
	double Eeff = Ee;
	double Eeff_int = pow(Emax/Ee,1.0/nEeff);
	double sum_tau = 0.0;
	while ((sum_tau < t) && (Eeff < Emax)) {
		double dEeff = Eeff/sqrt(Eeff_int)*(Eeff_int-1.0);
		double dtau_eff = dEeff/losses(Eeff,p,st,i);
		sum_tau += dtau_eff;
		Eeff *= Eeff_int;
	}
	return Eeff;
}

void distributionFast(Particle& p, State& st)
{
	if (p.id == "ntElectron") show_message(msgStart,Module_electronDistribution);
	else if (p.id == "ntProton") show_message(msgStart,Module_protonDistribution);
	
	p.ps.iterate([&](const SpaceIterator& itR) {
		
		const double r = itR.val(DIM_R);
		double delta_r = (r/sqrt(paso_r))*(paso_r-1.0);
		double tcell = delta_r/(-radialVel(r));
		
		ParamSpaceValues Nle(p.ps);
		p.ps.iterate([&](const SpaceIterator& itRR) {
			p.ps.iterate([&](const SpaceIterator& itRRE) {
				if (itRR.coord[DIM_R] == itR.coord[DIM_R]) {
					double r = itRR.val(DIM_R);
					double rB1 = r/sqrt(paso_r);
					double rB2 = rB1*paso_r;
					double vol = (4.0/3.0)*pi*costhetaH(r)*(rB2*rB2*rB2-rB1*rB1*rB1);
					double energy = itRRE.val(DIM_E);
					double tcool = energy/losses(energy,p,st,itRRE);
					if (p.id == "ntProton") {
						Nle.set(itRRE,p.injection.get(itRRE)*min(tcell,tcool)*vol);
					}
					else if (p.id == "ntElectron") {
						if (tcell < tcool) {
							Nle.set(itRRE,p.injection.get(itRRE)*tcell*vol);
						} else {
							double integ = RungeKuttaSimple(energy,p.emax()*0.99,[&](double e) {
							return p.injection.interpolate({{0,e}},&itRR.coord);});
							Nle.set(itRRE,integ/losses(energy,p,st,itRRE)*vol);
						}
					}
					/*if (p.id == "ntElectron") {
						double integ = RungeKuttaSimple(energy,p.emax()*0.99,[&](double e) {
							return p.injection.interpolate({{0,e}},&itRR.coord);});
						Nle.set(itRRE,integ/losses(energy,p,st,itRRE)*vol);
					}*/
				} else
					Nle.set(itRRE,0.0);
			},{-1,itRR.coord[DIM_R],0});
		},{0,-1,0});
		
		if (p.id == "ntProton") {
			for (int itRR=itR.coord[DIM_R]-1;itRR>=0;itRR--) {
				double rprim = p.ps[DIM_R][itRR];
				double delta_r = (rprim/sqrt(paso_r))*(paso_r-1.0);
				double vel = -radialVel(rprim);
				double tcell = delta_r/vel;
				p.ps.iterate([&](const SpaceIterator& itRRE) {  // para cada energía
					const double E = itRRE.val(DIM_E);
					double Emax = p.emax();
					double Eeff = Emax;
					if (itRRE.its[0].canPeek(+1)) { 
						Eeff = effectiveE(E,Emax,tcell,p,st,itRRE);
					}
					SpaceCoord itRRE_plus_1 = itRRE.moved({0,+1,0});
					double dist = (Eeff < p.ps[DIM_E].last()) ? 
									Nle.interpolate({{DIM_E,Eeff}},&itRRE_plus_1) : 0.0;
					double ratioLosses = losses(Eeff,p,st,itRRE_plus_1)/losses(E,p,st,itRRE);
					double dist2 = dist*ratioLosses;
					Nle.set(itRRE,dist2);
				},{-1,itRR,0});
			}
		}
		
		// REVISAR SI ESTO QUE SIGUE VA TMB PARA ELECTRONES
		if (itR.coord[DIM_R] == nR-1) {
			p.ps.iterate([&](const SpaceIterator& itRR) {
				double r = itRR.val(DIM_R);
				double rB1 = r/sqrt(paso_r);
				double rB2 = rB1*paso_r;
				double vol = (4.0/3.0)*pi*costhetaH(r)*(rB2*rB2*rB2-rB1*rB1*rB1);
				p.ps.iterate([&](const SpaceIterator& itRRE) {
					Nle.set(itRRE,Nle.get(itRRE)/vol);
				},{-1,itRR.coord[DIM_R],0});
			},{0,-1,0});
			if (p.id == "ntProton")	writeEandRParamSpace("linearEmitter_p",Nle,0);
			if (p.id == "ntElectron")	writeEandRParamSpace("linearEmitter_e",Nle,0);
		}

		p.ps.iterate([&](const SpaceIterator& itRR) {
			double r = itRR.val(DIM_R);
			double rB1 = r/sqrt(paso_r);
			double rB2 = rB1*paso_r;
			double vol = (4.0/3.0)*pi*costhetaH(r)*(rB2*rB2*rB2-rB1*rB1*rB1);
			if (itRR.coord[DIM_R] != nR-1) {
				p.ps.iterate([&](const SpaceIterator& itRRE) {
					p.distribution.set(itRRE,p.distribution.get(itRRE)+Nle.get(itRRE)/vol);
					//p.distribution.set(itRRE,p.distribution.get(itRRE)+Nle.get(itRRE));
				},{-1,itRR.coord[DIM_R],0});
			} else {
				p.ps.iterate([&](const SpaceIterator& itRRE) {
					p.distribution.set(itRRE,p.distribution.get(itRRE)+Nle.get(itRRE));
					//p.distribution.set(itRRE,p.distribution.get(itRRE)+Nle.get(itRRE)*vol);
				},{-1,itRR.coord[DIM_R],0});
			}
		},{0,-1,0});
	},{0,-1,0});
	
	///////////////////////////////////////////////////////////
	/*if (p.id == "ntElectron" || p.id == "ntPair") {
		p.ps.iterate([&](const SpaceIterator& i) {
			double tcool = i.val(DIM_E)/losses(i.val(DIM_E),p,st,i);
			p.distribution.set(i,p.injection.get(i)*tcool);
		},{-1,-1,0});
	}*/
	///////////////////////////////////////////////////////////
	if (p.id == "ntChargedPion" || p.id == "ntMuon") {
		p.ps.iterate([&](const SpaceIterator& i) {
			double r = i.val(DIM_R);
			double rB1 = r/sqrt(paso_r);
			double rB2 = rB1*paso_r;
			double tcell = (rB2-rB1)/(-radialVel(r));
			double tcool = i.val(DIM_E)/losses(i.val(DIM_E),p,st,i);
			double gamma = i.val(DIM_E)/(p.mass*cLight2);
			double tdecay;
			if (p.id == "ntChargedPion") tdecay = gamma*chargedPionMeanLife;
			if (p.id == "ntMuon") tdecay = gamma*muonMeanLife;
			p.distribution.set(i,p.injection.get(i)*min(min(tcool,tcell),tdecay));
		},{-1,-1,0});
	}
	
	////////////////////////////////////////////////////////////////////////////
	/*if (p.id == "ntElectron") {
		p.ps.iterate([&](const SpaceIterator& iR) {
			double sum = 0.0;
			double pasoE = pow(p.emax()/p.emin(),1.0/(nE-1.0));
			p.ps.iterate([&](const SpaceIterator& iRE) {
				double E = iRE.val(DIM_E);
				double dE = E*(pasoE-1.0);
				sum += p.distribution.get(iRE)*E*dE;
			},{-1,iR.coord[DIM_R],0});
			
			double norm_temp = boltzmann*st.tempElectrons.get(iR)/(electronMass*cLight2);
			double aTheta = 3.0 - 6.0/(4.0+5.0*norm_temp);
			double uth = st.denf_e.get(iR)*norm_temp*aTheta*electronMass*cLight2;
			double eta = GlobalConfig.get<double>("nonThermal.injection.energyFraction");
			double upl = eta*uth;
			
			double sum2 = 0.0;
			p.ps.iterate([&](const SpaceIterator& iRE) {
				double E = iRE.val(DIM_E);
				double dE = E*(pasoE-1.0);
				double N = p.distribution.get(iRE)*(upl/sum);
				sum2 += N*E*dE;
				p.distribution.set(iRE,N);
			},{-1,iR.coord[DIM_R],0});
			cout << "upl/uth = " << sum2/uth << "\t sum/uth = " << sum/uth << endl;
		},{0,-1,0});
	}*/
	////////////////////////////////////////////////////////////////////////////
	
	if (p.id == "ntElectron" || p.id == "ntPair") show_message(msgEnd,Module_electronDistribution);
	else if (p.id == "ntProton") show_message(msgEnd,Module_protonDistribution);
}

void distributionDetailed(Particle& p, State& st)
{
	if (p.id == "ntElectron") show_message(msgStart,Module_electronDistribution);
	else if (p.id == "ntProton") show_message(msgStart,Module_protonDistribution);

	p.ps.iterate([&](const SpaceIterator& itR) {
		
		const double r = itR.val(DIM_R);
		const double v = abs(radialVel(r));
		double rB1 = r/sqrt(paso_r);
		double rB2 = rB1*paso_r;
		double vol = 4.0/3.0 * pi * cos(st.thetaH.get(itR)) * (rB2*rB2*rB2-rB1*rB1*rB1);
		size_t nPoints = 30;
		const double pasoRc = pow(p.ps[DIM_R].last()/r,1.0/nPoints);
		double sumNT = 0.0;
		p.ps.iterate([&](const SpaceIterator& itRE) {
			double e = itRE.val(DIM_E); // REVISAR SI ESTA E SIRVE PARA LAS PARTICULAS
			Vector Rc(nPoints,r);
			Vector Ec(nPoints,e);
			
			// CHARACTERISTIC CURVE ///////////////////////////
			for (size_t j=0;j<nPoints-1;j++) {
				double vr = 1.0;
				double be = 0.0;
				double deriv = 0.0;
				if (rangeR(Rc[j],p)) {
					vr = radialVel(Rc[j]);
					be = b(Ec[j],Rc[j],p,st,itRE);
					deriv = be/vr;
				}
				double dRc = Rc[j]*(pasoRc-1.0);
				Ec[j+1] = Ec[j] + deriv*dRc;
				Rc[j+1] = Rc[j]*pasoRc;
			}
			///////////////////////////////////////////////////
			
			double N = 0.0;
			for (size_t j=0;j<nPoints;j++) {
				double Qinj = (rangeE(Ec[j],p) && rangeR(Rc[j],p)) ? 
						p.injection.interpolate({{DIM_E,Ec[j]},{DIM_R,Rc[j]}},&itRE.coord) : 0.0;
				double mu = 0.0;
				double pasoEaux = 1.00001;
				for (size_t jj=0;jj<j;jj++) {
					double dEe = Ec[jj]*(sqrt(pasoEaux)-1.0/sqrt(pasoEaux));
					double dbde = ( (rangeE(Ec[jj]*sqrt(pasoEaux),p) &&
							rangeE(Ec[jj]/sqrt(pasoEaux),p)) && rangeR(Rc[jj],p) ) ?
							(b(Ec[jj]*sqrt(pasoEaux),Rc[j],p,st,itRE)-
									b(Ec[jj]/sqrt(pasoEaux),Rc[j],p,st,itRE))/dEe : 0.0;
					double vR = abs(radialVel(Rc[jj]));
					double dRR = Rc[jj]*(pasoRc-1.0);
					double magfield = rangeR(Rc[jj],p) ? st.magf.interpolate({{DIM_R,Rc[jj]}},&itRE.coord) : 0.0;
					double f = fAux(Rc[jj],pasoRc,Ec[jj],magfield);
					mu += f + dbde*(dRR/vR);
				}
				mu = exp(-mu);
				double dRc = Rc[j]*(pasoRc-1.0);
				N += Qinj * mu * dRc;
			}
			N /= v;
			double pasoE = pow(p.emax()/p.emin(),1.0/nE);
			p.distribution.set(itRE,N);
			double dE = e*(pasoE-1.0);
			sumNT += N*e*dE;
		},{-1,itR.coord[DIM_R],0});
		double dens = st.denf_i.get(itR);
		double norm_temp = boltzmann * st.tempIons.get(itR) / (p.mass*cLight2);
		double aTheta = 3.0 - 6.0/(4.0+5.0*norm_temp);
		double uth = dens * norm_temp * p.mass * cLight2 * aTheta;
		cout << "uThermal = " << uth << "\t uNonThermal = " << sumNT << endl;
	},{0,-1,0});
}

double gammaMin(double g, void *params)
{
  struct one_d_params *pa = (struct one_d_params *) params;
  double norm_temp = pa->p;
  double pIndex = GlobalConfig.get<double>("nonThermal.injection.primaryIndex");
  double eta = GlobalConfig.get<double>("nonThermal.injection.energyFraction");
  double beta = sqrt(1.0-1.0/(g*g));
  double aTheta = 3.0 - 6.0/(4.0+5.0*norm_temp); // Gammie & Popham (1998)
  double bessel = gsl_sf_bessel_Kn(2,1.0/norm_temp);
  return pow(g,4.0)*beta*exp(-g/norm_temp)-(pIndex-2.0)*eta*aTheta*norm_temp*norm_temp*bessel;
}

void distributionSimplified(Particle& p, State& st)
{
	if (p.id == "ntElectron") show_message(msgStart,Module_electronDistribution);

	p.ps.iterate([&](const SpaceIterator& itR) {
		
		const double r = itR.val(DIM_R);
		const double v = abs(radialVel(r));
		double tacc = r/v;
		double magf = st.magf.get(itR);
		double gamma_cool = 6.0*pi*electronMass*cLight/ (thomson*tacc*magf*magf);
		double norm_temp = boltzmann*st.tempElectrons.get(itR)/(electronMass*cLight2);
		double gamma_min = 1.0;
		double pIndex = GlobalConfig.get<double>("nonThermal.injection.primaryIndex");
		double eta = GlobalConfig.get<double>("nonThermal.injection.energyFraction");
		double accE = GlobalConfig.get<double>("nonThermal.injection.accEfficiency");
		double gamma_max = sqrt(accE*6.0*pi*electronCharge/(thomson*magf));
		
		double nThermal = st.denf_e.get(itR);
		double aTheta = 3.0 - 6.0/(4.0+5.0*norm_temp); // Gammie & Popham (1998)
		double uThermal = aTheta*nThermal*electronMass*cLight2*norm_temp;
		double uPl = eta*uThermal;
		double int_E = 0.0;
		int_E = (gamma_cool < gamma_max) ? 
			(pow(gamma_min,2.0-pIndex)-pow(gamma_cool,2.0-pIndex))/(pIndex-2.0) +
				(pow(gamma_cool,2.0-pIndex)-gamma_cool*pow(gamma_max,1.0-pIndex))/(pIndex-1.0) :
			(pow(gamma_min,2.0-pIndex)-pow(gamma_max,2.0-pIndex))/(pIndex-2.0);
		double nPl = uPl / int_E / electronMass / cLight2;
		
		double pasoE = pow(p.emax()/p.emin(),1.0/(nE-1.0));
		double sumNT = 0.0;
		p.ps.iterate([&](const SpaceIterator& itRE) {
			double e = itRE.val(DIM_E);
			double gamma = e/(electronMass*cLight2);
			double N = 0.0;
			if (gamma > gamma_min) {
				if (gamma < gamma_cool) {
					N = exp(-gamma/gamma_max) * nPl * pow(gamma,-pIndex);
				} else {
					N = exp(-gamma/gamma_max) * nPl * gamma_cool * pow(gamma,-pIndex-1.0);
				}
				N /= (electronMass*cLight2);
			}
			p.distribution.set(itRE,N);
			double dE = e*(pasoE-1.0);
			sumNT += N*e*dE;
		},{-1,itR.coord[DIM_R],0});;
		cout << "uThermal = " << uThermal << "\t uNonThermal = " << sumNT << endl;
	},{0,-1,0});
}