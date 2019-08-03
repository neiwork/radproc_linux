#include "distribution.h"
#include "injection.h"
#include "messages.h"
#include "write.h"
#include "adafFunctions.h"
#include "globalVariables.h"
#include "losses.h"
#include <fmath/RungeKutta.h>
#include <flosses/nonThermalLosses.h>

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
	double Diffrate = diffusionRate(E,r,magfield);
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
					if (p.id == "ntProton") {
						double tcool = energy/losses(energy,p,st,itRRE);
						Nle.set(itRRE,p.injection.get(itRRE)*min(tcell,tcool)*vol);
					}
					if (p.id == "ntElectron") {
						double integ = RungeKuttaSimple(energy,p.emax()*0.99,[&](double e) {
							return p.injection.interpolate({{0,e}},&itRR.coord);});
						Nle.set(itRRE,integ/losses(energy,p,st,itRRE)*vol);
					}
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
		}

		p.ps.iterate([&](const SpaceIterator& itRR) {
			if (itRR.coord[DIM_R] != nR-1) {
				double r = itRR.val(DIM_R);
				double rB1 = r/sqrt(paso_r);
				double rB2 = rB1*paso_r;
				double vol = (4.0/3.0)*pi*costhetaH(r)*(rB2*rB2*rB2-rB1*rB1*rB1);
				p.ps.iterate([&](const SpaceIterator& itRRE) {
					p.distribution.set(itRRE,p.distribution.get(itRRE)+Nle.get(itRRE)/vol);
				},{-1,itRR.coord[DIM_R],0});
			} else {
				p.ps.iterate([&](const SpaceIterator& itRRE) {
					p.distribution.set(itRRE,p.distribution.get(itRRE)+Nle.get(itRRE));
				},{-1,itRR.coord[DIM_R],0});
			}
		},{0,-1,0});
	},{0,-1,0});
	
	if (p.id == "ntElectron") show_message(msgEnd,Module_electronDistribution);
	else if (p.id == "ntProton") show_message(msgEnd,Module_protonDistribution);
}

void distributionDetailed(Particle& p, State& st)
{
	if (p.id == "ntElectron") show_message(msgStart,Module_electronDistribution);
	else if (p.id == "ntProton") show_message(msgStart,Module_protonDistribution);

	p.ps.iterate([&](const SpaceIterator& itR) {
		
		const double r = itR.val(DIM_R);
		const double v = abs(radialVel(r));
		size_t nPoints = 50;
		const double pasoRc = pow(p.ps[DIM_R].last()/r,1.0/nPoints);
		
		p.ps.iterate([&](const SpaceIterator& itRE) {
			double e = itRE.val(DIM_E);
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
			p.distribution.set(itRE,N);
		},{-1,itR.coord[DIM_R],0});
	},{0,-1,0});
}