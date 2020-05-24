#include "NTdistribution.h"
#include "NTinjection.h"
#include "messages.h"
#include "write.h"
#include "adafFunctions.h"
#include <fmath/fbisection.h>
#include "globalVariables.h"
#include "losses.h"
#include <fmath/RungeKutta.h>
#include <flosses/nonThermalLosses.h>
#include <flosses/lossesSyn.h>
#include <nrMath/brent.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_randist.h>
#include <fmath/tridiagonal.h>

#include <fparameters/parameters.h>
#include <iostream>

#include <fparameters/Dimension.h>
#include <fparameters/SpaceIterator.h>

/*
int rangeE(double e, Particle& p)	{ return (e > p.emin() && e < p.emax()); }
int rangeR(double r, Particle& p)	{ return (r > p.ps[DIM_R].first() && r < p.ps[DIM_R].last()); }

double fAux(double r, double pasoR, double E, double magfield, Particle& p)
{
	double cos1 = costhetaH(r);
	double cos2 = costhetaH(r*pasoR);
	double cos0 = costhetaH(r*sqrt(pasoR));
	double dcos = cos2-cos1;
	double tadv = r/abs(radialVel(r));
	double height = height_fun(r);
	double Diffrate = 1.0/diffusionTimeTurbulence(E,height,p,magfield);
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
			if (itRR.coord[DIM_R] == itR.coord[DIM_R]) {
				double r = itRR.val(DIM_R);
				double rB1 = r/sqrt(paso_r);
				double rB2 = rB1*paso_r;
				double vol = volume(r);
				p.ps.iterate([&](const SpaceIterator& itRRE) {
					double energy = itRRE.val(DIM_E);
					double tcool = energy/losses(energy,p,st,itRRE);
					if (p.id == "ntProton") {
						Nle.set(itRRE,p.injection.get(itRRE)*min(tcell,tcool)*vol);
					}
					else if (p.id == "ntElectron" || p.id == "ntPair") {
						if (tcell < tcool) {
							Nle.set(itRRE,p.injection.get(itRRE)*tcell*vol);
						} else {
							//double integ = RungeKuttaSimple(energy,p.emax()*0.99,[&](double e) {
							//return p.injection.interpolate({{0,e}},&itRR.coord);});
							//Nle.set(itRRE,integ/losses(energy,p,st,itRRE)*vol);
							double integ = integSimpson(log(energy),log(p.emax()),[&](double loge)
								{
									double e = exp(loge);
									return p.injection.interpolate({{0,e}},&itRRE.coord)*e;
								},100);
							Nle.set(itRRE,integ/losses(energy,p,st,itRRE)*vol);
						}
					}
				},{-1,itRR.coord[DIM_R],0});
			} else
				p.ps.iterate([&](const SpaceIterator& itRRE) {
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
				double vol = volume(itRR.val(DIM_R));
				p.ps.iterate([&](const SpaceIterator& itRRE) {
					Nle.set(itRRE,Nle.get(itRRE)/vol);
				},{-1,itRR.coord[DIM_R],0});
			},{0,-1,0});
			if (p.id == "ntElectron") writeEandRParamSpace("linearEmitter_e",Nle,0,1);
			if (p.id == "ntProton")	writeEandRParamSpace("linearEmitter_p",Nle,0,1);
		}

		p.ps.iterate([&](const SpaceIterator& itRR) {
			double vol = volume(itRR.val(DIM_R));
			if (itRR.coord[DIM_R] != nR-1) {
				p.ps.iterate([&](const SpaceIterator& itRRE) {
					p.distribution.set(itRRE,p.distribution.get(itRRE)+Nle.get(itRRE)/vol);
				},{-1,itRR.coord[DIM_R],0});
			} else {
				p.ps.iterate([&](const SpaceIterator& itRRE) {
					p.distribution.set(itRRE,p.distribution.get(itRRE)+Nle.get(itRRE));
				},{-1,itRR.coord[DIM_R],0});
			}
		},{0,-1,0});
		double dens = st.denf_i.get(itR);
		double norm_temp = boltzmann * st.tempIons.get(itR) / (p.mass*cLight2);
		double aTheta = 3.0 - 6.0/(4.0+5.0*norm_temp);
		double uth = dens * norm_temp * p.mass * cLight2 * aTheta;
		double unth = integSimpson(log(p.emin()),log(p.emax()),[&](double loge)
						{
							double e = exp(loge);
							return e*e*p.distribution.interpolate({{DIM_E,e}},&itR.coord);
						},100);
		cout << "uThermal = " << uth << "\t uNonThermal = " << scientific << unth/uth << endl;
		
	},{0,-1,0});
	
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
					double f = fAux(Rc[jj],pasoRc,Ec[jj],magfield,p);
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
}*/

/*
double escapeRateNew(double E, double r, Particle& p)
{
	double B = magneticField(r);
	double h = height_fun(r);
	return 1.0/diffusionTimeTurbulence(E,h,p,B);
}

double logEaux(double logRaux, Vector logRc, Vector logEc)
{
	double aux = logRc[0];
	size_t pos = 0;
	while(aux < logRaux) {
		pos++;
		aux = logRc[pos];
	}
	double m = (logEc[pos]-logEc[pos-1])/(logRc[pos]-logRc[pos-1]);
	return m*(logRaux-logRc[pos-1]) + logEc[pos-1];
}

double f1(double logr, double loge, Particle& p)
{
	double r = exp(logr);
	double e = exp(loge);
	double tacc = r/(-radialVel(r));
	double tcool = t_cool(e,r,p);
	return tacc/tcool;
}

double f2(double logr, Vector logRc, Vector logEc, Particle& p)
{
	double r = exp(logr);
	double E = exp(logEaux(logr,logRc,logEc));
	double tacc = r/(-radialVel(r));
	double tesc = 1.0/escapeRateNew(E,r,p);
	return tacc/tesc;
}

void distributionNew(Particle& p, State& st)
{
	if (p.id == "ntElectron") show_message(msgStart,Module_electronDistribution);
	else if (p.id == "ntProton") show_message(msgStart,Module_protonDistribution);

	p.ps.iterate([&](const SpaceIterator& itR) {
		
		const double logr = log(itR.val(DIM_R));
		double r = exp(logr);
		double h = height_fun(r);
		size_t nPoints = 50;
		const double pasoRc = pow(p.ps[DIM_R].last()/r,1.0/nPoints);
		double dlogRc = pasoRc-1.0;
		p.ps.iterate([&](const SpaceIterator& itRE) {
			double e = itRE.val(DIM_E);
			double loge = log(e);
			Vector logRc(nPoints,logr);
			Vector logEc(nPoints,loge);
			
			// CHARACTERISTIC CURVE ///////////////////////////
			for (size_t j=0;j<nPoints;j++) {
				double d_logEc = 0.0;
				if (rangeR(exp(logRc[j]),p))
					d_logEc = RungeKuttaStep([&p](double logx, double logy)
								{return f1(logx,logy,p);},logRc[j],logEc[j],dlogRc);
				logEc[j+1] = logEc[j] + d_logEc;
				logRc[j+1] = logRc[j] + dlogRc;
			}
			///////////////////////////////////////////////////
			
			double dist = 0.0;
			for (size_t j=0;j<nPoints;j++) {
				double Ec = exp(logEc[j]);
				Ec = exp(loge);
				double Rc = exp(logRc[j]);
				if (rangeE(Ec,p) && rangeR(Rc,p)) {
					double Qinj = p.injection.interpolate({{DIM_E,Ec},{DIM_R,Rc}},&itRE.coord);
					double mu = 0.0;
					size_t nAux = 30;
					mu = integSimpson(logr,logRc[j],[&logRc,&logEc,&p](double logx)
								{return f2(logx,logRc,logEc,p);},nAux);
					mu = exp(-mu);
					double Hc = height_fun(Rc);
					double tacc = Rc/(-radialVel(Rc));
					dist += Qinj * tacc * pow(Rc/r,1.0-s) * (Hc/h) * mu * dlogRc;
				}
			}
			p.distribution.set(itRE,dist);
		},{-1,itR.coord[DIM_R],0});
		double dens = st.denf_i.get(itR);
		double norm_temp = boltzmann * st.tempIons.get(itR) / (p.mass*cLight2);
		double aTheta = 3.0 - 6.0/(4.0+5.0*norm_temp);
		double uth = dens * norm_temp * p.mass * cLight2 * aTheta;
		double unth = integSimpson(log(p.emin()),log(p.emax()),[&](double loge)
						{
							double e = exp(loge);
							return e*e*p.distribution.interpolate({{DIM_E,e}},&itR.coord);
						},100);
		cout << "uNTh/uTh = " << unth/uth << endl;
	},{0,-1,0});
}

double logRaux(double logEaux, Vector logEc, Vector logRc)
{
	double aux = logEc[0];
	size_t pos = 0;
	while(aux < logEaux) {
		pos++;
		aux = logEc[pos];
	}
	double m = (logRc[pos]-logRc[pos-1])/(logEc[pos]-logEc[pos-1]);
	return m*(logEaux-logEc[pos-1]) + logRc[pos-1];
}


double f1_e(double loge, double logr, Particle& p)
{
	double r = exp(logr);
	double e = exp(loge);
	double tacc = r/(-radialVel(r));
	double tcool = t_cool(e,r,p);
	return tcool/tacc;
}

double f2_e(double loge, Vector logEc, Vector logRc, Particle& p)
{
	double e = exp(loge);
	double R = exp(logRaux(loge,logEc,logRc));
	double tcool = t_cool(e,R,p);
	double tesc = 1.0/escapeRateNew(e,R,p);
	return tcool/tesc;
}

void distributionNewE(Particle& p, State& st)
{
	if (p.id == "ntElectron") show_message(msgStart,Module_electronDistribution);
	else if (p.id == "ntProton") show_message(msgStart,Module_protonDistribution);

	p.ps.iterate([&](const SpaceIterator& itR) {
		
		const double logr = log(itR.val(DIM_R));
		double r = exp(logr);
		double h = height_fun(r);
		size_t nPoints = 30;
		p.ps.iterate([&](const SpaceIterator& itRE) {
			double e = itRE.val(DIM_E);
			const double pasoEc = pow(p.emax()/e,1.0/(nPoints-1));
			double dlogEc = pasoEc-1.0;
			double loge = log(e);
			Vector logRc(nPoints,logr);
			Vector logEc(nPoints,loge);
			
			// CHARACTERISTIC CURVE ///////////////////////////
			for (size_t j=0;j<nPoints;j++) {
				double d_logRc = 0.0;
				if (rangeE(exp(logEc[j]),p))
					d_logRc = RungeKuttaStep([&p](double logx, double logy)
								{return f1_e(logx,logy,p);},logEc[j],logRc[j],dlogEc);
				logEc[j+1] = logEc[j] + dlogEc;
				logRc[j+1] = logRc[j] + d_logRc;
			}
			///////////////////////////////////////////////////
			
			double dist = 0.0;
			for (size_t j=0;j<nPoints;j++) {
				double Ec = exp(logEc[j]);
				double Rc = exp(logRc[j]);
				if (rangeE(Ec,p) && rangeR(Rc,p)) {
					double Qinj = p.injection.interpolate({{DIM_E,Ec},{DIM_R,Rc}},&itRE.coord);
					double mu = 0.0;
					size_t nAux = 30;
					mu = integSimpson(loge,logEc[j],[&logEc,&logRc,&p](double logx)
								{return f2_e(logx,logEc,logRc,p);},nAux);
					mu = exp(-mu);
					double Hc = height_fun(Rc);
					double tcool = t_cool(Ec,Rc,p);
					dist += Qinj * tcool * pow(Rc/r,1.0-s) * (Hc/h) * mu * dlogEc;
				}
			}
			p.distribution.set(itRE,dist);
		},{-1,itR.coord[DIM_R],0});
		double dens = st.denf_i.get(itR);
		double norm_temp = boltzmann * st.tempIons.get(itR) / (p.mass*cLight2);
		double aTheta = 3.0 - 6.0/(4.0+5.0*norm_temp);
		double uth = dens * norm_temp * p.mass * cLight2 * aTheta;
		double unth = integSimpson(log(p.emin()),log(p.emax()),[&](double loge)
						{
							double e = exp(loge);
							return e*e*p.distribution.interpolate({{DIM_E,e}},&itR.coord);
						},100);
		cout << "uNTh/uTh = " << unth/uth << endl;
	},{0,-1,0});
}*/

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

/*
void distributionOneZone_analytical(Particle& p, State& st)
{
	if (p.id == "ntElectron") show_message(msgStart,Module_electronDistribution);
	if (p.id == "ntChargedPion") show_message(msgStart,Module_pionDistribution);
	if (p.id == "ntMuon") show_message(msgStart,Module_muonDistribution);
	if (p.id == "ntPair") show_message(msgStart,Module_pairDistribution);
	p.ps.iterate([&](const SpaceIterator& iR) {
		double tacc = accretionTime(iR.val(DIM_R));
		double magf = st.magf.get(iR);
		double height = height_fun(iR.val(DIM_R));
		p.ps.iterate([&](const SpaceIterator& iRE) {
			double E = iRE.val(DIM_E);
			double dEdt = abs(losses(E,p,st,iRE));
			double tCool = E/dEdt;
			double tEscape = pow(1.0/diffusionTimeTurbulence(E,height,p,magf)+1.0/tacc,-1);
			double tDecayE = 1.0e99;
			if (p.id == "ntChargedPion") tDecayE = chargedPionMeanLife;
			if (p.id == "ntMuon") tDecayE = muonMeanLife;
			tDecayE *= (E/(p.mass*cLight2));
			if (tDecayE < 0.1*tCool && tDecayE < 0.1*tEscape) {
				p.distribution.set(iRE,p.injection.get(iRE)*tDecayE);
			} else if (tEscape < 0.1*tDecayE && tEscape < 0.1*tCool) {
				p.distribution.set(iRE,p.injection.get(iRE)*tEscape);
			} else {
				double integ = qImpropLog(E,p.emax(),[E,magf,height,tacc,&st,&p,&iRE](double Ep)
						{
							double mu = qImpropLog(E,Ep,[magf,height,tacc,&st,&p,&iRE](double Epp)
							{
								double tdiff = diffusionTimeTurbulence(Epp,height,p,magf);
								double rateEscape = 1.0/tacc + 1.0/tdiff;
								double tDecay = 1.0e99;
								if (p.id == "ntChargedPion") tDecay = chargedPionMeanLife;
								if (p.id == "ntMuon") tDecay = muonMeanLife;
								tDecay *= (Epp/(p.mass*cLight2));
								double rateDecay = 1.0/tDecay;
								double tcool = Epp/abs(losses(Epp,p,st,iRE));
								return tcool * (rateEscape+rateDecay) / Epp;
							},1.0e-3);
							double Q = p.injection.interpolate({{DIM_E,Ep}},&iRE.coord);
							return Q*exp(-mu);
						},1.0e-3);
				p.distribution.set(iRE,integ/dEdt);
			}
		},{-1,iR.coord[DIM_R],0});
	},{0,-1,0});
	if (p.id == "ntElectron") show_message(msgEnd,Module_electronDistribution);
	if (p.id == "ntChargedPion") show_message(msgEnd,Module_pionDistribution);
	if (p.id == "ntMuon") show_message(msgEnd,Module_muonDistribution);
	if (p.id == "ntPair") show_message(msgEnd,Module_pairDistribution);
}



void distributionMultiZone_analytical(Particle& p, State& st)
{
	if (p.id == "ntProton") show_message(msgStart,Module_protonDistribution);
	if (p.id == "ntElectron") show_message(msgStart,Module_electronDistribution);
	p.ps.iterate([&](const SpaceIterator& iR) {
		int jR = nR-1-iR.coord[DIM_R];
		SpaceCoord jRsc = {0,jR,0};
		double r = p.ps[DIM_R][jR];
		double vol = volume(r);
		double rB2 = r*sqrt(paso_r);
		double rB1 = rB2/paso_r;
		double rho = massDensityADAF(r);
		double rateAccretion = abs(radialVel(rB1))*4.0*pi*height_fun(rB1)*rB1/vol;
		double rateWind = (accRateADAF(rB2)-accRateADAF(rB1)) / (rho*vol);
		double magf = st.magf.get(jRsc);
		double height = height_fun(r);
		p.ps.iterate([&](const SpaceIterator& jRE) {
			double E = jRE.val(DIM_E);
			double dEdt = losses(E,p,st,jRE);
			double rateCool = dEdt/E;
			double rateDiff = 1.0/diffusionTimeTurbulence(E,height,p,magf);
			rateDiff = 0.0;
			double rateEscape = rateDiff + rateWind + rateAccretion;
			if (rateEscape/rateCool > 10) {
				double Q1(0.0), Q2(0.0);
				Q1 = p.injection.get(jRE);
				double rateEff = rateEscape + rateCool;
				if (jR < nR-1) {
					SpaceCoord jRplus = {jRE.coord[DIM_E],jR+1,0};
					double areaPlus = 4.0*pi*height_fun(rB2)*rB2;
					Q2 = p.distribution.get(jRplus) * areaPlus * (-radialVel(rB2)) / vol;
					p.distribution.set(jRE, (Q1+Q2) / rateEff);
				} else
					p.distribution.set(jRE, Q1 / rateEff);
			} else {
				double integ = integSimpsonLog(E,p.emax(),
						[E,magf,height,vol,rB2,jR,rateAccretion,rateWind,&st,&p,&jRE](double Ep)
						{
							double mu = integSimpsonLog(E,Ep,
							[magf,height,rateAccretion,rateWind,&st,&p,&jRE] (double Epp)
							{
								double rateDiff_pp = 1.0/diffusionTimeTurbulence(Epp,height,p,magf);
								rateDiff_pp = 0.0;
								double rateEscape_pp = rateDiff_pp + rateWind + rateAccretion;
								double rateCool_pp = losses(Epp,p,st,jRE)/Epp;
								return rateEscape_pp / Epp / rateCool_pp;
							},100);
							double Q1 = p.injection.interpolate({{DIM_E,Ep}},&jRE.coord);
							double Q2 = 0.0;
							if (jR < nR-1) {
								SpaceCoord jRplus = {jRE.coord[DIM_E],jR+1,0};
								double areaPlus = 4.0*pi*height_fun(rB2)*rB2;
								Q2 = p.distribution.interpolate({{DIM_E,Ep}},&jRplus)*areaPlus*(-radialVel(rB2)) / vol;
							}
							return (Q1+Q2)*exp(-mu);
						},30);
				p.distribution.set(jRE,integ/dEdt);
			}
		},{-1,jR,0});
		
		double dens = st.denf_i.get(jRsc);
		double norm_temp = boltzmann * st.tempIons.get(jRsc) / (p.mass*cLight2);
		double aTheta = 3.0 - 6.0/(4.0+5.0*norm_temp);
		double uth = dens * norm_temp * p.mass * cLight2 * aTheta;
		double unth = integSimpsonLog(p.emin(),p.emax(),[&](double e)
						{
							return e*p.distribution.interpolate({{DIM_E,e}},&jRsc);
						},100);
		cout << "uNTh/uTh = " << unth/uth << endl;
	},{0,-1,0});
	if (p.id == "ntProton") show_message(msgEnd,Module_protonDistribution);
}

void distributionFokkerPlanck_momentum(Particle& particle, State& st)
{
	particle.ps.iterate([&](const SpaceIterator& iR) {
		
		double r = iR.val(DIM_R);
		double vol = volume(r);
		double height = height_fun(r);
		double B = st.magf.get(iR);
		double rho = massDensityADAF(r);
		double rateAccretion = 1.0/accretionTime(r);
		double rateWind = s*rateAccretion;
		
		cout << "vA = " << B/sqrt(4.0*pi*rho)/cLight << endl;
		fun1 g_eq_fun = [rateAccretion,height,B,rho,&particle](double g)
				{
					double tAcceleration = g*g/diffCoeff_g(g,particle,height,B,rho);
					return tAcceleration - 1.0/rateAccretion;
				};
		double g_eq;
		int bis = Bisect(g_eq_fun,1.0e4,1.0e9,g_eq);
		cout << "Gamma equilibrium = " << safeLog10(g_eq) << endl;
		
		// We define a new mesh of points
		size_t M = 99;
		double p_min = particle.emin()/(particle.mass*cLight2);
		double p_max = particle.emax()/(particle.mass*cLight2);		// In units of mc  --> p = gamma*beta
		
		double paso_p = pow(p_max/p_min,1.0/M);
		Vector p_au(M+3,p_min/paso_p);
		Vector delta_p_m_au(M+2,0.0);
		Vector delta_p(M+1,0.0);
		
		double p_inj = 2.0*p_min;
		double sigma = p_inj * (paso_p-1.0)/10;
		double F0 = 1.0;
		
		for (size_t jp=1;jp<M+3;jp++)	p_au[jp] = p_au[jp-1]*paso_p;
		for (size_t jp=0;jp<M+2;jp++)	delta_p_m_au[jp] = p_au[jp+1]-p_au[jp];
		for (size_t jp=0;jp<M+1;jp++)	delta_p[jp] = 0.5*(p_au[jp+2]-p_au[jp]);
		
		fun1 Afun = [](double p) { return p*p; };
		
		fun1 Bfun = [&particle,&st,&iR] (double p) 
						{
							double E = sqrt(p*p+1.0)*particle.mass*cLight2;
							double tCool = E/abs(losses(E,particle,st,iR));
							return p*p*p/tCool;
						};

		fun1 Cfun = [&particle,B,height,rho] (double p)
						{
							double E = sqrt(p*p+1.0)*particle.mass*cLight2;
							double Dp = diffCoeff_p(E,particle,height,B,rho)/gsl_pow_2(particle.mass*cLight);
							return p*p*Dp;
						};

		fun1 Tfun = [r,&particle,height,B,rateAccretion,rateWind] (double p) 
						{
							double E = sqrt(p*p+1.0)*particle.mass*cLight2;
							return pow(rateAccretion+rateWind+1.0/diffusionTimeTurbulence(E,height,particle,B),-1);
						};
		
		fun1 Qfun = [F0,p_inj,sigma] (double p)
						{ 
							return F0*gsl_ran_gaussian_pdf(p-p_inj,sigma);
						};

		Vector d(M+1,0.0);
		ofstream file;
		file.open("fokkerPlanck.dat");
		for (size_t m=0;m<=M;m++)
			file << safeLog10(p_au[m]) << "\t" << 4.0*pi*gsl_pow_2(p_au[m])*Qfun(p_au[m]) << endl;
		file.close();

		Vector a(M+1,0.0), b(M+1,0.0), c(M+1,0.0);
		for (size_t m=0;m<M+1;m++) {
			double Am = p_au[m+1]*p_au[m+1];
			double Bm_minusHalf = 0.5 * (Bfun(p_au[m]) + Bfun(p_au[m+1]));
			double Bm_plusHalf = 0.5 * (Bfun(p_au[m+1]) + Bfun(p_au[m+2]));
			double Cm_minusHalf = 0.5 * (Cfun(p_au[m]) + Cfun(p_au[m+1]));
			double Cm_plusHalf = 0.5 * (Cfun(p_au[m+1]) + Cfun(p_au[m+2]));
			double Qm = Qfun(p_au[m+1]);
			double Tm = Tfun(p_au[m+1]);
			double wm_minusHalf = Bm_minusHalf/Cm_minusHalf * delta_p_m_au[m];
			double wm_plusHalf = Bm_plusHalf/Cm_plusHalf * delta_p_m_au[m+1];
			double Wm_minusHalf = (abs(wm_minusHalf) < 1.0e-3) ? 
					pow(1.0+gsl_pow_2(wm_minusHalf)/24+gsl_pow_4(wm_minusHalf)/1920,-1) :
					abs(wm_minusHalf)*exp(-0.5*abs(wm_minusHalf))/(1.0-exp(-abs(wm_minusHalf)));
			double Wm_plusHalf = (abs(wm_plusHalf) < 1.0e-3) ? 
					pow(1.0+gsl_pow_2(wm_plusHalf)/24+gsl_pow_4(wm_plusHalf)/1920,-1) :
					abs(wm_plusHalf)*exp(-0.5*abs(wm_plusHalf))/(1.0-exp(-abs(wm_plusHalf)));
			double Wminus_m_minusHalf = Wm_minusHalf*exp(-0.5*wm_minusHalf);
			double Wminus_m_plusHalf = Wm_plusHalf*exp(-0.5*wm_plusHalf);
			double Wplus_m_minusHalf = Wm_minusHalf*exp(0.5*wm_minusHalf);
			double Wplus_m_plusHalf = Wm_plusHalf*exp(0.5*wm_plusHalf);
			
			if (m >= 1)
				a[m] = - Cm_minusHalf * Wminus_m_minusHalf / 
							( Am * delta_p[m] * delta_p_m_au[m] );
			if (m <= M-1)
				c[m] = - Cm_plusHalf * Wplus_m_plusHalf / 
							( Am * delta_p[m] * delta_p_m_au[m+1] );
			b[m] = 1.0 / (delta_p[m] * Am) * 
					( Cm_minusHalf * Wplus_m_minusHalf / delta_p_m_au[m] +
						Cm_plusHalf * Wminus_m_plusHalf / delta_p_m_au[m+1] ) + 1.0/Tm;
			d[m] = Qm;
		}
		TriDiagSys(a,b,c,d,M);
		size_t m = 0;
		particle.ps.iterate([&](const SpaceIterator &iRE) {
			particle.distribution.set(iRE,4.0*pi*gsl_pow_2(p_au[m+1])*d[m]/(particle.mass*cLight2));
			m++;
		},{-1,iR.coord[DIM_R],0});
		
		file.open("fokkerPlanck.dat");
		for (size_t m=0;m<=M;m++) {
			double E = particle.mass*cLight2 * p_au[m+1];
			double E_GeV = E / 1.602e-3;
			file << safeLog10(E_GeV) << "\t" << 4.0*pi*gsl_pow_2(p_au[m+1])*d[m]/cLight*vol *E*E << endl;
		}
		file.close();
		
		double aux=1.0;
		double aux1 = aux;
	},{0,-1,0});
	
	// NORMALIZATION
	double sum = 0.0;
	particle.ps.iterate([&](const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		double height = height_fun(r);
		double B = st.magf.get(iR);
		double vol = volume(r);
		double tAccretion = accretionTime(r);
		double integ = integSimpsonLog(particle.emin(),particle.emax(),[&particle,height,B,tAccretion,&iR](double E)
						{
							double tDiff = diffusionTimeTurbulence(E,height,particle,B);
							double rateEscape = 1.0/tAccretion + 1.0/tDiff;
							return particle.distribution.interpolate({{DIM_E,E}},&iR.coord)*E*rateEscape;
						},50);
		sum += vol*integ;
	},{0,-1,0});
	double factor = etaInj*accRateOut*cLight2 / sum;
	particle.ps.iterate([&](const SpaceIterator& i) {
		particle.distribution.set(i,particle.distribution.get(i)*factor);
	},{-1,-1,0});
}*/

void distributionMultiZone(Particle& particle, State& st)
{
	particle.ps.iterate([&](const SpaceIterator& iR) {
		
		SpaceCoord jR = {0,nR-1-iR.coord[DIM_R],0};
		double r = particle.ps[DIM_R][jR[DIM_R]];
		double rB2 = r*sqrt(paso_r);
		double rB1 = rB2/paso_r;
		double vR = radialVel(r);
		double area = 4.0*pi*height_fun(r)*r;
		double vol = volume(r);
		double height = height_fun(r);
		double B = st.magf.get(jR);
		double rho = massDensityADAF(r);
		double rateAccretion = abs(radialVel(rB1))*4.0*pi*height_fun(rB1)*rB1/vol;
		double rateWind = (accRateADAF(rB2)-accRateADAF(rB1)) / (rho*vol);
		
		// We define a new mesh of points
		size_t M = 199;
		double g_min = particle.emin()/(particle.mass*cLight2);
		double g_max = particle.emax()/(particle.mass*cLight2);
		
		double paso_g = pow(g_max/g_min,1.0/M);
		Vector g_au(M+3,g_min/paso_g);
		Vector delta_g_m_au(M+2,0.0);
		Vector delta_g(M+1,0.0);
		
		for (size_t jg=1;jg<M+3;jg++)	g_au[jg] = g_au[jg-1]*paso_g;
		for (size_t jg=0;jg<M+2;jg++)	delta_g_m_au[jg] = g_au[jg+1]-g_au[jg];
		for (size_t jg=0;jg<M+1;jg++)	delta_g[jg] = 0.5*(g_au[jg+2]-g_au[jg]);
		
		fun1 Bfun = [&particle,&st,&jR,B,rho,height] (double g) 
						{
							double E = g*particle.mass*cLight2;
							double tCool = E/losses(E,particle,st,jR);
							return g/tCool;
						};

		fun1 Tfun = [r,&particle,height,B,rateAccretion,&rateWind] (double g) 
						{
							double E = g*particle.mass*cLight2;
							double rateDiff = 1.0/diffusionTimeTurbulence(E,height,particle,B);
							double rateDecay = 0.0;
							if (particle.id == "ntChargedPion")
								rateDecay = 1.0 / (g*chargedPionMeanLife);
							else if (particle.id == "ntMuon")
								rateDecay = 1.0 / (g*muonMeanLife);
	
							return pow(rateAccretion+rateWind+rateDiff,-1);
						};
		
		fun1 Qfun = [&particle,&jR,rB2,vol] (double g)
						{
							double E = g*particle.mass*cLight2;
							double Qplus = 0.0;
							if (jR[DIM_R] < nR-1) {
								SpaceCoord jRplus = {0,jR[DIM_R]+1,0};
								double areaPlus = 4.0*pi*height_fun(rB2)*rB2;
								double Nplus = particle.distribution.interpolate({{DIM_E,E}},&jRplus)*(particle.mass*cLight2);
								Qplus = Nplus * areaPlus * abs(radialVel(rB2)) / vol;
							}
							return particle.injection.interpolate({{DIM_E,E}},&jR)*(particle.mass*cLight2) + Qplus;
						};

		Vector d(M+1,0.0);
		Vector a(M+1,0.0), b(M+1,0.0), c(M+1,0.0);
		for (size_t m=0;m<M+1;m++) {
			double Bm_minusHalf = 0.5 * (Bfun(g_au[m]) + Bfun(g_au[m+1]));
			double Bm_plusHalf = 0.5 * (Bfun(g_au[m+1]) + Bfun(g_au[m+2]));
			double Qm = Qfun(g_au[m+1]);
			double Tm = Tfun(g_au[m+1]);

			if (m <= M-1)
				c[m] = - abs(Bm_plusHalf) / delta_g[m];
			b[m] = 1.0/delta_g[m] * abs(Bm_minusHalf) + 1.0/Tm;
			d[m] = Qm;
		}
		
		TriDiagSys(a,b,c,d,M);
		size_t m = 0;
		particle.ps.iterate([&](const SpaceIterator &jRE) {
			double logg = log10(jRE.val(DIM_E)/(particle.mass*cLight2));
			while (log10(g_au[m]) < logg)
				m++;
			double slope = (safeLog10(d[m])-safeLog10(d[m-1]))/log10(g_au[m]/g_au[m-1]);
			double dist = slope * (logg-log10(g_au[m-1])) + safeLog10(d[m-1]);
			dist = pow(10,dist);
			particle.distribution.set(jRE,dist/(particle.mass*cLight2));
		},{-1,jR[DIM_R],0});
	},{0,-1,0});
	
	// TESTS
	
	// CONSERVATION OF FLUX OF PARTICLES CROSSING EACH SHELL
	if (particle.id == "ntElectron" || particle.id == "ntProton") {
		particle.ps.iterate([&](const SpaceIterator& iR) {
			double rB1 = iR.val(DIM_R) / sqrt(paso_r);
			double area = 4.0*pi*rB1*height_fun(rB1);
			double vR = abs(radialVel(rB1));
			double flux = area*vR*integSimpsonLog(particle.emin(),particle.emax(),[&iR,&particle] (double e)
							{
								return particle.distribution.interpolate({{DIM_E,e}},&iR.coord);
							},100);
			cout << iR.coord[DIM_R] << "\t flux = " << flux << endl;
		},{0,-1,0});
		
		// FLUX OF ENERGY
		particle.ps.iterate([&](const SpaceIterator& iR) {
			double rB1 = iR.val(DIM_R) / sqrt(paso_r);
			double area = 4.0*pi*rB1*height_fun(rB1);
			double vR = abs(radialVel(rB1));
			double flux = area*vR*integSimpsonLog(particle.emin(),particle.emax(),[&iR,&particle] (double e)
							{
								return particle.distribution.interpolate({{DIM_E,e}},&iR.coord)*e;
							},100);
			cout << iR.coord[DIM_R] << "\t flux of energy = " << flux << endl;
		},{0,-1,0});
		
		// COSMIC RAY PRESSURE << THERMAL PRESSURE
		particle.ps.iterate([&](const SpaceIterator& iR) {
			double r = iR.val(DIM_R);
			double pCR = integSimpsonLog(particle.emin(),particle.emax(),[&iR,&particle] (double e)
							{
								return particle.distribution.interpolate({{DIM_E,e}},&iR.coord)*e/3.0;
							},100);
			double pFluid = massDensityADAF(r)*sqrdSoundVel(r);
			cout << iR.coord[DIM_R] << "\t pCR/pgas = " << pCR/pFluid << endl;
		},{0,-1,0});
		
		// THERMAL AND NONTHERMAL TOTAL ENERGY
		particle.ps.iterate([&](const SpaceIterator& iR) {
			double r = iR.val(DIM_R);
			double U_CR = integSimpsonLog(particle.emin(),particle.emax(),[&iR,&particle] (double e)
							{
								return particle.distribution.interpolate({{DIM_E,e}},&iR.coord)*e;
							},100);
			double U_Th = 1.5*st.denf_i.get(iR)*boltzmann*st.tempIons.get(iR);
			double Press = massDensityADAF(r)*sqrdSoundVel(r);
			cout << iR.coord[DIM_R] << "\t U_Nth/U_Th = " << U_CR/U_Th << "\t U_Nth/P_gas = " << U_Th/Press << endl;
		},{0,-1,0});
	}
	
}


void distributionSecondaries(Particle& particle, State& st)
{
	particle.ps.iterate([&](const SpaceIterator& iR) {
		
		SpaceCoord jR = {0,nR-1-iR.coord[DIM_R],0};
		double r = particle.ps[DIM_R][jR[DIM_R]];
		double height = height_fun(r);
		double magf = magneticField(r);
		double rB2 = r*sqrt(paso_r);
		double rB1 = rB2/paso_r;
		double vol = volume(r);
		double rho = massDensityADAF(r);
		double rateAccretion = abs(radialVel(rB1))*4.0*pi*height_fun(rB1)*rB1/vol;
		double rateWind = (accRateADAF(rB2)-accRateADAF(rB1)) / (rho*vol);
		rateWind = 0.0;
		
		particle.ps.iterate([&](const SpaceIterator& jRE) {
			
			double E = jRE.val(DIM_E);
			double g = E / (particle.mass*cLight2);
			double rateDiff = 1.0/diffusionTimeTurbulence(E,height,particle,magf);
			double rateCool = losses(E,particle,st,jR) / E;
			double rateDecay = 0.0;
			
			if (particle.id == "ntChargedPion") rateDecay = 1.0 / (chargedPionMeanLife*g);
			if (particle.id == "ntMuon") rateDecay = 1.0 / (muonMeanLife*g);
			
			double rateTotal = rateCool + rateDecay + rateWind + rateAccretion + rateDiff;

			double Qinj = 0.0;
			double Qplus = 0.0;
			if (jR[DIM_R] < nR-1) {
				SpaceCoord jRplus = {0,jR[DIM_R]+1,0};
				double areaPlus = 4.0*pi*height_fun(rB2)*rB2;
				double Nplus = particle.distribution.interpolate({{DIM_E,E}},&jRplus);
				Qplus = Nplus * areaPlus * abs(radialVel(rB2)) / vol;
			}
			Qinj = particle.injection.interpolate({{DIM_E,E}},&jR) + Qplus;
			
			particle.distribution.set(jRE,Qinj/rateTotal);
		},{-1,jR[DIM_R],0});
	},{0,-1,0});
}


void distributionFokkerPlanckMultiZone(Particle& particle, State& st)
{
	particle.ps.iterate([&](const SpaceIterator& iR) {
		
		SpaceCoord jR = {0,nR-1-iR.coord[DIM_R],0};
		double r = particle.ps[DIM_R][jR[DIM_R]];
		double rB2 = r*sqrt(paso_r);
		double rB1 = rB2/paso_r;
		double vR = radialVel(r);
		double area = 4.0*pi*height_fun(r)*r;
		double tCell = (rB2-rB1)/abs(vR);
		double vol = volume(r);
		double height = height_fun(r);
		double B = st.magf.get(jR);
		double rho = massDensityADAF(r);
		double tAccretion = accretionTime(r);
		double rateAccretion = 1.0/min(tCell,tAccretion);
		rateAccretion = abs(radialVel(rB1))*4.0*pi*height_fun(rB1)*rB1/vol;
		//double rateWind = accRateADAF(r)*(hAcc(rB2)-hAcc(rB1)) / (rho*vol);
		double rateWind = 2.0*s*abs(radialVel(r))/r;
		
		// We define a new mesh of points
		size_t M = 199;
		double g_min = particle.emin()/(particle.mass*cLight2);
		double g_max = particle.emax()/(particle.mass*cLight2);
		
		double paso_g = pow(g_max/g_min,1.0/M);
		Vector g_au(M+3,g_min/paso_g);
		Vector delta_g_m_au(M+2,0.0);
		Vector delta_g(M+1,0.0);
		
		double g_inj = 2.0*g_min;
		double E_inj = g_inj*particle.mass*cLight2;
		
		double norm_temp, dens;
		if (particle.id == "ntProton") {
			norm_temp = boltzmann*st.tempIons.get(jR)/(particle.mass*cLight2);
			dens = st.denf_i.get(jR);
		} else if (particle.id == "ntElectron") {
			norm_temp = boltzmann*st.tempElectrons.get(jR)/(particle.mass*cLight2);
			dens = st.denf_e.get(jR);
		}
		
		// If NORMALIZATION A PRIORI:
		double Q0 = Ainjection * dens * cLight / schwRadius / g_inj;
		
		double sigma = g_inj / 10.0;
		
		for (size_t jg=1;jg<M+3;jg++)	g_au[jg] = g_au[jg-1]*paso_g;
		for (size_t jg=0;jg<M+2;jg++)	delta_g_m_au[jg] = g_au[jg+1]-g_au[jg];
		for (size_t jg=0;jg<M+1;jg++)	delta_g[jg] = 0.5*(g_au[jg+2]-g_au[jg]);
		
		fun1 Bfun = [&particle,&st,&jR,B,rho,height,vR] (double g) 
						{
							double E = g*particle.mass*cLight2;
							double tCool = E/abs(losses(E,particle,st,jR));
							double Dg = diffCoeff_g(g,particle,height,B,rho);
							return g/tCool - 2.0*Dg/g; // - g*abs(vR)/(2.0*r);
						};

		fun1 Cfun = [&particle,B,height,rho,&jR] (double g)
						{
							double Dg = diffCoeff_g(g,particle,height,B,rho);
							return Dg;
						};

		fun1 Tfun = [r,&particle,height,B,rateAccretion,&rateWind] (double g) 
						{
							double E = g*particle.mass*cLight2;
							double rateDiff = 1.0/diffusionTimeTurbulence(E,height,particle,B);
							return pow(rateAccretion+rateWind+rateDiff,-1);
						};
		
		fun1 Qfun = [Q0,g_inj,sigma,&particle,&jR,tCell,rB2,area,vR,vol] (double g)
						{
							double Qplus = 0.0;
							if (jR[DIM_R] < nR-1) {
								double E = g*particle.mass*cLight2;
								SpaceCoord jRplus = {0,jR[DIM_R]+1,0};
								double rPlus = particle.ps[DIM_R][jRplus[DIM_R]];
								double areaPlus = 4.0*pi*height_fun(rB2)*rB2;
								double vRplus = radialVel(rPlus);
								double Nplus = particle.distribution.interpolate({{DIM_E,E}},&jRplus)*(particle.mass*cLight2);
								Qplus = Nplus * areaPlus * abs(radialVel(rB2)) / vol;
							}
							return Q0*gsl_ran_gaussian_pdf(g-g_inj,sigma) + Qplus;
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
		particle.ps.iterate([&](const SpaceIterator &jRE) {
			double logg = log10(jRE.val(DIM_E)/(particle.mass*cLight2));
			while (log10(g_au[m]) < logg)
				m++;
			double slope = (safeLog10(d[m])-safeLog10(d[m-1]))/log10(g_au[m]/g_au[m-1]);
			double dist = slope * (logg-log10(g_au[m-1])) + safeLog10(d[m-1]);
			dist = pow(10,dist);
			particle.distribution.set(jRE,dist/(particle.mass*cLight2));
		},{-1,jR[DIM_R],0});
	},{0,-1,0});
	
	// NORMALIZATION A POSTERIORI
	/*
	double sum = 0.0;
	particle.ps.iterate([&](const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		double rB2 = r*sqrt(paso_r);
		double rB1 = rB2/paso_r;
		double vol = volume(r);
		double massDens = massDensityADAF(r);
		double massLostWinds = accRateADAF(rB2)-accRateADAF(rB1);
		double magf = st.magf.get(iR);
		double height = height_fun(r);
		double u_nth = integSimpsonLog(particle.emin(),particle.emax(),
					[&particle,&iR] (double e)
					{
						return particle.distribution.interpolate({{DIM_E,e}},&iR.coord)*e;
					},100);
		double q_losses = integSimpsonLog(particle.emin(),particle.emax(),
					[&particle,&iR,&height,&magf,&st] (double e)
					{
						double rateDiff = 1.0 / diffusionTimeTurbulence(e,height,particle,magf);
						rateDiff = 0.0;
						double rateCooling = losses(e,particle,st,iR)/e;
						double rateLosses = rateDiff + rateCooling;
						return particle.distribution.interpolate({{DIM_E,e}},&iR.coord)*e*rateLosses;
					},100);
		sum += vol*q_losses + u_nth * massLostWinds/massDens;
	},{0,-1,0});
	SpaceCoord iRcoord = {0,1,0};
	double u_nth_0 = integSimpsonLog(particle.emin(),particle.emax(),
					[&particle,&iRcoord] (double e)
					{
						return particle.distribution.interpolate({{DIM_E,e}},&iRcoord)*e;
					},100);
	sum += accRateADAF(particle.ps[DIM_R][1])/massDensityADAF(particle.ps[DIM_R][1])*u_nth_0;
	
	double normFactor = etaInj * accRateOut * cLight2 / sum;
	
	particle.ps.iterate([&](const SpaceIterator& i) {
		particle.distribution.set(i,particle.distribution.get(i)*normFactor);
	},{-1,-1,0});
	*/
	
	// TESTS
	
	// CONSERVATION OF FLUX OF PARTICLES CROSSING EACH SHELL
	particle.ps.iterate([&](const SpaceIterator& iR) {
		double rB1 = iR.val(DIM_R) / sqrt(paso_r);
		double area = 4.0*pi*rB1*height_fun(rB1);
		double vR = abs(radialVel(rB1));
		double flux = area*vR*integSimpsonLog(particle.emin(),particle.emax(),[&iR,&particle] (double e)
						{
							return particle.distribution.interpolate({{DIM_E,e}},&iR.coord);
						},100);
		cout << iR.coord[DIM_R] << "\t flux = " << flux << endl;
	},{0,-1,0});
	
	// FLUX OF ENERGY
	particle.ps.iterate([&](const SpaceIterator& iR) {
		double rB1 = iR.val(DIM_R) / sqrt(paso_r);
		double area = 4.0*pi*rB1*height_fun(rB1);
		double vR = abs(radialVel(rB1));
		double flux = area*vR*integSimpsonLog(particle.emin(),particle.emax(),[&iR,&particle] (double e)
						{
							return particle.distribution.interpolate({{DIM_E,e}},&iR.coord)*e;
						},100);
		cout << iR.coord[DIM_R] << "\t flux of energy = " << flux << endl;
	},{0,-1,0});
	
	// COSMIC RAY PRESSURE << THERMAL PRESSURE
	particle.ps.iterate([&](const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		double pCR = integSimpsonLog(particle.emin(),particle.emax(),[&iR,&particle] (double e)
						{
							return particle.distribution.interpolate({{DIM_E,e}},&iR.coord)*e/3.0;
						},100);
		double pFluid = massDensityADAF(r)*sqrdSoundVel(r);
		cout << iR.coord[DIM_R] << "\t pCR/pgas = " << pCR/pFluid << endl;
	},{0,-1,0});
	
}

/*
void distributionFokkerPlanckMultiZoneTimeDependent(Particle& particle, State& st)
{
	particle.ps.iterate([&](const SpaceIterator& iR) {
		
		SpaceCoord jR = {0,nR-1-iR.coord[DIM_R],0};
		double r = particle.ps[DIM_R][jR[DIM_R]];
		double rB2 = r*sqrt(paso_r);
		double rB1 = rB2/paso_r;
		double vR = radialVel(r);
		double area = 4.0*pi*height_fun(r)*r;
		double tCell = (rB2-rB1)/abs(vR);
		double vol = volume(r);
		double height = height_fun(r);
		double B = st.magf.get(jR);
		double rho = massDensityADAF(r);
		double tAccretion = accretionTime(r);
		double rateAccretion = 1.0/min(tCell,tAccretion);
		rateAccretion = abs(radialVel(rB1))*4.0*pi*height_fun(rB1)*rB1/vol;
		double rateWind = (accRateADAF(rB2)-accRateADAF(rB1)) / (rho*vol);
		
		// We define a new mesh of points
		size_t M = 99;
		double g_min = particle.emin()/(particle.mass*cLight2);
		double g_max = particle.emax()/(particle.mass*cLight2);
		
		double paso_g = pow(g_max/g_min,1.0/M);
		Vector g_au(M+3,g_min/paso_g);
		Vector delta_g_m_au(M+2,0.0);
		Vector delta_g(M+1,0.0);
		
		double g_inj = 2.0*g_min;
		double E_inj = g_inj*particle.mass*cLight2;
		
		double norm_temp, dens;
		if (particle.id == "ntProton") {
			norm_temp = boltzmann*st.tempIons.get(jR)/(particle.mass*cLight2);
			dens = st.denf_i.get(jR);
		} else if (particle.id == "ntElectron") {
			norm_temp = boltzmann*st.tempElectrons.get(jR)/(particle.mass*cLight2);
			dens = st.denf_e.get(jR);
		}
		double aTheta = 3.0 - 6.0/(4.0+5.0*norm_temp); // Gammie & Popham (1998)
		double uth = dens*norm_temp*(particle.mass*cLight2)*aTheta;   // erg cm^-3
		
		// If NORMALIZATION A PRIORI:
		double Q0 = Ainjection *dens * norm_temp;
		
		double sigma = g_inj * (paso_g-1.0) / 10;
		
		for (size_t jg=1;jg<M+3;jg++)	g_au[jg] = g_au[jg-1]*paso_g;
		for (size_t jg=0;jg<M+2;jg++)	delta_g_m_au[jg] = g_au[jg+1]-g_au[jg];
		for (size_t jg=0;jg<M+1;jg++)	delta_g[jg] = 0.5*(g_au[jg+2]-g_au[jg]);
		
		fun1 Bfun = [&particle,&st,&jR,B,rho,height] (double g) 
						{
							double E = g*particle.mass*cLight2;
							double tCool = E/abs(losses(E,particle,st,jR));
							double Dg = diffCoeff_g(g,particle,height,B,rho);
							return g/tCool - 2.0*Dg/g;
						};

		fun1 Cfun = [&particle,B,height,rho,&jR] (double g)
						{
							double Dg = diffCoeff_g(g,particle,height,B,rho);
							return Dg;
						};

		fun1 Tfun = [r,&particle,height,B,rateAccretion,&rateWind] (double g) 
						{
							double E = g*particle.mass*cLight2;
							double rateDiff = 1.0/diffusionTimeTurbulence(E,height,particle,B);
							rateDiff = 0.0;
							return pow(rateAccretion+rateWind+rateDiff,-1);
						};
		
		fun1 Qfun = [Q0,g_inj,sigma,&particle,&jR,tCell,rB2,area,vR,vol] (double g)
						{
							double Qplus = 0.0;
							//double Qplus = Q0*gsl_ran_gaussian_pdf(g-g_inj,sigma);
							if (jR[DIM_R] < nR-1) {
								double E = g*particle.mass*cLight2;
								SpaceCoord jRplus = {0,jR[DIM_R]+1,0};
								double rPlus = particle.ps[DIM_R][jRplus[DIM_R]];
								double areaPlus = 4.0*pi*height_fun(rB2)*rB2;
								double vRplus = radialVel(rPlus);
								double factor = (areaPlus/area) * abs(vRplus/vR);
								double Nplus = particle.distribution.interpolate({{DIM_E,E}},&jRplus)*(particle.mass*cLight2);
								Qplus = Nplus * areaPlus * abs(radialVel(rB2)) / vol;
							}
							return Q0*gsl_ran_gaussian_pdf(g-g_inj,sigma) + Qplus;
							//return Qplus;
						};
							
		size_t Ntime = 200;
		double dt = g_min*g_min/diffCoeff_g(g_min,particle,height,B,rho);
		Vector time(Ntime+1,0.0);
		time[0] = 0.0;
		time[1] = dt;
		double pasot = pow(10000*tCell/time[1],1.0/(Ntime-1));
		Vector deltat(Ntime,0.0);
		for (size_t jt=2;jt<Ntime;jt++) time[jt] = time[jt-1]*pasot;
		for (size_t jt=0;jt<Ntime;jt++) deltat[jt] = (time[jt+1]-time[jt]);

		Vector d(M+1,0.0);
		ofstream file;
		file.open("fokkerPlanck.dat");
		for (size_t m=0;m<=M;m++)
			file << safeLog10(g_au[m]) << "\t" << Qfun(g_au[m]) << endl;
		file.close();

		for (size_t jt=0;jt<Ntime;jt++) {
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
					a[m] = - deltat[jt] * Cm_minusHalf * Wminus_m_minusHalf / ( delta_g[m] * delta_g_m_au[m] );
				if (m <= M-1)
					c[m] = - deltat[jt] * Cm_plusHalf * Wplus_m_plusHalf / ( delta_g[m] * delta_g_m_au[m+1] );
				b[m] = 1.0 + deltat[jt]/delta_g[m] * ( Cm_minusHalf * Wplus_m_minusHalf / delta_g_m_au[m] +
							Cm_plusHalf * Wminus_m_plusHalf / delta_g_m_au[m+1] ) + deltat[jt]/Tm;
				d[m] = Qm*deltat[jt] + d[m];
				
			}
			TriDiagSys(a,b,c,d,M);
		}
		
		size_t m = 0;
		particle.ps.iterate([&](const SpaceIterator &jRE) {
			double logg = log10(jRE.val(DIM_E)/(particle.mass*cLight2));
			while (log10(g_au[m]) < logg)
				m++;
			double slope = (safeLog10(d[m])-safeLog10(d[m-1]))/log10(g_au[m]/g_au[m-1]);
			double dist = slope * (logg-log10(g_au[m-1])) + safeLog10(d[m-1]);
			dist = pow(10,dist);
			particle.distribution.set(jRE,dist/(particle.mass*cLight2));
		},{-1,jR[DIM_R],0});
		
		m = 0;
		particle.ps.iterate([&](const SpaceIterator &jRE) {
			particle.distribution.set(jRE,d[m]/(particle.mass*cLight2));
			m++;
		},{-1,jR[DIM_R],0});
		
		file.open("fokkerPlanck.dat");
		for (size_t m=0;m<=M;m++) {
			double E = g_au[m+1]*particle.mass*cLight2;
			double E_GeV = E / 1.602e-3;
			file << safeLog10(E_GeV) << "\t" << vol*d[m]*E*g_au[m+1] << endl;
		}
		file.close();

	},{0,-1,0});
	
	// NORMALIZATION
	/*
	double sum = 0.0;
	particle.ps.iterate([&](const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		double rB2 = r*sqrt(paso_r);
		double rB1 = rB2/paso_r;
		double vol = volume(r);
		double massDens = massDensityADAF(r);
		double massLostWinds = accRateADAF(rB2)-accRateADAF(rB1);
		double magf = st.magf.get(iR);
		double height = height_fun(r);
		double u_nth = integSimpsonLog(particle.emin(),particle.emax(),
					[&particle,&iR] (double e)
					{
						return particle.distribution.interpolate({{DIM_E,e}},&iR.coord)*e;
					},100);
		double q_losses = integSimpsonLog(particle.emin(),particle.emax(),
					[&particle,&iR,&height,&magf,&st] (double e)
					{
						double rateDiff = 1.0 / diffusionTimeTurbulence(e,height,particle,magf);
						rateDiff = 0.0;
						double rateCooling = losses(e,particle,st,iR)/e;
						double rateLosses = rateDiff + rateCooling;
						return particle.distribution.interpolate({{DIM_E,e}},&iR.coord)*e*rateLosses;
					},100);
		sum += vol*q_losses + u_nth * massLostWinds/massDens;
	},{0,-1,0});
	SpaceCoord iRcoord = {0,1,0};
	double u_nth_0 = integSimpsonLog(particle.emin(),particle.emax(),
					[&particle,&iRcoord] (double e)
					{
						return particle.distribution.interpolate({{DIM_E,e}},&iRcoord)*e;
					},100);
	sum += accRateADAF(particle.ps[DIM_R][1])/massDensityADAF(particle.ps[DIM_R][1])*u_nth_0;
	
	double normFactor = etaInj * accRateOut * cLight2 / sum;
	
	particle.ps.iterate([&](const SpaceIterator& i) {
		particle.distribution.set(i,particle.distribution.get(i)*normFactor);
	},{-1,-1,0});
	*/
	
	// TESTS
	/*
	// CONSERVATION OF FLUX OF PARTICLES CROSSING EACH SHELL
	particle.ps.iterate([&](const SpaceIterator& iR) {
		double rB1 = iR.val(DIM_R) / sqrt(paso_r);
		double area = 4.0*pi*rB1*height_fun(rB1);
		double vR = abs(radialVel(rB1));
		double flux = area*vR*integSimpsonLog(particle.emin(),particle.emax(),[&iR,&particle] (double e)
						{
							return particle.distribution.interpolate({{DIM_E,e}},&iR.coord);
						},100);
		cout << iR.coord[DIM_R] << "\t flux = " << flux << endl;
	},{0,-1,0});
	
	// FLUX OF ENERGY
	particle.ps.iterate([&](const SpaceIterator& iR) {
		double rB1 = iR.val(DIM_R) / sqrt(paso_r);
		double area = 4.0*pi*rB1*height_fun(rB1);
		double vR = abs(radialVel(rB1));
		double flux = area*vR*integSimpsonLog(particle.emin(),particle.emax(),[&iR,&particle] (double e)
						{
							return particle.distribution.interpolate({{DIM_E,e}},&iR.coord)*e;
						},100);
		cout << iR.coord[DIM_R] << "\t flux of energy = " << flux << endl;
	},{0,-1,0});
	
	// COSMIC RAY PRESSURE << THERMAL PRESSURE
	particle.ps.iterate([&](const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		double pCR = integSimpsonLog(particle.emin(),particle.emax(),[&iR,&particle] (double e)
						{
							return particle.distribution.interpolate({{DIM_E,e}},&iR.coord)*e/3.0;
						},100);
		double pFluid = massDensityADAF(r)*sqrdSoundVel(r);
		cout << iR.coord[DIM_R] << "\t pCR/pgas = " << pCR/pFluid << endl;
	},{0,-1,0});
}*/



void distributionFokkerPlanckRadial(Particle& particle, State& st)
{
	
	// We define a new mesh of points
	size_t M = 199;
	double g_min = particle.emin()/(particle.mass*cLight2);
	double g_max = particle.emax()/(particle.mass*cLight2);
	
	double paso_g = pow(g_max/g_min,1.0/M);
	Vector g_au(M+3,g_min/paso_g);
	Vector delta_g_m_au(M+2,0.0);
	Vector delta_g(M+1,0.0);
	
	double g_inj = 2.0*g_min;
	double E_inj = g_inj*particle.mass*cLight2;
	
	double sigma = g_inj / 10.0;
	
	for (size_t jg=1;jg<M+3;jg++)	g_au[jg] = g_au[jg-1]*paso_g;
	for (size_t jg=0;jg<M+2;jg++)	delta_g_m_au[jg] = g_au[jg+1]-g_au[jg];
	for (size_t jg=0;jg<M+1;jg++)	delta_g[jg] = 0.5*(g_au[jg+2]-g_au[jg]);
	
	fun2 Bfun = [&particle,&st] (double g, double r) 
					{
						double E = g*particle.mass*cLight2;
						double rad = r;
						if (r > particle.ps[DIM_R][nR-1]) rad = particle.ps[DIM_R][nR-1]/1.1;
						if (r < particle.ps[DIM_R][0]) rad = particle.ps[DIM_R][0]*1.1;
						
						double vR = radialVel(rad);
						double B = magneticField(rad);
						double height = height_fun(rad);
						double rho = massDensityADAF(rad);
						
						
						size_t jjR=0;
						while (particle.ps[DIM_R][jjR] < rad) jjR++;
						
						SpaceCoord jR1 = {0,jjR-1,0};
						SpaceCoord jR2 = {0,jjR,0};
						double logdEdt1 = safeLog10(abs(losses(E,particle,st,jR1)));
						double logdEdt2 = safeLog10(abs(losses(E,particle,st,jR2)));
						double deltaLogR = safeLog10(particle.ps[DIM_R][jjR]/particle.ps[DIM_R][jjR-1]);
						double slope = (logdEdt2 - logdEdt1)/deltaLogR;
						double logdEdt = logdEdt1 + slope * (safeLog10(rad) - safeLog10(particle.ps[DIM_R][jjR-1]));
						
						double dgdt = - pow(10,logdEdt) / (particle.mass*cLight2);
						
						/*
						size_t jjR=0;
						particle.ps.iterate([&](const SpaceIterator& iR) {
							double rAux = iR.val(DIM_R);
							if (rAux < rad) jjR++;
						},{0,-1,0});
						SpaceCoord jR = {0,jjR,0};
						double dgdt = - abs(losses(E,particle,st,jR)) / (particle.mass*cLight2);
						*/
						
						double Dg = diffCoeff_g(g,particle,height,B,rho);
						return - dgdt/vR - 2.0*Dg/vR / g + 0.5*g/r;
						//return 0.0;
					};
					
	fun2 Cfun = [&particle] (double g, double r)
					{
						double rad = r;
						if (r > particle.ps[DIM_R][nR-1]) rad = particle.ps[DIM_R][nR-1]/1.1;
						if (r < particle.ps[DIM_R][0]) rad = particle.ps[DIM_R][0]*1.1;
						
						double B = magneticField(rad);
						double height = height_fun(rad);
						double rho = massDensityADAF(rad);
						double vR = radialVel(rad);
						double Dg = diffCoeff_g(g,particle,height,B,rho);
						return Dg / vR;
					};

	fun2 Rfun = [&particle] (double g, double r) 
					{
						double rad = r;
						if (r > particle.ps[DIM_R][nR-1]) rad = particle.ps[DIM_R][nR-1]/1.11;
						if (r < particle.ps[DIM_R][0]) rad = particle.ps[DIM_R][0]*1.1;
						
						double E = g*particle.mass*cLight2;
						double vR = radialVel(rad);
						double rateWind = 2.0*s*abs(vR) / rad;
						double B = magneticField(rad);
						double height = height_fun(rad);
						double rateDiff = 1.0/(diffusionTimeTurbulence(E,height,particle,B));
						//return 1.0e99;
						return pow(rateWind + rateDiff,-1) * vR;
					};
	
	fun2 Qfun = [g_inj,sigma,&particle] (double g, double r)
					{
						double rad = r;
						if (r > particle.ps[DIM_R][nR-1]) rad = particle.ps[DIM_R][nR-1]/1.11;
						if (r < particle.ps[DIM_R][0]) rad = particle.ps[DIM_R][0]*1.1;
						
						double norm_temp, dens;
						if (particle.id == "ntProton") {
							norm_temp = boltzmann*ionTemp(rad)/(particle.mass*cLight2);
							dens = electronDensity(rad)*eMeanMolecularWeight/iMeanMolecularWeight;
						} else if (particle.id == "ntElectron") {
							norm_temp = boltzmann*electronTemp(rad)/(particle.mass*cLight2);
							dens = electronDensity(rad);
						}
						double vR = radialVel(rad);
						double vA = magneticField(rad)/sqrt(4.0*pi*massDensityADAF(rad));
						double height = height_fun(rad);
						double Q0 = 1.0e-8*Ainjection * dens * cLight/schwRadius / g_inj;
						//return (rad*sqrt(paso_r) > 70*schwRadius) ? rad*height*Q0*gsl_ran_gaussian_pdf(g-g_inj,sigma) : 0.0;
						return rad*height*Q0*gsl_ran_gaussian_pdf(g-g_inj,sigma);
					};
		
	size_t Nrad = 100;
	Vector radius(Nrad+1,0.0);
	radius[0] = particle.ps[DIM_R][nR-1];
	double pasor = pow(radius[0]/particle.ps[DIM_R][0],1.0/(Nrad-1));
	Vector deltar(Nrad,0.0);
	for (size_t jr=1;jr<Nrad+1;jr++) radius[jr] = radius[jr-1]/pasor;
	for (size_t jr=0;jr<Nrad;jr++) deltar[jr] = (radius[jr+1]-radius[jr]);

	Vector d(M+1,0.0);
	Matrix dist;
	matrixInit(dist,Nrad,M+1,0.0);
	for (size_t jr=0;jr<Nrad;jr++) {
		Vector a(M+1,0.0), b(M+1,0.0), c(M+1,0.0);
		for (size_t m=0;m<M+1;m++) {
			double Bm_minusHalf = 0.5 * (Bfun(g_au[m],radius[jr]) + Bfun(g_au[m+1],radius[jr]));
			double Bm_plusHalf = 0.5 * (Bfun(g_au[m+1],radius[jr]) + Bfun(g_au[m+2],radius[jr]));
			double Cm_minusHalf = 0.5 * (Cfun(g_au[m],radius[jr]) + Cfun(g_au[m+1],radius[jr]));
			double Cm_plusHalf = 0.5 * (Cfun(g_au[m+1],radius[jr]) + Cfun(g_au[m+2],radius[jr]));
			double Qm = Qfun(g_au[m+1],radius[jr]);
			double Rm = Rfun(g_au[m+1],radius[jr]);
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
				a[m] = - deltar[jr] * Cm_minusHalf * Wminus_m_minusHalf / ( delta_g[m] * delta_g_m_au[m] );
			if (m <= M-1)
				c[m] = - deltar[jr] * Cm_plusHalf * Wplus_m_plusHalf / ( delta_g[m] * delta_g_m_au[m+1] );
			b[m] = 1.0 + deltar[jr]/delta_g[m] * ( Cm_minusHalf * Wplus_m_minusHalf / delta_g_m_au[m] +
						Cm_plusHalf * Wminus_m_plusHalf / delta_g_m_au[m+1] ) + deltar[jr]/Rm;
			d[m] = Qm*deltar[jr] + d[m];
		}
		TriDiagSys(a,b,c,d,M);
		for (size_t m=0;m<M+1;m++) dist[jr][m] = d[m]/(radius[jr]*height_fun(radius[jr])*radialVel(radius[jr]));
	}
	
	size_t jr(Nrad);
	particle.ps.iterate([&](const SpaceIterator& jR) {
		double logr = log10(jR.val(DIM_R));
		while (log10(radius[jr]) < logr) jr--;
		if (jr >= Nrad-1) jr = Nrad-2;
		size_t m(1);
		particle.ps.iterate([&](const SpaceIterator &jRE) {
			double logg = log10(jRE.val(DIM_E)/(particle.mass*cLight2));
			while (log10(g_au[m+1]) < logg) m++;
			
			double slope1 = (safeLog10(dist[jr][m])-safeLog10(dist[jr][m-1]))/log10(g_au[m+1]/g_au[m]);
			double dist_1 = slope1 * (logg-log10(g_au[m])) + safeLog10(dist[jr][m-1]);
			double slope2 = (safeLog10(dist[jr+1][m])-safeLog10(dist[jr+1][m-1]))/log10(g_au[m+1]/g_au[m]);
			double dist_2 = slope2 * (logg-log10(g_au[m])) + safeLog10(dist[jr+1][m-1]);
			double slope = (safeLog10(dist_2)-safeLog10(dist_1))/log10(radius[jr+1]/radius[jr]);
			double distr = slope * (logr-log10(radius[jr+1])) + dist_2;

			if (m < M+1) distr = pow(10,distr);
			else distr = 0.0;
			
			particle.distribution.set(jRE,distr/(particle.mass*cLight2));
		},{-1,jR.coord[DIM_R],0});
	},{0,-1,0});
	
	// NORMALIZATION
	double sum = 0.0;
	particle.ps.iterate([&](const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		double rB2 = r*sqrt(paso_r);
		double rB1 = rB2/paso_r;
		double vol = volume(r);
		double massDens = massDensityADAF(r);
		double massLostWinds = accRateADAF(r) * 2.0*s* (paso_r - 1.0);
		double magf = st.magf.get(iR);
		double height = height_fun(r);
		double u_nth = integSimpsonLog(particle.emin(),particle.emax(),
					[&particle,&iR] (double e)
					{
						return particle.distribution.interpolate({{DIM_E,e}},&iR.coord)*e;
					},100);
		double q_losses = integSimpsonLog(particle.emin(),particle.emax(),
					[&particle,&iR,&height,&magf,&st] (double e)
					{
						double rateDiff = 1.0 / diffusionTimeTurbulence(e,height,particle,magf);
						rateDiff = 0.0;
						double rateCooling = losses(e,particle,st,iR)/e;
						double rateLosses = rateDiff + rateCooling;
						return particle.distribution.interpolate({{DIM_E,e}},&iR.coord)*e*rateLosses;
					},100);
		sum += vol*q_losses + u_nth * massLostWinds/massDens;
	},{0,-1,0});
	SpaceCoord iRcoord = {0,1,0};
	double u_nth_0 = integSimpsonLog(particle.emin(),particle.emax(),
					[&particle,&iRcoord] (double e)
					{
						return particle.distribution.interpolate({{DIM_E,e}},&iRcoord)*e;
					},100);
	sum += accRateADAF(particle.ps[DIM_R][1])/massDensityADAF(particle.ps[DIM_R][1])*u_nth_0;
	
	double normFactor = etaInj * accRateOut * cLight2 / sum;
	
	particle.ps.iterate([&](const SpaceIterator& i) {
		particle.distribution.set(i,particle.distribution.get(i)*normFactor);
	},{-1,-1,0});
	
	
	// TESTS
	
	// CONSERVATION OF FLUX OF PARTICLES CROSSING EACH SHELL
	particle.ps.iterate([&](const SpaceIterator& iR) {
		double rB1 = iR.val(DIM_R) / sqrt(paso_r);
		double area = 4.0*pi*rB1*height_fun(rB1);
		double vR = abs(radialVel(rB1));
		double flux = area*vR*integSimpsonLog(particle.emin(),particle.emax(),[&iR,&particle] (double e)
						{
							return particle.distribution.interpolate({{DIM_E,e}},&iR.coord);
						},100);
		cout << iR.coord[DIM_R] << "\t flux = " << flux << endl;
	},{0,-1,0});
	
	// FLUX OF ENERGY
	particle.ps.iterate([&](const SpaceIterator& iR) {
		double rB1 = iR.val(DIM_R) / sqrt(paso_r);
		double area = 4.0*pi*rB1*height_fun(rB1);
		double vR = abs(radialVel(rB1));
		double flux = area*vR*integSimpsonLog(particle.emin(),particle.emax(),[&iR,&particle] (double e)
						{
							return particle.distribution.interpolate({{DIM_E,e}},&iR.coord)*e;
						},100);
		cout << iR.coord[DIM_R] << "\t flux of energy = " << flux << endl;
	},{0,-1,0});
	
	// COSMIC RAY PRESSURE << THERMAL PRESSURE
	particle.ps.iterate([&](const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		double pCR = integSimpsonLog(particle.emin(),particle.emax(),[&iR,&particle] (double e)
						{
							return particle.distribution.interpolate({{DIM_E,e}},&iR.coord)*e/3.0;
						},100);
		double pFluid = massDensityADAF(r)*sqrdSoundVel(r);
		double uTh = 1.5 * st.denf_i.get(iR) * boltzmann * st.tempIons.get(iR);
		cout << iR.coord[DIM_R] << "\t pCR/pgas = " << pCR/pFluid << "\t uNth/uth = " << 3.0*pCR/uTh << endl;
	},{0,-1,0});
	
}



void distributionMultiZoneRadial(Particle& particle, State& st)
{
	
	// We define a new mesh of points
	size_t M = 299;
	double g_min = particle.emin()/(particle.mass*cLight2);
	double g_max = particle.emax()/(particle.mass*cLight2);
	
	double paso_g = pow(g_max/g_min,1.0/M);
	Vector g_au(M+3,g_min/paso_g);
	Vector delta_g_m_au(M+2,0.0);
	Vector delta_g(M+1,0.0);
	
	for (size_t jg=1;jg<M+3;jg++)	g_au[jg] = g_au[jg-1]*paso_g;
	for (size_t jg=0;jg<M+2;jg++)	delta_g_m_au[jg] = g_au[jg+1]-g_au[jg];
	for (size_t jg=0;jg<M+1;jg++)	delta_g[jg] = 0.5*(g_au[jg+2]-g_au[jg]);
	
	fun2 Bfun = [&particle,&st] (double g, double r) 
					{
						double E = g*particle.mass*cLight2;
						double rad = r;
						if (r > particle.ps[DIM_R][nR-1]) rad = particle.ps[DIM_R][nR-1]/1.1;
						if (r < particle.ps[DIM_R][0]) rad = particle.ps[DIM_R][0]*1.1;
						
						size_t jjR=0;
						while (particle.ps[DIM_R][jjR] < rad) jjR++;
						
						SpaceCoord jR1 = {0,jjR-1,0};
						SpaceCoord jR2 = {0,jjR,0};
						double vR = radialVel(rad);
						double logdEdt1 = safeLog10(abs(losses(E,particle,st,jR1)));
						double logdEdt2 = safeLog10(abs(losses(E,particle,st,jR2)));
						double deltaLogR = safeLog10(particle.ps[DIM_R][jjR]/particle.ps[DIM_R][jjR-1]);
						double slope = (logdEdt2 - logdEdt1)/deltaLogR;
						double logdEdt = logdEdt1 + slope * (safeLog10(rad) - safeLog10(particle.ps[DIM_R][jjR-1]));
						
						double dgdt = - pow(10,logdEdt) / (particle.mass*cLight2);
						return -  dgdt / vR;
					};

	fun2 Rfun = [&particle] (double g, double r) 
					{
						double rad = r;
						if (r > particle.ps[DIM_R][nR-1]) rad = particle.ps[DIM_R][nR-1]/1.11;
						if (r < particle.ps[DIM_R][0]) rad = particle.ps[DIM_R][0]*1.1;
						
						double E = g*particle.mass*cLight2;
						double vR = radialVel(rad);
						double rateWind = s/rad;
						double B = magneticField(rad);
						double height = height_fun(rad);
						double rateDiff = 1.0/(diffusionTimeTurbulence(E,height,particle,B));
						return pow(rateDiff + rateWind,-1) * vR;
					};
	
	fun2 Qfun = [&particle] (double g, double r)
					{
						double E = g * particle.mass * cLight2;
						double rad = r;
						if (r > particle.ps[DIM_R][nR-1]) rad = particle.ps[DIM_R][nR-1]/1.11;
						if (r < particle.ps[DIM_R][0]) rad = particle.ps[DIM_R][0]*1.1;
						
						size_t jjR=0;
						particle.ps.iterate([&](const SpaceIterator& iR) {
							double rAux = iR.val(DIM_R);
							if (rAux < rad) jjR++;
						},{0,-1,0});
						SpaceCoord jR = {0,jjR,0};
						
						double norm_temp, dens;
						if (particle.id == "ntProton") {
							norm_temp = boltzmann*ionTemp(rad)/(particle.mass*cLight2);
							dens = electronDensity(rad)*eMeanMolecularWeight/iMeanMolecularWeight;
						} else if (particle.id == "ntElectron") {
							norm_temp = boltzmann*electronTemp(rad)/(particle.mass*cLight2);
							dens = electronDensity(rad);
						}
						double height = height_fun(rad);
						return (rad > 70*schwRadius) ? particle.injection.interpolate({{DIM_E,E},{DIM_R,rad}},&jR)*(particle.mass*cLight2)*rad*height : 0.0;
					};
		
	size_t Nrad = 50;
	Vector radius(Nrad+1,0.0);
	radius[0] = particle.ps[DIM_R][nR-1];
	double pasor = pow(radius[0]/particle.ps[DIM_R][0],1.0/(Nrad-1));
	Vector deltar(Nrad,0.0);
	for (size_t jr=1;jr<Nrad+1;jr++) radius[jr] = radius[jr-1]/pasor;
	for (size_t jr=0;jr<Nrad;jr++) deltar[jr] = (radius[jr+1]-radius[jr]);

	Vector d(M+1,0.0);
	Matrix dist;
	matrixInit(dist,Nrad,M+1,0.0);
	for (size_t jr=0;jr<Nrad;jr++) {
		Vector a(M+1,0.0), b(M+1,0.0), c(M+1,0.0);
		for (size_t m=0;m<M+1;m++) {
			double Bm_minusHalf = 0.5 * (Bfun(g_au[m],radius[jr]) + Bfun(g_au[m+1],radius[jr]));
			double Bm_plusHalf = 0.5 * (Bfun(g_au[m+1],radius[jr]) + Bfun(g_au[m+2],radius[jr]));

			double Qm = Qfun(g_au[m+1],radius[jr]);
			double Rm = Rfun(g_au[m+1],radius[jr]);
			
			if (m <= M-1)
				c[m] = - deltar[jr] * Bm_plusHalf / delta_g[m];
			b[m] = 1.0 + deltar[jr] * Bm_minusHalf / delta_g[m] + deltar[jr]/Rm;
			d[m] = Qm*deltar[jr] + d[m];
		}
		TriDiagSys(a,b,c,d,M);
		for (size_t m=0;m<M+1;m++) dist[jr][m] = d[m]/(radius[jr]*height_fun(radius[jr])*radialVel(radius[jr]));
	}
	
	size_t jr(Nrad);
	particle.ps.iterate([&](const SpaceIterator& jR) {
		double logr = log10(jR.val(DIM_R));
		while (log10(radius[jr]) < logr) jr--;
		if (jr >= Nrad-1) jr = Nrad-2;
		size_t m(1);
		particle.ps.iterate([&](const SpaceIterator &jRE) {
			double logg = log10(jRE.val(DIM_E)/(particle.mass*cLight2));
			while (log10(g_au[m+1]) < logg) m++;
			
			double slope1 = (safeLog10(dist[jr][m])-safeLog10(dist[jr][m-1]))/log10(g_au[m+1]/g_au[m]);
			double dist_1 = slope1 * (logg-log10(g_au[m])) + safeLog10(dist[jr][m-1]);
			double slope2 = (safeLog10(dist[jr+1][m])-safeLog10(dist[jr+1][m-1]))/log10(g_au[m+1]/g_au[m]);
			double dist_2 = slope2 * (logg-log10(g_au[m])) + safeLog10(dist[jr+1][m-1]);
			double slope = (safeLog10(dist_2)-safeLog10(dist_1))/log10(radius[jr+1]/radius[jr]);
			double distr = slope * (logr-log10(radius[jr+1])) + dist_2;

			if (m < M+1) distr = pow(10,distr);
			else distr = 0.0;
			
			particle.distribution.set(jRE,distr/(particle.mass*cLight2));
		},{-1,jR.coord[DIM_R],0});
	},{0,-1,0});
	
	// TESTS
	
	// CONSERVATION OF FLUX OF PARTICLES CROSSING EACH SHELL
	particle.ps.iterate([&](const SpaceIterator& iR) {
		double rB1 = iR.val(DIM_R) / sqrt(paso_r);
		double area = 4.0*pi*rB1*height_fun(rB1);
		double vR = abs(radialVel(rB1));
		double flux = area*vR*integSimpsonLog(particle.emin(),particle.emax(),[&iR,&particle] (double e)
						{
							return particle.distribution.interpolate({{DIM_E,e}},&iR.coord);
						},100);
		cout << iR.coord[DIM_R] << "\t flux = " << flux << endl;
	},{0,-1,0});
	
	// FLUX OF ENERGY
	particle.ps.iterate([&](const SpaceIterator& iR) {
		double rB1 = iR.val(DIM_R) / sqrt(paso_r);
		double area = 4.0*pi*rB1*height_fun(rB1);
		double vR = abs(radialVel(rB1));
		double flux = area*vR*integSimpsonLog(particle.emin(),particle.emax(),[&iR,&particle] (double e)
						{
							return particle.distribution.interpolate({{DIM_E,e}},&iR.coord)*e;
						},100);
		cout << iR.coord[DIM_R] << "\t flux of energy = " << flux << endl;
	},{0,-1,0});
	
	// COSMIC RAY PRESSURE << THERMAL PRESSURE
	particle.ps.iterate([&](const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		double pCR = integSimpsonLog(particle.emin(),particle.emax(),[&iR,&particle] (double e)
						{
							return particle.distribution.interpolate({{DIM_E,e}},&iR.coord)*e/3.0;
						},100);
		double pFluid = massDensityADAF(r)*sqrdSoundVel(r);
		cout << iR.coord[DIM_R] << "\t pCR/pgas = " << pCR/pFluid << endl;
	},{0,-1,0});
	
}


/*
#include "Particles.h"
#include "Radiation.h"
void distributionGAMERA(Particle& p, State& st)
{
	ofstream file;
	file.open("distributionGAMERA.dat");
	
	Vector t(nR,0.001);
	Vector B(nR,0.0);
	Vector Q(nR,0.0);
	Vector e(nE,0.0);
	
	p.ps.iterate([&](const SpaceIterator& iE) {
		e[iRE.coord[DIM_E]] = iE.val(DIM_E);
	},{-1,0,0});
	
	for (iR=0;iR<nR-1;iR++) {
		size_t jR = nR-1-iR;
		double r = p.ps[DIM_R][jR];
		double dr = r*(paso_r-1.0);
		double dt = dr/abs(radialVel(r));
		t[iR+1] = t[iR]+dT; 
		B[iR+1] = magneticField(r);
		p.ps.iterate([&](const SpaceIterator& iRE) {
			Q[][iRE.coord[DIM_E]] = p.injection.get(iRE);
		},{-1,iR,0});
	},{0,-1,0});
	
	p.ps.iterate([&](const SpaceIterator& iR) {
		Vector e(nE,0.0);
		Matrix q;
		matrixInit(q,nE,)
		Particles *proton = new Particles();
		Matrix eq;
		matrixInit(eq,nE,2,0.0);
		matrixInitTwoVec(eq,nE,e,q);
		electron -> SetCustomInjectionSpectrum(eq);
		electron -> SetBField(st.magf.get(iR));
		//electron -> SetConstantEscapeTime(tacc);
		electron -> SetTmin(1e-10);
		electron -> SetAge(1.0e-8);
		electron -> SetSolverMethod(1);
		electron -> CalculateElectronSpectrum();
		Matrix sed = electron -> GetParticleSpectrum();
		p.ps.iterate([&](const SpaceIterator& iRE) {
			file << safeLog10(sed[iRE.coord[DIM_E]][0]) << "\t" << safeLog10(sed[iRE.coord[DIM_E]][1]) << endl;
		},{-1,iR.coord[DIM_R],0});
	},{0,-1,0});
	file.close();
}
*/