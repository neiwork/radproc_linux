#include "NTdistribution.h"
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
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_randist.h>

#include <fparameters/parameters.h>
#include <fmath/tridiagonal.h>
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
/*
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
}*/

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
// This routine solves the transport equation by considering separated one zones at each
// cell in the RIAF. The escape time includes the time the fluid takes to cross the cell
// and the injection includes the flux of particles coming from the immediately outer cell.
// It is solved backward (from the outermost radius to the innermost one) and so it only
// takes into account spatial advection toward the center and it does not handle spatial
// diffusion. At each cell, the steady transport equation is a simple ordinary diff. equation
// in the energy, and it is solved by finite differences (CC70, PP96).
// Problems: it does not conserve the number of particles.
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
		double rateWind = s*abs(vR)/r;
		
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
		
		fun1 Bfun = [&particle,&st,&jR,vR,r] (double g) 
						{
							double E = g*particle.mass*cLight2;
							double tCool = E/losses(E,particle,st,jR);
							return g/tCool;
						};

		fun1 Tfun = [r,&particle,height,B,&rateAccretion,&rateWind] (double g) 
						{
							double E = g*particle.mass*cLight2;
							double rateDiff = 1.0/diffusionTimeTurbulence(E,height,particle,B);
							double rateDecay = 0.0;
							if (particle.id == "ntChargedPion")
								rateDecay = 1.0 / (g*chargedPionMeanLife);
							else if (particle.id == "ntMuon")
								rateDecay = 1.0 / (g*muonMeanLife);
							//return 1.0e99;
							rateAccretion = 0.0;
							return pow(rateAccretion+rateWind+rateDiff,-1);
						};
		
		fun1 Qfun = [&particle,&jR,rB2,vol] (double g)
						{
							double E = g*particle.mass*cLight2;
							double Qplus = 0.0;
							if (jR[DIM_R] < nR-1) {
								SpaceCoord jRplus = {0,jR[DIM_R]+1,0};
								double areaPlus = 4.0*pi*height_fun(rB2)*rB2;
								//double Nplus = particle.distribution.interpolate({{DIM_E,E}},&jRplus)*(particle.mass*cLight2);
								//Qplus = Nplus * areaPlus * abs(radialVel(rB2)) / vol;
							}
							double Qlocal = (E > particle.emin() && E < particle.emax() ?
									particle.injection.interpolate({{DIM_E,E}},&jR)*(particle.mass*cLight2) : 0.0);
							return Qlocal + Qplus;
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
		
		TriDiagSys(a,b,c,d,M+1);
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
		if (particle.id == "ntProton") {
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
		
		// THERMAL AND NONTHERMAL TOTAL ENERGY
		particle.ps.iterate([&](const SpaceIterator& iR) {
			double r = iR.val(DIM_R);
			double U_Nth = integSimpsonLog(particle.emin(),particle.emax(),[&iR,&particle] (double e)
							{
								return particle.distribution.interpolate({{DIM_E,e}},&iR.coord)*e;
							},100);
			double temp = (particle.id == "ntElectron") ? st.tempElectrons.get(iR) : st.tempIons.get(iR);
			double dens = (particle.id == "ntElectron") ? st.denf_e.get(iR) : st.denf_i.get(iR);
			double norm_temp = boltzmann*temp / (particle.mass*cLight2);
			double aTheta = 3.0 - 6.0/(4.0+5.0*norm_temp); // Gammie & Popham (1998)
			double U_Th = aTheta * dens * particle.mass*cLight2 * norm_temp;
			double Press = massDensityADAF(r)*sqrdSoundVel(r);
			cout << iR.coord[DIM_R] << "\t U_Nth/U_Th = " << U_Nth/U_Th << endl;
		},{0,-1,0});
	}
	
}

void distributionSecondaries(Particle& particle, State& st)
// This routine solves the transport equation for secondary particlesby considering separated 
// one zones at each cell in the RIAF. The escape time includes the decay timescale of the
// particles, and it consider the losses as cathastropic. The injection includes the flux of 
// particles coming from the immediately outer cell, though it plays no role.
// It is solved backward (from the outermost radius to the innermost one) and so it only
// takes into account spatial advection toward the center and it does not handle spatial
// diffusion. At each cell, the steady transport equation is an algebraic equation and it is
// easily solved.
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
			double rateAnn = 0.0;
			
			if (particle.id == "ntChargedPion") rateDecay = 1.0 / (chargedPionMeanLife*g);
			if (particle.id == "ntMuon") rateDecay = 1.0 / (muonMeanLife*g);
			if (particle.id == "ntPair") rateAnn = annihilationRate(E,st.ntPair,iR);
			
			double rateTotal = rateCool + rateDecay + rateWind + rateAccretion + rateDiff + rateAnn;

			double Qinj = 0.0;
			double Qplus = 0.0;
			if (jR[DIM_R] < nR-1) {
				SpaceCoord jRplus = {0,jR[DIM_R]+1,0};
				double areaPlus = 4.0*pi*height_fun(rB2)*rB2;
				double Nplus = particle.distribution.interpolate({{DIM_E,E}},&jRplus);
				Qplus = Nplus * areaPlus * abs(radialVel(rB2)) / vol;
				Qplus = 0.0;
			}
			Qinj = particle.injection.interpolate({{DIM_E,E}},&jR) + Qplus;
			
			particle.distribution.set(jRE,Qinj/rateTotal);
		},{-1,jR[DIM_R],0});
	},{0,-1,0});
}


void distributionFokkerPlanckMultiZone(Particle& particle, State& st)
// This routine solves the transport equation by considering separated one zones at each
// cell in the RIAF. The escape time includes the time the fluid takes to cross the cell
// and the injection includes the flux of particles coming from the immediately outer cell.
// It is solved backward (from the outermost radius to the innermost one) and so it only
// takes into account spatial advection toward the center and it does not handle spatial
// diffusion. It considers stochastic acceleration by turbulence at each cell, where 
// the steady transport equation is an ordinary Fokker-Plack diff. equation in the energy.
// It is solved by finite differences (CC70, PP96).
// Problems: It does not conserve the number of particles.
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
							//return pow(rateAccretion+rateWind+rateDiff,-1);
							return pow(rateAccretion,-1);
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

void distributionFokkerPlanckSpatialDiffusion(Particle& particle, State& st)
// This routine solves the transport equation by considering separated one zones at each
// energy. The escape time includes the time that takes to the particles to "cross" the energy-cell
// by losing energy, and the injection includes the flux of particles coming from the immediately 
// higher-energy cell. It is solved backward (from the highest energy to the lowest one) and so 
// it can only handle systematic energy losses (cooling) and not acceleration.
// At each cell, the steady transport equation is an ordinary diff. equation
// in the radial variable, it considers spatial advection and diffusion and it is solved by 
// finite differences (CC70, PP96).
{
	double paso_E = pow(particle.emax()/particle.emin(),1.0/(nE-1));
	double restEnergy = particle.mass * cLight2;
	double factor = 1.0;
	particle.ps.iterate([&](const SpaceIterator& iE) {
		
		SpaceCoord jE = {nE-1-iE.coord[DIM_E],0,0};
		double E = particle.ps[DIM_E][jE[DIM_E]];
		double g = E / restEnergy;
		
		// We define a new mesh of points
		size_t M = 120;
		double r_min = st.denf_e.ps[DIM_R][0];// / schwRadius;
		double r_max = st.denf_e.ps[DIM_R][nR-1];// / schwRadius;
		
		double paso_r = pow(r_max/r_min,1.0/M);
		Vector r_au(M+3,r_min/paso_r);
		Vector delta_r_m_au(M+2,0.0);
		Vector delta_r(M+1,0.0);
		
		for (size_t jr=1;jr<M+3;jr++)	r_au[jr] = r_au[jr-1]*paso_r;
		for (size_t jr=0;jr<M+2;jr++)	delta_r_m_au[jr] = r_au[jr+1]-r_au[jr];
		for (size_t jr=0;jr<M+1;jr++)	delta_r[jr] = 0.5*(r_au[jr+2]-r_au[jr]);
		
		fun1 Bfun = [&particle,&st,&jE,g] (double rTrue) 
						{
							double height = height_fun(rTrue);
							double B = magneticField(rTrue);
							double vR = radialVel(rTrue);
							double k_rr = diffCoeff_r(g,particle,height,B);
							return - vR * ( 2.0*k_rr/(rTrue*vR) + 1.0 );
						};

		fun1 Cfun = [&particle,&jE,g] (double rTrue)
						{
							double height = height_fun(rTrue);
							double B = magneticField(rTrue);
							double vR = radialVel(rTrue);
							double k_rr = diffCoeff_r(g,particle,height,B);
							return k_rr;
						};

		fun1 Tfun = [&st,&particle,&jE,g,E,paso_E,restEnergy] (double rTrue) 
						{
							size_t jjR=0;
							particle.ps.iterate([&](const SpaceIterator& iR) {
								if (iR.val(DIM_R) < rTrue) jjR++;
							},{0,-1,0});
							SpaceCoord jEaux = {jE[DIM_E],jjR,0};
							double height = height_fun(rTrue);
							double B = magneticField(rTrue);
							double vR = radialVel(rTrue);
							double rateDiff = diffCoeff_r(g,particle,height,B) / (height*height);
							double rateWind = 2.0*s*abs(vR)/rTrue;
							double loss = abs(losses(E,particle,st,jEaux));
							double rateCool = loss / (E*(sqrt(paso_E)-1.0/sqrt(paso_E)));
							double rateAccretion = (rTrue > st.denf_e.ps[DIM_R][0] && rTrue < 2*schwRadius) ? 
														pow(accretionTime(rTrue),-1) : 0.0;
							return pow(rateDiff+rateCool,-1);
						};
		
		fun1 Qfun = [&st,&particle,&jE,E,paso_E,restEnergy,factor] (double rTrue)
					{
						size_t jjR=0;
						particle.ps.iterate([&](const SpaceIterator& iR) {
							if (iR.val(DIM_R) < rTrue) jjR++;
						},{0,-1,0});
						SpaceCoord jEaux = {jE[DIM_E],jjR,0};
						double height = height_fun(rTrue);
						double Qplus = 0.0;
						if (jE[DIM_E] < nE-1) {
							SpaceCoord jEplus = {jE[DIM_E]+1,jjR,0};
							double Nplus = (rTrue > particle.ps[DIM_R].first() &&
											rTrue < particle.ps[DIM_R].last()) ?
								particle.distribution.interpolate({{DIM_R,rTrue}},&jEplus)*rTrue*height : 0.0;
							double loss = losses(st.denf_e.ps[DIM_E][jE[DIM_E]],particle,st,jEplus);
							Qplus = Nplus * abs(loss) / (E*(sqrt(paso_E)-1.0/sqrt(paso_E)));
						}
						double vR = radialVel(rTrue);
						double Qlocal = (rTrue > particle.ps[DIM_R].first() && 
								rTrue < particle.ps[DIM_R].last()) ?
								particle.injection.interpolate({{DIM_R,rTrue}},&jEaux)*rTrue*height: 0.0;
						//return factor * (rTrue > 100*schwRadius ? (Qlocal+Qplus) : Qplus);
						return factor * (Qlocal + Qplus);
					};
		
		Vector d(M+1,0.0);
		Vector a(M+1,0.0), b(M+1,0.0), c(M+1,0.0);
		for (size_t m=0;m<M+1;m++) {
			size_t m2 = m+1;
			double Bm_minusHalf = 0.5 * (Bfun(r_au[m2]) + Bfun(r_au[m2+1]));
			double Bm_plusHalf = 0.5 * (Bfun(r_au[m2+1]) + Bfun(r_au[m2+2]));
			double Cm_minusHalf = 0.5 * (Cfun(r_au[m2]) + Cfun(r_au[m2+1]));
			double Cm_plusHalf = 0.5 * (Cfun(r_au[m2+1]) + Cfun(r_au[m2+2]));
			double Qm = Qfun(r_au[m2+1]);
			double Tm = Tfun(r_au[m2+1]);
			double wm_minusHalf = Bm_minusHalf/Cm_minusHalf * delta_r_m_au[m2];
			double wm_plusHalf = Bm_plusHalf/Cm_plusHalf * delta_r_m_au[m2+1];
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
			
			a[m] = - Cm_minusHalf * Wminus_m_minusHalf / ( delta_r[m2] * delta_r_m_au[m2] );
			if (m <= M-2)
				c[m] = - Cm_plusHalf * Wplus_m_plusHalf / ( delta_r[m2] * delta_r_m_au[m2+1] );
			b[m] = 1.0/delta_r[m2] * ( Cm_minusHalf * Wplus_m_minusHalf / delta_r_m_au[m2] +
						Cm_plusHalf * Wminus_m_plusHalf / delta_r_m_au[m2+1] ) + 1.0/Tm;
			d[m] = Qm;
		}
		TriDiagSys(a,b,c,d,M-1);
		Vector dNew(M+1,0.0);
		for (size_t j=1;j<M+1;j++) dNew[j] = d[j-1];
		
		size_t m = 1;
		particle.ps.iterate([&](const SpaceIterator &jER) {
			double rTrue = jER.val(DIM_R);
			double height = height_fun(rTrue);
			double vR = radialVel(rTrue);
			double logrTrue = log10(rTrue);// / schwRadius);
			while (log10(r_au[m]) < logrTrue)
				m++;
			//double slope = (m > 1) ? safeLog10(d[m]/d[m-1])/safeLog10(r_au[m]/r_au[m-1]) : 0.0;
			double slope = safeLog10(dNew[m]/dNew[m-1])/safeLog10(r_au[m]/r_au[m-1]);
			double dist = dNew[m-1] * pow(rTrue / r_au[m-1], slope);
			particle.distribution.set(jER,dist/(rTrue*height*factor));//*schwRadius*timeMeasure*restEnergy));
		},{jE[DIM_E],-1,0});
		
		/*
		Vector d(M+1,0.0);
		Vector a(M+1,0.0), b(M+1,0.0), c(M+1,0.0);
		for (size_t m=0;m<M+1;m++) {
			double Bm_minusHalf = 0.5 * (Bfun(r_au[m]) + Bfun(r_au[m+1]));
			double Bm_plusHalf = 0.5 * (Bfun(r_au[m+1]) + Bfun(r_au[m+2]));
			double Cm_minusHalf = 0.5 * (Cfun(r_au[m]) + Cfun(r_au[m+1]));
			double Cm_plusHalf = 0.5 * (Cfun(r_au[m+1]) + Cfun(r_au[m+2]));
			double Qm = Qfun(r_au[m+1]);
			double Tm = Tfun(r_au[m+1]);
			double wm_minusHalf = Bm_minusHalf/Cm_minusHalf * delta_r_m_au[m];
			double wm_plusHalf = Bm_plusHalf/Cm_plusHalf * delta_r_m_au[m+1];
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
				a[m] = - Cm_minusHalf * Wminus_m_minusHalf / ( delta_r[m] * delta_r_m_au[m] );
			if (m <= M-1)
				c[m] = - Cm_plusHalf * Wplus_m_plusHalf / ( delta_r[m] * delta_r_m_au[m+1] );
			b[m] = 1.0/delta_r[m] * ( Cm_minusHalf * Wplus_m_minusHalf / delta_r_m_au[m] +
						Cm_plusHalf * Wminus_m_plusHalf / delta_r_m_au[m+1] ) + 1.0/Tm;
			d[m] = Qm;
		}
		TriDiagSys(a,b,c,d,M);
		
		
		size_t m = 1;
		particle.ps.iterate([&](const SpaceIterator &jER) {
			double rTrue = jER.val(DIM_R);
			double height = height_fun(rTrue);
			double vR = radialVel(rTrue);
			double logrTrue = log10(rTrue);// / schwRadius);
			while (log10(r_au[m]) < logrTrue)
				m++;
			//double slope = (m > 1) ? safeLog10(d[m]/d[m-1])/safeLog10(r_au[m]/r_au[m-1]) : 0.0;
			double slope = safeLog10(d[m]/d[m-1])/safeLog10(r_au[m]/r_au[m-1]);
			double dist = d[m-1] * pow(rTrue / r_au[m-1], slope);
			particle.distribution.set(jER,dist/(rTrue*height*factor));//*schwRadius*timeMeasure*restEnergy));
		},{jE[DIM_E],-1,0});
		 */ 
	},{-1,0,0});

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







/*
void distributionFokkerPlanckER(Particle& particle, State& st)
// This routine solves the transport equation in two dimensions: radius and energy.
// The time-dependence in the usual Fokker-Planck equation is replaced via changes of variables
// by a radial dependence (This is only possible is spatial diffusion is neglected, as it is the case).
// It evolves the distribution in radius considering stochastic acceleration by turbulence; problem 
// that it is solved by finite differences (CC70, PP96).
// Problems: not always work well. Sometimes unphysical peaks appear.
{
	double factor = 1.0;
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
						
						double Dg = diffCoeff_g(g,particle,height,B,rho);
						return - dgdt/vR - 2.0*Dg/vR / g + 0.5*g/r;
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
						//double kRR = diffCoeff_r(g,particle,height,B);
						// vR = vR + 2.0*kRR / rad;
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
						
						//double kRR = diffCoeff_r(g,particle,height,B);
						// vR = vR + 2.0*kRR / rad;
						return 1e99*vR;
						//return pow(rateWind + rateDiff,-1) * vR;
					};
	
	fun2 Qfun = [g_inj,sigma,&particle,factor] (double g, double r)
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
						//double kRR = diffCoeff_r(g,particle,height,B);
						// vR = vR + 2.0*kRR / rad;
						double vA = magneticField(rad)/sqrt(4.0*pi*massDensityADAF(rad));
						double height = height_fun(rad);
						double Q0 = factor*Ainjection * dens * cLight/schwRadius / g_inj;
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
		for (size_t m=0;m<M+1;m++) { 
			double vR = radialVel(radius[jr]);
			//double kRR = diffCoeff_r(g,particle,height,B);
			// vR = vR + 2.0*kRR / rad;
			dist[jr][m] = d[m]/(radius[jr]*height_fun(radius[jr])*vR*factor);
		}
	}
	
	size_t jr(Nrad);
	particle.ps.iterate([&](const SpaceIterator& jR) {
		double logr = log10(jR.val(DIM_R));
		while (log10(radius[jr]) < logr) jr--;
		if (jr >= Nrad-1) jr = Nrad-2;
		size_t m(1);
		particle.ps.iterate([&](const SpaceIterator &jRE) {
			double logg = log10(jRE.val(DIM_E)/(particle.mass*cLight2));
			while (log10(g_au[m]) < logg) m++;
			
			double slope1 = safeLog10(dist[jr][m] / dist[jr][m-1]) / safeLog10(g_au[m]/g_au[m-1]);
			double dist_1 = slope1 * (logg-log10(g_au[m-1])) + safeLog10(dist[jr][m-1]);
			double slope2 = safeLog10(dist[jr+1][m]/dist[jr+1][m-1]) / safeLog10(g_au[m]/g_au[m-1]);
			double dist_2 = slope2 * (logg-log10(g_au[m-1])) + safeLog10(dist[jr+1][m-1]);
			double slope = safeLog10(dist_1/dist_2) / safeLog10(radius[jr]/radius[jr+1]);
			double distr = slope * (logr-log10(radius[jr+1])) + dist_2;

			if (m < M+1) 
				distr = pow(10,distr);
			else
				distr = 0.0;
			
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
// This routine solves the transport equation in two dimensions: radius and energy.
// The time-dependence in the usual Fokker-Planck equation is replaced via changes of variables
// by a radial dependence (This is only possible is spatial diffusion is neglected, as it is the case).
// It evolves the distribution in radius without considering acceleration; the problem 
// is solved by finite differences (CC70, PP96).
// Problems: not always work well. Sometimes unphysical peaks appear.
{
	
	// We define a new mesh of points
	size_t M = 199;
	double factor = 1e10;
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
						double rateWind = 2.0*s*abs(vR)/rad;
						double B = magneticField(rad);
						double height = height_fun(rad);
						double rateDiff = 1.0/(diffusionTimeTurbulence(E,height,particle,B));
						rateDiff = 0.0;
						return pow(rateDiff + rateWind,-1) * vR;
					};
	
	fun2 Qfun = [&particle,&factor] (double g, double r)
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
						double height = height_fun(rad);
						return (1.0/factor)*particle.injection.interpolate({{DIM_E,E},{DIM_R,rad}},&jR)*(particle.mass*cLight2)*rad*height;
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
			
			particle.distribution.set(jRE,factor*distr/(particle.mass*cLight2));
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





void distributionFokkerPlanckComplete(Particle& particle, State& st)
// This routine solves the complete time-dependent transport equation considering turbulence
// acceleration and spatial diffusion.
{
	double factor = 1.0e-20;
	
	size_t J = 50;
	double r_min = sqrt(st.denf_e.ps[DIM_R][0] / schwRadius);
	double r_max = (st.denf_e.ps[DIM_R][nR-1]/schwRadius)*paso_r;
	r_max = 1.0e2*paso_r;
	double paso_r_local = pow(r_max/r_min,1.0/J);
	
	// We define a mesh of points for the radial coordinate
	
	Vector r_local(J+1,r_min);
	Vector r_local_extended(J+3,r_min/paso_r_local);
	Vector delta_r_j_plushalf(J+1,0.0);
	Vector delta_r_j_minushalf(J+1,0.0);
	Vector delta_r(J+1,0.0);
		
	for (size_t j=1;j<J+1;j++)	r_local[j] = r_local[j-1]*paso_r_local;
	for (size_t j=1;j<J+3;j++)	r_local_extended[j] = r_local_extended[j-1]*paso_r_local;
	for (size_t j=0;j<J+1;j++)	delta_r_j_plushalf[j] = r_local[j]*(paso_r_local-1.0);
	for (size_t j=0;j<J+1;j++)	delta_r_j_minushalf[j] = r_local[j]*(1.0-1.0/paso_r_local);
	for (size_t j=0;j<J+1;j++)	delta_r[j] = 0.5*r_local[j]*(paso_r_local-1.0/paso_r_local);
	
	// We define a mesh of points for the Lorentz factor coordinate
	size_t M = 80;
	double g_min = 0.9 * (particle.emin() / (particle.mass*cLight2));
	double g_max = 1.1 * (particle.emax() / (particle.mass*cLight2));
	
	double paso_g_local = pow(g_max/g_min,1.0/M);
	Vector g_local(M+1,g_min);
	Vector g_local_extended(M+3,g_min/paso_g_local);
	Vector delta_g_m_plushalf(M+1,0.0);
	Vector delta_g_m_minushalf(M+1,0.0);
	Vector delta_g(M+1,0.0);
		
	for (size_t m=1;m<M+1;m++)	g_local[m] = g_local[m-1]*paso_g_local;
	for (size_t m=1;m<M+3;m++)	g_local_extended[m] = g_local_extended[m-1]*paso_g_local;
	for (size_t m=0;m<M+1;m++)	delta_g_m_plushalf[m] = g_local[m]*(paso_g_local-1.0);
	for (size_t m=0;m<M+1;m++)	delta_g_m_minushalf[m] = g_local[m]*(1.0-1.0/paso_g_local);
	for (size_t m=0;m<M+1;m++)	delta_g[m] = 0.5*g_local[m]*(paso_g_local-1.0/paso_g_local);
	
	// We define a vector for the time evolution
	
	double naturalTimescale = schwRadius/cLight;
	for (size_t j=2;j<J-1;j++) {
		double rTrue = r_local[j]*schwRadius;
		double deltar = delta_r[j]*schwRadius;
		double magf = magneticField(rTrue);
		double rho = massDensityADAF(rTrue);
		double height = height_fun(rTrue);
		double vR = radialVel(rTrue);
		double advTime = deltar / abs(vR);
		for (size_t m=2;m<M-1;m++) {
			double localAccTime = P2(delta_g[m]) / diffCoeff_g(g_local[m],particle,height,magf,rho);
			double localDiffTime = P2(deltar) / diffCoeff_r(g_local[j],particle,height,magf);
			naturalTimescale = min(naturalTimescale,min(localDiffTime,min(advTime,localAccTime)));
		}
	}
	
	cout << "Acceleration timescale = " << naturalTimescale << endl;
	cout << "Natural timescale = " << schwRadius/cLight << endl;
	
	naturalTimescale = schwRadius/cLight;
	
	size_t Ntime = 5e2;
	size_t NdeltaTime = 1e2;
	
	Vector time(Ntime+1,1.0);
	double dt = naturalTimescale;
	//double paso_time = pow(maxTime/naturalTimescale,1.0/(Ntime-1));
	//for (size_t t=1;t<Ntime+1;t++) time[t] = time[t-1]*paso_time;
	for (size_t t=1;t<Ntime+1;t++) time[t] = time[t-1] + dt;
	Vector deltat(Ntime,1.0);
	for (size_t t=0;t<Ntime;t++) deltat[t] = (time[t+1]-time[t])/NdeltaTime;
	
	fun2 Afun = [&particle,naturalTimescale] (double r, double g) 
				{
					double rTrue = r * schwRadius;
					double height = height_fun(rTrue);
					double B = magneticField(rTrue);
					double vR = radialVel(rTrue);
					double k_rr = diffCoeff_r(g,particle,height,B);
					return - (vR/schwRadius) * naturalTimescale * ( 2.0*k_rr/(rTrue*vR) + 1.0 );
				};

	fun2 Bfun = [&particle,naturalTimescale] (double r, double g)
				{
					double rTrue = r * schwRadius;
					double height = height_fun(rTrue);
					double B = magneticField(rTrue);
					double k_rr = diffCoeff_r(g,particle,height,B);
					return k_rr / (schwRadius*schwRadius) * naturalTimescale;
				};

	fun2 Cfun = [&particle,&st,paso_r_local,naturalTimescale] (double r, double g) 
				{
					double rTrue = r * schwRadius;
					double E = g*particle.mass*cLight2;
					
					size_t jjR=0;
					if (rTrue > particle.ps[DIM_R][0]) {
						if (rTrue < particle.ps[DIM_R][nR-1]) {
							while (particle.ps[DIM_R][jjR] < rTrue) jjR++;
						} else
							jjR = nR-1;
					}
					
					double dgdt = 0.0;
					if (jjR > 0 && jjR < nR-1) {
						SpaceCoord jR1 = {0,jjR-1,0};
						SpaceCoord jR2 = {0,jjR,0};
						double r1 = st.denf_e.ps[DIM_R][jjR-1];
						double r2 = st.denf_e.ps[DIM_R][jjR];
						double dEdt_1 = losses(E,particle,st,jR1);
						double dEdt_2 = losses(E,particle,st,jR2);
						double slope = safeLog10(dEdt_2/dEdt_1)/safeLog10(r2/r1);
						dgdt = - dEdt_1 * pow(rTrue/r1, slope) / (particle.mass*cLight2);
						
					} else if (jjR == 0) {
						SpaceCoord jR = {0,0,0};
						dgdt = - losses(E,particle,st,jR) / (particle.mass*cLight2);
					} else if (jjR == nR-1) {
						SpaceCoord jR = {0,nR-1,0};
						dgdt = - losses(E,particle,st,jR) / (particle.mass*cLight2);
					}
					double height = height_fun(rTrue);
					double B = magneticField(rTrue);
					double rho = massDensityADAF(rTrue);
					double Dg = diffCoeff_g(g,particle,height,B,rho);
					double vR = radialVel(rTrue);
					double dvRdr = (radialVel(rTrue*paso_r_local) - radialVel(rTrue/paso_r_local)) / 
									(rTrue*(paso_r_local-1.0/paso_r_local));
					
					dvRdr = (dvRdr > 0.0) ? dvRdr : 0.0;
					return (-2.0*Dg/g - dgdt + 1.0/3.0 * (2.0*vR/rTrue + dvRdr) * g ) * naturalTimescale;
				};
	
	fun2 Dfun = [&particle,naturalTimescale] (double r, double g)
				{
					double rTrue = r * schwRadius;
					double height = height_fun(rTrue);
					double B = magneticField(rTrue);
					double rho = massDensityADAF(rTrue);
					double Dg = diffCoeff_g(g,particle,height,B,rho);
					return naturalTimescale * Dg;
				};
	
	fun2 Tfun = [&particle,naturalTimescale] (double r, double g) 
				{
					double rTrue = r * schwRadius;
					double height = height_fun(rTrue);
					double B = magneticField(rTrue);
					double vR = radialVel(rTrue);
					double rateDiff = diffCoeff_r(g,particle,height,B) / (height*height);
					double rateWinds = 2.0*s*abs(vR)/rTrue;
					return 1e99;
					return pow(rateDiff+rateWinds,-1) / naturalTimescale;
				};
		
	fun2 Qfun = [&st,&particle,paso_g_local,factor,naturalTimescale,r_local,J] (double r, double g)
				{
					double rTrue = r * schwRadius;
					double E = g * particle.mass*cLight2;
					double g_inj = 2.0;
					double sigma = g_inj * (paso_g_local-1.0)/2.0;
					double Qlocal = (rTrue > particle.ps[DIM_R][0] && 
							rTrue < particle.ps[DIM_R][nR-1]) ?
							rTrue*rTrue * Ainjection * ionDensity(rTrue) : 0.0;
					double result = Qlocal * gsl_ran_gaussian_pdf(g-g_inj,sigma);
					double smoother = 0.25 * (1.0 + gsl_sf_erf((r_local[J-5]-r)/(r_local[J-5]/10.0))) 
											* (1.0 + gsl_sf_erf(r-r_local[0]));
					return (r>80) ? result * naturalTimescale * factor * smoother : 0.0;
				};
	
	Matrix Bj_minusHalf_m;				matrixInit(Bj_minusHalf_m,M+1,J+1,0.0);
	Matrix Bj_plusHalf_m;				matrixInit(Bj_plusHalf_m,M+1,J+1,0.0);
	Matrix Djm_minusHalf;				matrixInit(Djm_minusHalf,J+1,M+1,0.0);
	Matrix Djm_plusHalf;				matrixInit(Djm_plusHalf,J+1,M+1,0.0);
	Matrix Wminus_j_minusHalf;			matrixInit(Wminus_j_minusHalf,M+1,J+1,0.0);
	Matrix Wminus_j_plusHalf;			matrixInit(Wminus_j_plusHalf,M+1,J+1,0.0);
	Matrix Wplus_j_minusHalf;			matrixInit(Wplus_j_minusHalf,M+1,J+1,0.0);
	Matrix Wplus_j_plusHalf;			matrixInit(Wplus_j_plusHalf,M+1,J+1,0.0);
	Matrix Vminus_m_minusHalf;			matrixInit(Vminus_m_minusHalf,J+1,M+1,0.0);
	Matrix Vminus_m_plusHalf;			matrixInit(Vminus_m_plusHalf,J+1,M+1,0.0);
	Matrix Vplus_m_minusHalf;			matrixInit(Vplus_m_minusHalf,J+1,M+1,0.0);
	Matrix Vplus_m_plusHalf;			matrixInit(Vplus_m_plusHalf,J+1,M+1,0.0);
	Matrix Tjm;							matrixInit(Tjm,J+1,M+1,0.0);
	Matrix Qjm;							matrixInit(Qjm,J+1,M+1,0.0);
	
	cout << "TRANSPORT EQUATION: Starting matrix element calculation" << endl << endl;
	
	#pragma omp parallel for
	for (size_t j=0;j<J+1;j++) {
		for (size_t m=0;m<M+1;m++) {
			size_t l = j*(M+1)+m;
			double Ajm = Afun(r_local[j],g_local[m]);
			double Aj_minusHalf_m = 0.5 * (Afun(r_local_extended[j],g_local[m]) + Ajm);
			double Aj_plusHalf_m = 0.5 * (Afun(r_local_extended[j+2],g_local[m]) + Ajm);
					
			double Bjm = Bfun(r_local[j],g_local[m]);
			Bj_minusHalf_m[m][j] = 0.5 * (Bfun(r_local_extended[j],g_local[m]) + Bjm);
			Bj_plusHalf_m[m][j] = 0.5 * (Bfun(r_local_extended[j+2],g_local[m]) + Bjm);
			
			double Cjm = Cfun(r_local[j],g_local[m]);
			double Cjm_minusHalf = 0.5 * (Cfun(r_local[j],g_local_extended[m]) + Cjm);
			double Cjm_plusHalf = 0.5 * (Cfun(r_local[j],g_local_extended[m+2]) + Cjm);
			
			double Djm = Dfun(r_local[j],g_local[m]);
			Djm_minusHalf[j][m] = 0.5 * (Dfun(r_local[j],g_local_extended[m]) + Djm);
			Djm_plusHalf[j][m] = 0.5 * (Dfun(r_local[j],g_local_extended[m+2]) + Djm);
			
			Qjm[j][m] = Qfun(r_local[j],g_local[m]);
			Tjm[j][m] = Tfun(r_local[j],g_local[m]);
					
			double wj_minusHalf = Aj_minusHalf_m/Bj_minusHalf_m[m][j] * delta_r_j_minushalf[j];
			double wj_plusHalf = Aj_plusHalf_m/Bj_plusHalf_m[m][j] * delta_r_j_plushalf[j];
			double Wj_minusHalf(0.0), Wj_plusHalf(0.0);
			
			if ( abs(wj_minusHalf) < 1.0e-3 ) {
				Wj_minusHalf = pow(1.0+gsl_pow_2(wj_minusHalf)/24+gsl_pow_4(wj_minusHalf)/1920,-1);
				Wplus_j_minusHalf[m][j] = Wj_minusHalf * exp(0.5*wj_minusHalf);
				Wminus_j_minusHalf[m][j] = Wj_minusHalf * exp(-0.5*wj_minusHalf);
			} else {
				Wj_minusHalf = abs(wj_minusHalf)*exp(-0.5*abs(wj_minusHalf)) /
								(1.0-exp(-abs(wj_minusHalf)));
				if (wj_minusHalf > 0.0) {
					Wplus_j_minusHalf[m][j] = abs(wj_minusHalf) / (1.0-exp(-abs(wj_minusHalf)));
					Wminus_j_minusHalf[m][j] = Wj_minusHalf * exp(-0.5*wj_minusHalf);
				} else {
					Wplus_j_minusHalf[m][j] = Wj_minusHalf * exp(0.5*wj_minusHalf);
					Wminus_j_minusHalf[m][j] = abs(wj_minusHalf) / (1.0-exp(-abs(wj_minusHalf)));
				}
			}
			if ( abs(wj_plusHalf) < 1.0e-3 ) {
				Wj_plusHalf = pow(1.0+gsl_pow_2(wj_plusHalf)/24+gsl_pow_4(wj_plusHalf)/1920,-1);
				Wplus_j_plusHalf[m][j] = Wj_plusHalf * exp(0.5*wj_plusHalf);
				Wminus_j_plusHalf[m][j] = Wj_plusHalf * exp(-0.5*wj_plusHalf);
			} else {
				Wj_plusHalf = abs(wj_plusHalf)*exp(-0.5*abs(wj_plusHalf)) /
								(1.0-exp(-abs(wj_plusHalf)));
				if (wj_plusHalf > 0.0) {
					Wplus_j_plusHalf[m][j] = abs(wj_plusHalf) / (1.0-exp(-abs(wj_plusHalf)));
					Wminus_j_plusHalf[m][j] = Wj_plusHalf * exp(-0.5*wj_plusHalf);
				} else {
					Wplus_j_plusHalf[m][j] = Wj_plusHalf * exp(0.5*wj_plusHalf);
					Wminus_j_plusHalf[m][j] = abs(wj_plusHalf) / (1.0-exp(-abs(wj_plusHalf)));
				}
			}
			
			double vm_minusHalf = Cjm_minusHalf/Djm_minusHalf[j][m] * delta_g_m_minushalf[m];
			double vm_plusHalf = Cjm_plusHalf/Djm_plusHalf[j][m] * delta_g_m_plushalf[m];
			double Vm_minusHalf(0.0), Vm_plusHalf(0.0);
			
			if ( abs(vm_minusHalf) < 1.0e-3 ) {
				Vm_minusHalf = pow(1.0+gsl_pow_2(vm_minusHalf)/24.0+gsl_pow_4(vm_minusHalf)/1920.0,-1);
				Vplus_m_minusHalf[j][m] = Vm_minusHalf * exp(0.5*vm_minusHalf);
				Vminus_m_minusHalf[j][m] = Vm_minusHalf * exp(-0.5*vm_minusHalf);
			} else {
				Vm_minusHalf = abs(vm_minusHalf)*exp(-0.5*abs(vm_minusHalf)) /
								(1.0-exp(-abs(vm_minusHalf)));
				if (vm_minusHalf > 0.0) {
					Vplus_m_minusHalf[j][m] = abs(vm_minusHalf) / (1.0-exp(-abs(vm_minusHalf)));
					Vminus_m_minusHalf[j][m] = Vm_minusHalf * exp(-0.5*vm_minusHalf);
				} else {
					Vplus_m_minusHalf[j][m] = Vm_minusHalf * exp(0.5*vm_minusHalf);
					Vminus_m_minusHalf[j][m] = abs(vm_minusHalf) / (1.0-exp(-abs(vm_minusHalf)));
				}
			}
			if ( abs(vm_plusHalf) < 1.0e-3 ) {
				Vm_plusHalf = pow(1.0+gsl_pow_2(vm_plusHalf)/24+gsl_pow_4(vm_plusHalf)/1920,-1);
				Vplus_m_plusHalf[j][m] = Vm_plusHalf * exp(0.5*vm_plusHalf);
				Vminus_m_plusHalf[j][m] = Vm_plusHalf * exp(-0.5*vm_plusHalf);
			} else {
				Vm_plusHalf = abs(vm_plusHalf)*exp(-0.5*abs(vm_plusHalf)) /
								(1.0-exp(-abs(vm_plusHalf)));
				if (vm_plusHalf > 0.0) {
					Vplus_m_plusHalf[j][m] = abs(vm_plusHalf) / (1.0-exp(-abs(vm_plusHalf)));
					Vminus_m_plusHalf[j][m] = Vm_plusHalf * exp(-0.5*vm_plusHalf);
				} else {
					Vplus_m_plusHalf[j][m] = Vm_plusHalf * exp(0.5*vm_plusHalf);
					Vminus_m_plusHalf[j][m] = abs(vm_plusHalf) / (1.0-exp(-abs(vm_plusHalf)));
				}
			}
		}
	}
	
	ofstream fileMatrixDistEne, fileMatrixDistPos, fileMatrixDistQT;
	fileMatrixDistEne.open("matrixCoeffTransportEq_Ene.txt");
	fileMatrixDistPos.open("matrixCoeffTransportEq_Pos.txt");
	fileMatrixDistQT.open("matrixCoeffTransportEq_QT.txt");
	fileMatrixDistEne 	<< "j \t m \t Dj- \t Dj+ \t V-- \t V-+ \t V+- \t V++" << endl;
	fileMatrixDistPos 	<< "m \t j \t Bj- \t Bj+ \t W-- \t W-+ \t W+- \t W++" << endl;
	fileMatrixDistQT 	<< "j \t m \t T \t Q \t" << endl;
	
	for (size_t j=0;j<J+1;j++) {
		for (size_t m=0;m<M+1;m++) {
			fileMatrixDistEne 	<< j << "\t" << m << "\t"
								<< Djm_minusHalf[j][m] << "\t" << Djm_plusHalf[j][m] << "\t"
								<< Vminus_m_minusHalf[j][m] << "\t" << Vminus_m_plusHalf[j][m] << "\t"
								<< Vplus_m_minusHalf[j][m] << "\t" << Vplus_m_plusHalf[j][m] << endl;
			fileMatrixDistQT	<< j << "\t" << m << "\t"
								<< Tjm[j][m] << "\t" << Qjm[j][m] << endl;
		}
	}
	for (size_t m=0;m<M+1;m++) {
		for (size_t j=0;j<J+1;j++) {
			fileMatrixDistPos	<< m << "\t" << j << "\t"
								<< Bj_minusHalf_m[m][j] << "\t" << Bj_plusHalf_m[m][j] << "\t"
								<< Wminus_j_minusHalf[m][j] << "\t" << Wminus_j_plusHalf[m][j] << "\t"
								<< Wplus_j_minusHalf[m][j] << "\t" << Wplus_j_plusHalf[m][j] << endl;
		}
	}
 
	fileMatrixDistEne.close();
	fileMatrixDistPos.close();
	fileMatrixDistQT.close();
	
	cout << "TRANSPORT EQUATION: Finished matrix calculation" << endl << endl;
	
	Vector rr((J+1)*(M+1),0.0);
	Matrix rrOld;
	matrixInit(rrOld,J+1,M+1,0.0);
	Matrix dist;
	
	cout << "TRANSPORT EQUATION: Starting time evolution" << endl << endl;
	
	for (size_t t=0;t<Ntime;t++) {
		cout << "t = " << t << endl;
		for (size_t jt=0;jt<NdeltaTime;jt++) {
			
			Vector a((J+1)*(M+1),0.0);
			Vector b((J+1)*(M+1),0.0);
			Vector c((J+1)*(M+1),0.0);
			Vector d((J+1)*(M+1),0.0);
			Vector e((J+1)*(M+1),0.0);
			Vector f((J+1)*(M+1),0.0);
			
			// BLOCK 1: Advances in half a time step the energy evolution operator.
			
			#pragma omp parallel for
			for (size_t j=0;j<J+1;j++) {
				
				for (size_t m=0;m<M+1;m++) {
					size_t l = j*(M+1)+m;
				
					if (j > 0 && j < J) {
					
						if (m > 0)
							a[l] = - (0.5*deltat[t])/delta_g[m] * Djm_minusHalf[j][m] * Vminus_m_minusHalf[j][m] / delta_g_m_minushalf[m];
						if (m < M)
							c[l] = - (0.5*deltat[t])/delta_g[m] * Djm_plusHalf[j][m] * Vplus_m_plusHalf[j][m] / delta_g_m_plushalf[m];
						b[l] = 1.0 + (0.5*deltat[t])/delta_g[m] * ( Djm_minusHalf[j][m] * Vplus_m_minusHalf[j][m] / delta_g_m_minushalf[m] +
									Djm_plusHalf[j][m] * Vminus_m_plusHalf[j][m] / delta_g_m_plushalf[m] ) + 0.5*deltat[t]/Tjm[j][m];
						
						d[l] = (0.5*deltat[t])/delta_r[j] * Bj_minusHalf_m[m][j] * Wminus_j_minusHalf[m][j] / delta_r_j_minushalf[j];
						e[l] = 1.0 - 0.5*(deltat[t])/delta_r[j] * ( Bj_minusHalf_m[m][j] * Wplus_j_minusHalf[m][j] / delta_r_j_minushalf[j] +
									Bj_plusHalf_m[m][j] * Wminus_j_plusHalf[m][j] / delta_r_j_plushalf[j] );
						f[l] = (j < J) ? 
							(0.5*deltat[t])/delta_r[j] * Bj_plusHalf_m[m][j] * Wplus_j_plusHalf[m][j] / delta_r_j_plushalf[j] : 0.0;
					} else {
						b[l] = 1.0;
					}
					rr[l] = rrOld[j][m] + (j > 0 && j < J ? Qjm[j][m] * deltat[t] : 0.0);
				}
			}
			TriDiagSys(a,b,c,rr,(J+1)*(M+1));
			
			#pragma omp parallel for
			for (size_t j=0;j<J+1;j++) {
				for (size_t m=0;m<M+1;m++) {
					size_t l = j*(M+1)+m;
					rrOld[j][m] = rr[l];
				}
			}
			
			// BLOCK 2: Advances in a time step the spatial evolution operator.
			
			#pragma omp parallel for
			for (size_t m=0;m<M+1;m++){
				for (size_t j=0;j<J+1;j++) {
					size_t lp = m*(J+1)+j;
					size_t l = j*(M+1)+m;
					a[lp] = - d[l];
					b[lp] = 2.0 - e[l];
					c[lp] = (j < J) ? - f[l] : 0.0;
					
					rr[lp] = rrOld[j][m];
				}
			}
			TriDiagSys(a,b,c,rr,(J+1)*(M+1));
			
			#pragma omp parallel for
			for (size_t j=0;j<J+1;j++) {
				for (size_t m=0;m<M+1;m++) {
					size_t lp = m*(J+1)+j;
					rrOld[j][m] = rr[lp];
				}
			}
			
			// BLOCK 3: Advances in half a time the energy evolution operator.
			fill(a.begin(),a.end(),0.0);
			fill(b.begin(),b.end(),0.0);
			fill(c.begin(),c.end(),0.0);
			
			#pragma omp parallel for
			for (size_t j=0;j<J+1;j++) {
				for (size_t m=0;m<M+1;m++) {
					size_t l = j*(M+1)+m;
					
					if (j > 0 && j < J) {
						
						if (m > 0)
							a[l] = - (0.5*deltat[t])/delta_g[m] * Djm_minusHalf[j][m] * Vminus_m_minusHalf[j][m] / delta_g_m_minushalf[m];
						if (m < M)
							c[l] = - (0.5*deltat[t])/delta_g[m] * Djm_plusHalf[j][m] * Vplus_m_plusHalf[j][m] / delta_g_m_plushalf[m];
						b[l] = 1.0 + (0.5*deltat[t])/delta_g[m] * ( Djm_minusHalf[j][m] * Vplus_m_minusHalf[j][m] / delta_g_m_minushalf[m] +
									Djm_plusHalf[j][m] * Vminus_m_plusHalf[j][m] / delta_g_m_plushalf[m] ) + 0.5*deltat[t]/Tjm[j][m];
					} else {
						b[l] = 1.0;
					}
					rr[l] = rrOld[j][m];
				}
			}
			TriDiagSys(a,b,c,rr,(J+1)*(M+1));
			
			#pragma omp parallel for
			for (size_t j=0;j<J+1;j++) {
				for (size_t m=0;m<M+1;m++) {
					size_t l = j*(M+1)+m;
					rrOld[j][m] = rr[l];
				}
			}
			
		}
		if (t % 100 == 0 || t == Ntime-2) {
			matrixInit(dist,J+1,M+1,0.0);
			ofstream fileDist;
			fileDist.open("dist"+to_string(t)+".txt");
			for (size_t j=0;j<J+1;j++) {
				for (size_t m=0;m<M+1;m++) {
					size_t l = j*(M+1)+m;
					dist[j][m] = rrOld[j][m] / P2(r_local[j]*schwRadius) / factor;
					fileDist << r_local[j] << "\t" << -radialVel(r_local[j]*schwRadius)/cLight
								   << "\t" << g_local[m] << "\t" << dist[j][m] << endl;
				}
			}
			fileDist.close();
		}
	}
	
	// TO IMPROVE THE SOLUTION
	/*
	Ntime = 101;
	double deltat_shorter = 1.6e-3 * naturalTimescale;
	NdeltaTime = (time[1]-time[0]) / deltat_shorter;
	Vector deltat_shorterVec(Ntime,deltat_shorter);
	cout << NdeltaTime << endl;
	
	for (size_t t=0;t<Ntime;t++) {
		cout << "t = " << t << endl;
		for (size_t jt=0;jt<NdeltaTime;jt++) {
			
			Vector a((J+1)*(M+1),0.0);
			Vector b((J+1)*(M+1),0.0);
			Vector c((J+1)*(M+1),0.0);
			Vector d((J+1)*(M+1),0.0);
			Vector e((J+1)*(M+1),0.0);
			Vector f((J+1)*(M+1),0.0);
			
			// BLOCK 1: Advances in half a time with the energy operators implicit and the spatial
			// operator explicit.
			
			#pragma omp parallel for
			for (size_t j=0;j<J+1;j++) {
				
				for (size_t m=0;m<M+1;m++) {
					size_t l = j*(M+1)+m;
				
					if (j > 0 && j < J) {
					
						if (m > 0)
							a[l] = - (0.5*deltat_shorterVec[t])/delta_g[m] * Djm_minusHalf[j][m] * Vminus_m_minusHalf[j][m] / delta_g_m_minushalf[m];
						if (m < M)
							c[l] = - (0.5*deltat_shorterVec[t])/delta_g[m] * Djm_plusHalf[j][m] * Vplus_m_plusHalf[j][m] / delta_g_m_plushalf[m];
						b[l] = 1.0 + (0.5*deltat_shorterVec[t])/delta_g[m] * ( Djm_minusHalf[j][m] * Vplus_m_minusHalf[j][m] / delta_g_m_minushalf[m] +
									Djm_plusHalf[j][m] * Vminus_m_plusHalf[j][m] / delta_g_m_plushalf[m] ) + 0.5*deltat_shorterVec[t]/Tjm[j][m];
						
						d[l] = (0.5*deltat_shorterVec[t])/delta_r[j] * Bj_minusHalf_m[m][j] * Wminus_j_minusHalf[m][j] / delta_r_j_minushalf[j];
						e[l] = 1.0 - 0.5*(deltat_shorterVec[t])/delta_r[j] * ( Bj_minusHalf_m[m][j] * Wplus_j_minusHalf[m][j] / delta_r_j_minushalf[j] +
									Bj_plusHalf_m[m][j] * Wminus_j_plusHalf[m][j] / delta_r_j_plushalf[j] );
						f[l] = (j < J) ? 
							(0.5*deltat_shorterVec[t])/delta_r[j] * Bj_plusHalf_m[m][j] * Wplus_j_plusHalf[m][j] / delta_r_j_plushalf[j] : 0.0;
						
						//rr[l] = (j2 > 0 ? rrOld[j2-1][m] : 0.0) * d[l] + rrOld[j2][m] * e[l] + (j < J ? rrOld[j2+1][m] * f[l] : 0.0) + Qjm[j2][m]*deltat;
					} else {
						b[l] = 1.0;
					}
					rr[l] = (j > 0 ? rrOld[j-1][m] : 0.0) * d[l] + rrOld[j][m] * e[l] 
							+ (j < J ? rrOld[j+1][m] * f[l] : 0.0) + (j > 0 && j < J ? Qjm[j][m] :0.0)*0.5*deltat_shorterVec[t];
				}
			}
			TriDiagSys(a,b,c,rr,(J+1)*(M+1));
			
			#pragma omp parallel for
			for (size_t j=0;j<J+1;j++) {
				for (size_t m=0;m<M+1;m++) {
					size_t l = j*(M+1)+m;
					rrOld[j][m] = rr[l];
				}
			}
		}
		if (t % 1 == 0) {
			matrixInit(dist,J+1,M+1,0.0);
			ofstream fileDist;
			fileDist.open("dist_improved_"+to_string(t)+".txt");
			for (size_t j=0;j<J+1;j++) {
				for (size_t m=0;m<M+1;m++) {
					size_t l = j*(M+1)+m;
					dist[j][m] = rrOld[j][m] / P2(r_local[j]*schwRadius) / factor;
					fileDist << r_local[j] << "\t" << g_local[m] << "\t" << dist[j][m] << endl;
				}
			}
			fileDist.close();
		}
	}
	
	for (size_t j=0;j<J+1;j++) {
		for (size_t m=0;m<M+1;m++) {
			size_t l = j*(M+1)+m;
			rr[l] = rrOld[j][m];
		}
	}
	*/
	
	matrixInit(dist,J+1,M+1,0.0);
	ofstream fileDist;
	fileDist.open("dist_last.txt");
	for (size_t j=0;j<J+1;j++) {
		for (size_t m=0;m<M+1;m++) {
			size_t l = j*(M+1)+m;
			dist[j][m] = rr[l] / P2(r_local[j]*schwRadius) / factor;
			fileDist << r_local[j] << "\t" << -radialVel(r_local[j]*schwRadius)/cLight
								   << "\t" << g_local[m] << "\t" << dist[j][m] << endl;
		}
	}
	fileDist.close();

	cout << endl << "TRANSPORT EQUATION: Finished time evolution" << endl << endl;
	
	particle.ps.iterate([&] (const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		if (r/schwRadius < r_local.back()) {
			size_t j = 1;
			while (r_local[j]*schwRadius < r) j++;
			particle.ps.iterate([&] (const SpaceIterator& iRE) {
				double g = iRE.val(DIM_E) / (particle.mass*cLight2);
				size_t m = 0;
				while (g_local[m] < g) m++;
				
				double N11 = dist[j-1][m-1];
				double N12 = dist[j-1][m];
				double N1 = 0.0;
				if (N11 > 0.0 && N12 > 0.0) {
					double s1 = safeLog10(N12/N11)/safeLog10(g_local[m]/g_local[m-1]);
					N1 = N11 * pow(g / g_local[m-1], s1);
				} else {
					if (N11 > 0.0)
						N1 = - N11 * (g - g_local[m-1])/(g_local[m]-g_local[m-1]) + N11;
					else if (N12 > 0.0)
						N1 = N12 * (g - g_local[m-1])/(g_local[m]-g_local[m-1]);
					else
						N1 = 0.0;
				}
				double N21 = dist[j][m-1];
				double N22 = dist[j][m];
				double N2 = 0.0;
				if (N21 > 0.0 && N22 > 0.0) {
					double s2 = safeLog10(N22/N21)/safeLog10(g_local[m]/g_local[m-1]);
					N2 = N21 * pow(g / g_local[m-1], s2);
				} else {
					if (N21 > 0.0)
						N2 = - N21 * (g - g_local[m-1])/(g_local[m]-g_local[m-1]) + N21;
					else if (N22 > 0.0)
						N2 = N22 * (g - g_local[m-1])/(g_local[m]-g_local[m-1]);
					else
						N2 = 0.0;
				}
				
				double N = 0.0;
				if (N1 > 0.0 && N2 > 0.0) {
					double s = safeLog10(N2/N1)/safeLog10(r_local[j]/r_local[j-1]);
					N = N1 * pow(r / (r_local[j-1]*schwRadius), s);
				} else {
					if (N1 > 0.0)
						N = - N1 * (r/schwRadius - r_local[j-1])/(r_local[j]-r_local[j-1]) + N1;
					else if (N2 > 0.0)
						N = N2 * (r/schwRadius - r_local[j-1])/(r_local[j]-r_local[j-1]);
					else
						N = 0.0;
				}
				particle.distribution.set(iRE, N / (particle.mass*cLight2));
			},{-1,iR.coord[DIM_R],0});
		}
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

void distributionFokkerPlanckComplete2(Particle& particle, State& st)
// This routine solves the transport equation by considering separated one zones at each
// energy. The escape time includes the time that takes to the particles to "cross" the energy-cell
// by losing energy, and the injection includes the flux of particles coming from the immediately 
// higher-energy cell. It is solved backward (from the highest energy to the lowest one) and so 
// it can only handle systematic energy losses (cooling) and not acceleration.
// At each cell, the steady transport equation is an ordinary diff. equation
// in the radial variable, it considers spatial advection and diffusion and it is solved by 
// finite differences (CC70, PP96).
{
	double factor = 1.0e-20;
	
	size_t J = 40;
	double r_min = sqrt(st.denf_e.ps[DIM_R][0] / schwRadius);
	double r_max = (st.denf_e.ps[DIM_R][nR-1]/schwRadius)*paso_r;
	double paso_r_local = pow(r_max/r_min,1.0/J);
	
	// We define a mesh of points for the radial coordinate
	
	Vector r_local(J+1,r_min);
	Vector r_local_extended(J+3,r_min/paso_r_local);
	Vector delta_r_j_plushalf(J+1,0.0);
	Vector delta_r_j_minushalf(J+1,0.0);
	Vector delta_r(J+1,0.0);
		
	for (size_t j=1;j<J+1;j++)	r_local[j] = r_local[j-1]*paso_r_local;
	for (size_t j=1;j<J+3;j++)	r_local_extended[j] = r_local_extended[j-1]*paso_r_local;
	for (size_t j=0;j<J+1;j++)	delta_r_j_plushalf[j] = r_local[j]*(paso_r_local-1.0);
	for (size_t j=0;j<J+1;j++)	delta_r_j_minushalf[j] = r_local[j]*(1.0-1.0/paso_r_local);
	for (size_t j=0;j<J+1;j++)	delta_r[j] = 0.5*r_local[j]*(paso_r_local-1.0/paso_r_local);
	
	// We define a mesh of points for the Lorentz factor coordinate
	size_t M = 110;
	double g_min = 0.9 * (particle.emin() / (particle.mass*cLight2));
	double g_max = 1.1 * (particle.emax() / (particle.mass*cLight2));
	
	double paso_g_local = pow(g_max/g_min,1.0/M);
	Vector g_local(M+1,g_min);
	Vector g_local_extended(M+3,g_min/paso_g_local);
	Vector delta_g_m_plushalf(M+1,0.0);
	Vector delta_g_m_minushalf(M+1,0.0);
	Vector delta_g(M+1,0.0);
		
	for (size_t m=1;m<M+1;m++)	g_local[m] = g_local[m-1]*paso_g_local;
	for (size_t m=1;m<M+3;m++)	g_local_extended[m] = g_local_extended[m-1]*paso_g_local;
	for (size_t m=0;m<M+1;m++)	delta_g_m_plushalf[m] = g_local[m]*(paso_g_local-1.0);
	for (size_t m=0;m<M+1;m++)	delta_g_m_minushalf[m] = g_local[m]*(1.0-1.0/paso_g_local);
	for (size_t m=0;m<M+1;m++)	delta_g[m] = 0.5*g_local[m]*(paso_g_local-1.0/paso_g_local);
	
	// We define a vector for the time evolution
	
	size_t Ntime = 30;
	size_t NdeltaTime = 20;
	double charEnergy = 1e5*particle.emin();
	double charGamma = charEnergy / (particle.mass*cLight2);
	double charHeight = 100*schwRadius;
	double charRho = massDensityADAF(10*schwRadius);
	double charMagf = magneticField(5*schwRadius);
	double timescaleAdvection = schwRadius / cLight;
	double timescaleDiffusion = diffCoeff_r(charGamma,particle,charHeight,charMagf);
	double timescaleRadial = min(timescaleAdvection, timescaleDiffusion);
	SpaceCoord iR = {0,3,0};
	double timescaleCooling = charEnergy / losses(charEnergy,particle,st,iR);
	double timescaleAcceleration = P2(2.0) / diffCoeff_g(2.0,particle,10*schwRadius,charMagf,charRho);
	double timescaleEnergy = min(timescaleCooling,timescaleAcceleration);
	double naturalTimescale = min(timescaleRadial, timescaleEnergy);
	double maxTime = particle.ps[DIM_R][nR-1]/abs(radialVel(particle.ps[DIM_R][nR-1]));
	maxTime = naturalTimescale*100;
	Vector time(Ntime+1,1.0);
	double paso_time = pow(maxTime/naturalTimescale,1.0/(Ntime-1));
	for (size_t t=1;t<Ntime+1;t++) time[t] = time[t-1]*paso_time;
	
	Vector deltat(Ntime,1.0);
	for (size_t t=0;t<Ntime;t++) deltat[t] = (time[t+1]-time[t])/NdeltaTime;
	
	cout << deltat[Ntime-1] << endl;
	fun2 Afun = [&particle,naturalTimescale] (double r, double g) 
				{
					double rTrue = r * schwRadius;
					double height = height_fun(rTrue);
					double B = magneticField(rTrue);
					double vR = radialVel(rTrue);
					double k_rr = diffCoeff_r(g,particle,height,B);
					return - (vR/schwRadius) * naturalTimescale * ( 2.0*k_rr/(rTrue*vR) + 1.0 );
				};

	fun2 Bfun = [&particle,naturalTimescale] (double r, double g)
				{
					double rTrue = r * schwRadius;
					double height = height_fun(rTrue);
					double B = magneticField(rTrue);
					double k_rr = diffCoeff_r(g,particle,height,B);
					return k_rr / (schwRadius*schwRadius) * naturalTimescale;
				};

	fun2 Cfun = [&particle,&st,paso_r_local,naturalTimescale] (double r, double g) 
				{
					double rTrue = r * schwRadius;
					double E = g*particle.mass*cLight2;
					
					size_t jjR=0;
					if (rTrue > particle.ps[DIM_R][0]) {
						if (rTrue < particle.ps[DIM_R][nR-1]) {
							while (particle.ps[DIM_R][jjR] < rTrue) jjR++;
						} else
							jjR = nR-1;
					}
					
					double dgdt = 0.0;
					if (jjR > 0 && jjR < nR-1) {
						SpaceCoord jR1 = {0,jjR-1,0};
						SpaceCoord jR2 = {0,jjR,0};
						double r1 = st.denf_e.ps[DIM_R][jjR-1];
						double r2 = st.denf_e.ps[DIM_R][jjR];
						double dEdt_1 = losses(E,particle,st,jR1);
						double dEdt_2 = losses(E,particle,st,jR2);
						double slope = safeLog10(dEdt_2/dEdt_1)/safeLog10(r2/r1);
						dgdt = - dEdt_1 * pow(rTrue/r1, slope) / (particle.mass*cLight2);
						
					} else if (jjR == 0) {
						SpaceCoord jR = {0,0,0};
						dgdt = - losses(E,particle,st,jR) / (particle.mass*cLight2);
					} else if (jjR == nR-1) {
						SpaceCoord jR = {0,nR-1,0};
						dgdt = - losses(E,particle,st,jR) / (particle.mass*cLight2);
					}
					double height = height_fun(rTrue);
					double B = magneticField(rTrue);
					double rho = massDensityADAF(rTrue);
					double Dg = diffCoeff_g(g,particle,height,B,rho);
					double vR = radialVel(rTrue);
					double dvRdr = (radialVel(rTrue*paso_r_local) - radialVel(rTrue/paso_r_local)) / 
									(rTrue*(paso_r_local-1.0/paso_r_local));
					
					dvRdr = (dvRdr > 0.0) ? dvRdr : 0.0;
					return (-2.0*Dg/g - dgdt + 1.0/3.0 * (2.0*vR/rTrue + dvRdr) * g ) * naturalTimescale;
				};
	
	fun2 Dfun = [&particle,naturalTimescale] (double r, double g)
				{
					double rTrue = r * schwRadius;
					double height = height_fun(rTrue);
					double B = magneticField(rTrue);
					double rho = massDensityADAF(rTrue);
					double Dg = diffCoeff_g(g,particle,height,B,rho);
					return naturalTimescale * Dg;
				};
	
	fun2 Tfun = [&particle,naturalTimescale] (double r, double g) 
				{
					double rTrue = r * schwRadius;
					double height = height_fun(rTrue);
					double B = magneticField(rTrue);
					double vR = radialVel(rTrue);
					double rateDiff = diffCoeff_r(g,particle,height,B) / (height*height);
					double rateWinds = 2.0*s*abs(vR)/rTrue;
					double rateAccretion = 1.0/accretionTime(rTrue);
					return pow(rateDiff+rateWinds,-1) / naturalTimescale;
				};
		
	fun2 Qfun = [&st,&particle,paso_g_local,factor,naturalTimescale] (double r, double g)
				{
					double rTrue = r * schwRadius;
					double E = g * particle.mass*cLight2;
					double g_inj = 10.0;
					double sigma = g_inj * (paso_g_local-1.0)/2.0;
					double Qlocal = (rTrue > particle.ps[DIM_R][0] && 
							rTrue < particle.ps[DIM_R][nR-1]) ?
							rTrue*rTrue* ionDensity(rTrue) * cLight/schwRadius: 0.0;
					double result = Qlocal * gsl_ran_gaussian_pdf(g-g_inj,sigma);
					return result * naturalTimescale * factor;
				};
	
	Matrix Bj_minusHalf_m;				matrixInit(Bj_minusHalf_m,J+1,M+1,0.0);
	Matrix Bj_plusHalf_m;				matrixInit(Bj_plusHalf_m,J+1,M+1,0.0);
	Matrix Djm_minusHalf;				matrixInit(Djm_minusHalf,J+1,M+1,0.0);
	Matrix Djm_plusHalf;				matrixInit(Djm_plusHalf,J+1,M+1,0.0);
	Matrix Wminus_j_minusHalf;			matrixInit(Wminus_j_minusHalf,J+1,M+1,0.0);
	Matrix Wminus_j_plusHalf;			matrixInit(Wminus_j_plusHalf,J+1,M+1,0.0);
	Matrix Wplus_j_minusHalf;			matrixInit(Wplus_j_minusHalf,J+1,M+1,0.0);
	Matrix Wplus_j_plusHalf;			matrixInit(Wplus_j_plusHalf,J+1,M+1,0.0);
	Matrix Vminus_m_minusHalf;			matrixInit(Vminus_m_minusHalf,J+1,M+1,0.0);
	Matrix Vminus_m_plusHalf;			matrixInit(Vminus_m_plusHalf,J+1,M+1,0.0);
	Matrix Vplus_m_minusHalf;			matrixInit(Vplus_m_minusHalf,J+1,M+1,0.0);
	Matrix Vplus_m_plusHalf;			matrixInit(Vplus_m_plusHalf,J+1,M+1,0.0);
	Matrix Tjm;							matrixInit(Tjm,J+1,M+1,0.0);
	Matrix Qjm;							matrixInit(Qjm,J+1,M+1,0.0);
	
	cout << "TRANSPORT EQUATION: Starting matrix element calculation" << endl << endl;
	
	for (size_t j=0;j<J+1;j++) {
		for (size_t m=0;m<M+1;m++) {
			size_t l = j*(M+1)+m;
			double Ajm = Afun(r_local[j],g_local[m]);
			double Aj_minusHalf_m = 0.5 * (Afun(r_local_extended[j],g_local[m]) + Ajm);
			double Aj_plusHalf_m = 0.5 * (Afun(r_local_extended[j+2],g_local[m]) + Ajm);
					
			double Bjm = Bfun(r_local[j],g_local[m]);
			Bj_minusHalf_m[j][m] = 0.5 * (Bfun(r_local_extended[j],g_local[m]) + Bjm);
			Bj_plusHalf_m[j][m] = 0.5 * (Bfun(r_local_extended[j+2],g_local[m]) + Bjm);
			
			double Cjm = Cfun(r_local[j],g_local[m]);
			double Cjm_minusHalf = 0.5 * (Cfun(r_local[j],g_local_extended[m]) + Cjm);
			double Cjm_plusHalf = 0.5 * (Cfun(r_local[j],g_local_extended[m+2]) + Cjm);
			
			double Djm = Dfun(r_local[j],g_local[m]);
			Djm_minusHalf[j][m] = 0.5 * (Dfun(r_local[j],g_local_extended[m]) + Djm);
			Djm_plusHalf[j][m] = 0.5 * (Dfun(r_local[j],g_local_extended[m+2]) + Djm);
			
			Qjm[j][m] = Qfun(r_local[j],g_local[m]);
			Tjm[j][m] = Tfun(r_local[j],g_local[m]);
					
			double wj_minusHalf = Aj_minusHalf_m/Bj_minusHalf_m[j][m] * delta_r_j_minushalf[j];
			double wj_plusHalf = Aj_plusHalf_m/Bj_plusHalf_m[j][m] * delta_r_j_plushalf[j];
			double Wj_minusHalf(0.0), Wj_plusHalf(0.0);
			
			if ( abs(wj_minusHalf) < 1.0e-3 ) {
				Wj_minusHalf = pow(1.0+gsl_pow_2(wj_minusHalf)/24+gsl_pow_4(wj_minusHalf)/1920,-1);
				Wplus_j_minusHalf[j][m] = Wj_minusHalf * exp(0.5*wj_minusHalf);
				Wminus_j_minusHalf[j][m] = Wj_minusHalf * exp(-0.5*wj_minusHalf);
			} else {
				Wj_minusHalf = abs(wj_minusHalf)*exp(-0.5*abs(wj_minusHalf)) /
								(1.0-exp(-abs(wj_minusHalf)));
				if (wj_minusHalf > 0.0) {
					Wplus_j_minusHalf[j][m] = abs(wj_minusHalf) / (1.0-exp(-abs(wj_minusHalf)));
					Wminus_j_minusHalf[j][m] = Wj_minusHalf * exp(-0.5*wj_minusHalf);
				} else {
					Wplus_j_minusHalf[j][m] = Wj_minusHalf * exp(0.5*wj_minusHalf);
					Wminus_j_minusHalf[j][m] = abs(wj_minusHalf) / (1.0-exp(-abs(wj_minusHalf)));
				}
			}
			if ( abs(wj_plusHalf) < 1.0e-3 ) {
				Wj_plusHalf = pow(1.0+gsl_pow_2(wj_plusHalf)/24+gsl_pow_4(wj_plusHalf)/1920,-1);
				Wplus_j_plusHalf[j][m] = Wj_plusHalf * exp(0.5*wj_plusHalf);
				Wminus_j_plusHalf[j][m] = Wj_plusHalf * exp(-0.5*wj_plusHalf);
			} else {
				Wj_plusHalf = abs(wj_plusHalf)*exp(-0.5*abs(wj_plusHalf)) /
								(1.0-exp(-abs(wj_plusHalf)));
				if (wj_plusHalf > 0.0) {
					Wplus_j_plusHalf[j][m] = abs(wj_plusHalf) / (1.0-exp(-abs(wj_plusHalf)));
					Wminus_j_plusHalf[j][m] = Wj_plusHalf * exp(-0.5*wj_plusHalf);
				} else {
					Wplus_j_plusHalf[j][m] = Wj_plusHalf * exp(0.5*wj_plusHalf);
					Wminus_j_plusHalf[j][m] = abs(wj_plusHalf) / (1.0-exp(-abs(wj_plusHalf)));
				}
			}
			
			double vm_minusHalf = Cjm_minusHalf/Djm_minusHalf[j][m] * delta_g_m_minushalf[m];
			double vm_plusHalf = Cjm_plusHalf/Djm_plusHalf[j][m] * delta_g_m_plushalf[m];
			double Vm_minusHalf(0.0), Vm_plusHalf(0.0);
			
			if ( abs(vm_minusHalf) < 1.0e-3 ) {
				Vm_minusHalf = pow(1.0+gsl_pow_2(vm_minusHalf)/24.0+gsl_pow_4(vm_minusHalf)/1920.0,-1);
				Vplus_m_minusHalf[j][m] = Vm_minusHalf * exp(0.5*vm_minusHalf);
				Vminus_m_minusHalf[j][m] = Vm_minusHalf * exp(-0.5*vm_minusHalf);
			} else {
				Vm_minusHalf = abs(vm_minusHalf)*exp(-0.5*abs(vm_minusHalf)) /
								(1.0-exp(-abs(vm_minusHalf)));
				if (vm_minusHalf > 0.0) {
					Vplus_m_minusHalf[j][m] = abs(vm_minusHalf) / (1.0-exp(-abs(vm_minusHalf)));
					Vminus_m_minusHalf[j][m] = Vm_minusHalf * exp(-0.5*vm_minusHalf);
				} else {
					Vplus_m_minusHalf[j][m] = Vm_minusHalf * exp(0.5*vm_minusHalf);
					Vminus_m_minusHalf[j][m] = abs(vm_minusHalf) / (1.0-exp(-abs(vm_minusHalf)));
				}
			}
			if ( abs(vm_plusHalf) < 1.0e-3 ) {
				Vm_plusHalf = pow(1.0+gsl_pow_2(vm_plusHalf)/24+gsl_pow_4(vm_plusHalf)/1920,-1);
				Vplus_m_plusHalf[j][m] = Vm_plusHalf * exp(0.5*vm_plusHalf);
				Vminus_m_plusHalf[j][m] = Vm_plusHalf * exp(-0.5*vm_plusHalf);
			} else {
				Vm_plusHalf = abs(vm_plusHalf)*exp(-0.5*abs(vm_plusHalf)) /
								(1.0-exp(-abs(vm_plusHalf)));
				if (vm_plusHalf > 0.0) {
					Vplus_m_plusHalf[j][m] = abs(vm_plusHalf) / (1.0-exp(-abs(vm_plusHalf)));
					Vminus_m_plusHalf[j][m] = Vm_plusHalf * exp(-0.5*vm_plusHalf);
				} else {
					Vplus_m_plusHalf[j][m] = Vm_plusHalf * exp(0.5*vm_plusHalf);
					Vminus_m_plusHalf[j][m] = abs(vm_plusHalf) / (1.0-exp(-abs(vm_plusHalf)));
				}
			}
		}
	}
	
	ofstream fileMatrixDistEne, fileMatrixDistPos, fileMatrixDistQT;
	fileMatrixDistEne.open("matrixCoeffTransportEq_Ene.txt");
	fileMatrixDistPos.open("matrixCoeffTransportEq_Pos.txt");
	fileMatrixDistQT.open("matrixCoeffTransportEq_QT.txt");
	fileMatrixDistEne 	<< "j \t m \t Dj- \t Dj+ \t V-- \t V-+ \t V+- \t V++" << endl;
	fileMatrixDistPos 	<< "m \t j \t Bj- \t Bj+ \t W-- \t W-+ \t W+- \t W++" << endl;
	fileMatrixDistQT 	<< "j \t m \t T \t Q \t" << endl;
	
	for (size_t j=0;j<J+1;j++) {
		for (size_t m=0;m<M+1;m++) {
			fileMatrixDistEne 	<< j << "\t" << m << "\t"
								<< Djm_minusHalf[j][m] << "\t" << Djm_plusHalf[j][m] << "\t"
								<< Vminus_m_minusHalf[j][m] << "\t" << Vminus_m_plusHalf[j][m] << "\t"
								<< Vplus_m_minusHalf[j][m] << "\t" << Vplus_m_plusHalf[j][m] << endl;
			fileMatrixDistQT	<< j << "\t" << m << "\t"
								<< Tjm[j][m] << "\t" << Qjm[j][m] << endl;
		}
	}
	for (size_t m=0;m<M+1;m++) {
		for (size_t j=0;j<J+1;j++) {
			fileMatrixDistPos	<< m << "\t" << j << "\t"
								<< Bj_minusHalf_m[j][m] << "\t" << Bj_plusHalf_m[j][m] << "\t"
								<< Wminus_j_minusHalf[j][m] << "\t" << Wminus_j_plusHalf[j][m] << "\t"
								<< Wplus_j_minusHalf[j][m] << "\t" << Wplus_j_plusHalf[j][m] << endl;
		}
	}
 
	fileMatrixDistEne.close();
	fileMatrixDistPos.close();
	fileMatrixDistQT.close();
	
	cout << "TRANSPORT EQUATION: Finished matrix calculation" << endl << endl;
	
	Vector rr(J*(M+1),0.0);
	Matrix rrOld;
	matrixInit(rrOld,J,M+1,0.0);
	Matrix dist;
	
	cout << "TRANSPORT EQUATION: Starting time evolution" << endl << endl;
	
	Vector d((J+1)*(M+1),0.0);
	for (size_t t=0;t<Ntime;t++) {
		cout << "t = " << t << endl;
		for (size_t jt=0;jt<NdeltaTime;jt++) {
			
			Vector a((J+1)*(M+1),0.0);
			Vector ba((J+1)*(M+1),0.0);
			Vector bb((J+1)*(M+1),0.0);
			Vector bc((J+1)*(M+1),0.0);
			Vector c((J+1)*(M+1),0.0);
			
			for (size_t j=0;j<J+1;j++) {
				for (size_t m=0;m<M+1;m++) {
					size_t l = j*(M+1)+m;
					if (j == 0) {
						bb[l] = 1.0;
					} else {
						if (j > 0)
							a[l] = - deltat[t]/delta_r[j] * Bj_minusHalf_m[j][m] * Wminus_j_minusHalf[j][m] / delta_r_j_minushalf[j];
						a[l] /= 2.0;
						
						if (j < J)
							c[l] = - deltat[t]/delta_r[j] * Bj_plusHalf_m[j][m] * Wplus_j_plusHalf[j][m] / delta_r_j_plushalf[j];
						c[l] /= 2.0;
						
						if (m > 0)
							ba[l] = - deltat[t]/delta_g[m] * Djm_minusHalf[j][m] * Vminus_m_minusHalf[j][m] / delta_g_m_minushalf[m];
				
						bb[l] = 1.0 + 0.5 * deltat[t]/delta_r[j] * ( Bj_plusHalf_m[j][m] * Wminus_j_plusHalf[j][m] / delta_r_j_plushalf[j]
								+ Bj_minusHalf_m[j][m] * Wplus_j_minusHalf[j][m] / delta_r_j_minushalf[j] )
						+ 1.0/delta_g[m] * ( Djm_plusHalf[j][m] * Vminus_m_plusHalf[j][m] / delta_g_m_plushalf[m] 
								+ Djm_minusHalf[j][m] * Vplus_m_minusHalf[j][m] / delta_g_m_minushalf[m] ) 
						+ deltat[t]/Tjm[j][m];
						if (m < M)
							bc[l] = - deltat[t]/delta_g[m] * Djm_plusHalf[j][m] * Vplus_m_plusHalf[j][m] / delta_g_m_plushalf[m];
						
						d[l] = Qjm[j][m]*deltat[t] + d[l] + d[l] - a[l] * (j>0 ? d[(j-1)*(M+1)+m] : 0.0) - 
									bb[l] * d[l] - c[l] * (j<J ? d[(j+1)*(M+1)+m] : 0.0);
					}
				}
			}
			TriBlockDiagSys2(a,ba,bb,bc,c,d,J+1,M+1);
		}
		
		if (t % 10 == 0) {
			matrixInit(dist,J+1,M+1,0.0);
			ofstream fileDist;
			fileDist.open("dist"+to_string(t)+".txt");
			for (size_t j=0;j<J+1;j++) {
				for (size_t m=0;m<M+1;m++) {
					size_t l = j*(M+1)+m;
					dist[j][m] = d[l] / P2(r_local[j]*schwRadius) / factor;
					fileDist << r_local[j] << "\t" << g_local[m] << "\t" << dist[j][m] << endl;
				}
			}
			fileDist.close();
		}
	}
	cout << endl << "TRANSPORT EQUATION: Finished time evolution" << endl << endl;
	
	particle.ps.iterate([&] (const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		size_t j = 1;
		while (r_local[j]*schwRadius < r) j++;
		particle.ps.iterate([&] (const SpaceIterator& iRE) {
			double g = iRE.val(DIM_E) / (particle.mass*cLight2);
			size_t m = 0;
			while (g_local[m] < g) m++;
			
			double N11 = dist[j-1][m-1];
			double N12 = dist[j-1][m];
			double N1 = 0.0;
			if (N11 > 0.0 && N12 > 0.0) {
				double s1 = safeLog10(N12/N11)/safeLog10(g_local[m]/g_local[m-1]);
				N1 = N11 * pow(g / g_local[m-1], s1);
			} else {
				if (N11 > 0.0)
					N1 = - N11 * (g - g_local[m-1])/(g_local[m]-g_local[m-1]) + N11;
				else if (N12 > 0.0)
					N1 = N12 * (g - g_local[m-1])/(g_local[m]-g_local[m-1]);
				else
					N1 = 0.0;
			}
			double N21 = dist[j][m-1];
			double N22 = dist[j][m];
			double N2 = 0.0;
			if (N21 > 0.0 && N22 > 0.0) {
				double s2 = safeLog10(N22/N21)/safeLog10(g_local[m]/g_local[m-1]);
				N2 = N21 * pow(g / g_local[m-1], s2);
			} else {
				if (N21 > 0.0)
					N2 = - N21 * (g - g_local[m-1])/(g_local[m]-g_local[m-1]) + N21;
				else if (N22 > 0.0)
					N2 = N22 * (g - g_local[m-1])/(g_local[m]-g_local[m-1]);
				else
					N2 = 0.0;
			}
			
			double N = 0.0;
			if (N1 > 0.0 && N2 > 0.0) {
				double s = safeLog10(N2/N1)/safeLog10(r_local[j]/r_local[j-1]);
				N = N1 * pow(r / (r_local[j-1]*schwRadius), s);
			} else {
				if (N1 > 0.0)
					N = - N1 * (r/schwRadius - r_local[j-1])/(r_local[j]-r_local[j-1]) + N1;
				else if (N2 > 0.0)
					N = N2 * (r/schwRadius - r_local[j-1])/(r_local[j]-r_local[j-1]);
				else
					N = 0.0;
			}
			particle.distribution.set(iRE, N / (particle.mass*cLight2));
		},{-1,iR.coord[DIM_R],0});
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
	
}*/



/*

void distributionFokkerPlanckCompleteSteadyState(Particle& particle, State& st)
// This routine solves the transport equation by considering separated one zones at each
// energy. The escape time includes the time that takes to the particles to "cross" the energy-cell
// by losing energy, and the injection includes the flux of particles coming from the immediately 
// higher-energy cell. It is solved backward (from the highest energy to the lowest one) and so 
// it can only handle systematic energy losses (cooling) and not acceleration.
// At each cell, the steady transport equation is an ordinary diff. equation
// in the radial variable, it considers spatial advection and diffusion and it is solved by 
// finite differences (CC70, PP96).
{
	double factor = 1.0e-35;
	
	size_t J = 20;
	double r_min = sqrt(st.denf_e.ps[DIM_R][0] / schwRadius);
	double r_max = (st.denf_e.ps[DIM_R][nR-1]/schwRadius)*paso_r;
	double paso_r_local = pow(r_max/r_min,1.0/J);
	
	// We define a mesh of points for the radial coordinate
	
	Vector r_local(J+1,r_min);
	Vector r_local_extended(J+3,r_min/paso_r_local);
	Vector delta_r_j_plushalf(J+1,0.0);
	Vector delta_r_j_minushalf(J+1,0.0);
	Vector delta_r(J+1,0.0);
		
	for (size_t j=1;j<J+1;j++)	r_local[j] = r_local[j-1]*paso_r_local;
	for (size_t j=1;j<J+3;j++)	r_local_extended[j] = r_local_extended[j-1]*paso_r_local;
	for (size_t j=0;j<J+1;j++)	delta_r_j_plushalf[j] = r_local[j]*(paso_r_local-1.0);
	for (size_t j=0;j<J+1;j++)	delta_r_j_minushalf[j] = r_local[j]*(1.0-1.0/paso_r_local);
	for (size_t j=0;j<J+1;j++)	delta_r[j] = 0.5*r_local[j]*(paso_r_local-1.0/paso_r_local);
	
	// We define a mesh of points for the Lorentz factor coordinate
	size_t M = 40;
	double g_min = 0.9 * (particle.emin() / (particle.mass*cLight2));
	double g_max = 1.1 * (particle.emax() / (particle.mass*cLight2));
	
	double paso_g_local = pow(g_max/g_min,1.0/M);
	Vector g_local(M+1,g_min);
	Vector g_local_extended(M+3,g_min/paso_g_local);
	Vector delta_g_m_plushalf(M+1,0.0);
	Vector delta_g_m_minushalf(M+1,0.0);
	Vector delta_g(M+1,0.0);
		
	for (size_t m=1;m<M+1;m++)	g_local[m] = g_local[m-1]*paso_g_local;
	for (size_t m=1;m<M+3;m++)	g_local_extended[m] = g_local_extended[m-1]*paso_g_local;
	for (size_t m=0;m<M+1;m++)	delta_g_m_plushalf[m] = g_local[m]*(paso_g_local-1.0);
	for (size_t m=0;m<M+1;m++)	delta_g_m_minushalf[m] = g_local[m]*(1.0-1.0/paso_g_local);
	for (size_t m=0;m<M+1;m++)	delta_g[m] = 0.5*g_local[m]*(paso_g_local-1.0/paso_g_local);
	
	fun2 Afun = [&particle] (double r, double g) 
				{
					double rTrue = r * schwRadius;
					double height = height_fun(rTrue);
					double B = magneticField(rTrue);
					double vR = radialVel(rTrue);
					double k_rr = diffCoeff_r(g,particle,height,B);
					return - (vR/cLight) * ( 2.0*k_rr/(rTrue*vR) + 1.0 );
				};

	fun2 Bfun = [&particle] (double r, double g)
				{
					double rTrue = r * schwRadius;
					double height = height_fun(rTrue);
					double B = magneticField(rTrue);
					double k_rr = diffCoeff_r(g,particle,height,B);
					return k_rr / (schwRadius*cLight);
				};

	fun2 Cfun = [&particle,&st,paso_r_local] (double r, double g) 
				{
					double rTrue = r * schwRadius;
					double E = g*particle.mass*cLight2;
					
					size_t jjR=0;
					if (rTrue > particle.ps[DIM_R][0]) {
						if (rTrue < particle.ps[DIM_R][nR-1]) {
							while (particle.ps[DIM_R][jjR] < rTrue) jjR++;
						} else
							jjR = nR-1;
					}
					
					double dgdt = 0.0;
					if (jjR > 0 && jjR < nR-1) {
						SpaceCoord jR1 = {0,jjR-1,0};
						SpaceCoord jR2 = {0,jjR,0};
						double r1 = st.denf_e.ps[DIM_R][jjR-1];
						double r2 = st.denf_e.ps[DIM_R][jjR];
						double dEdt_1 = losses(E,particle,st,jR1);
						double dEdt_2 = losses(E,particle,st,jR2);
						double slope = safeLog10(dEdt_2/dEdt_1)/safeLog10(r2/r1);
						dgdt = - dEdt_1 * pow(rTrue/r1, slope) / (particle.mass*cLight2);
						
					} else if (jjR == 0) {
						SpaceCoord jR = {0,0,0};
						dgdt = - losses(E,particle,st,jR) / (particle.mass*cLight2);
					} else if (jjR == nR-1) {
						SpaceCoord jR = {0,nR-1,0};
						dgdt = - losses(E,particle,st,jR) / (particle.mass*cLight2);
					}
					double height = height_fun(rTrue);
					double B = magneticField(rTrue);
					double rho = massDensityADAF(rTrue);
					double Dg = diffCoeff_g(g,particle,height,B,rho);
					double vR = radialVel(rTrue);
					double dvRdr = (radialVel(rTrue*paso_r_local) - radialVel(rTrue/paso_r_local)) / 
									(rTrue*(paso_r_local-1.0/paso_r_local));
					
					dvRdr = (dvRdr > 0.0) ? dvRdr : 0.0;
					return (-2.0*Dg/g - dgdt + 1.0/3.0 * (2.0*vR/rTrue + dvRdr) * g ) *(schwRadius/cLight);
				};
	
	fun2 Dfun = [&particle] (double r, double g)
				{
					double rTrue = r * schwRadius;
					double height = height_fun(rTrue);
					double B = magneticField(rTrue);
					double rho = massDensityADAF(rTrue);
					double Dg = diffCoeff_g(g,particle,height,B,rho);
					return (schwRadius/cLight) * Dg;
				};
	
	fun2 Tfun = [&particle] (double r, double g) 
				{
					double rTrue = r * schwRadius;
					double height = height_fun(rTrue);
					double B = magneticField(rTrue);
					double vR = radialVel(rTrue);
					double rateDiff = diffCoeff_r(g,particle,height,B) / (height*height);
					double rateWinds = 2.0*s*abs(vR)/rTrue;
					double rateAccretion = 1.0/accretionTime(rTrue);
					return pow(rateDiff+rateWinds,-1) / (schwRadius/cLight);
				};
		
	fun2 Qfun = [&st,&particle,paso_g_local,factor] (double r, double g)
				{
					double rTrue = r * schwRadius;
					double E = g * particle.mass*cLight2;
					double g_inj = 3.0;
					double sigma = g_inj * (paso_g_local-1.0)/2.0;
					double Qlocal = (rTrue > particle.ps[DIM_R][0] && 
							rTrue < particle.ps[DIM_R][nR-1]) ?
							rTrue*rTrue* ionDensity(rTrue) * cLight/schwRadius: 0.0;
					double result = Qlocal * gsl_ran_gaussian_pdf(g-g_inj,sigma);
					return result * (schwRadius/cLight) * factor;
					//SpaceCoord iAux = {0,0,0};
					//return factor * (schwRadius/cLight) * pow(g,-2) * Qlocal * exp(-g/1e7);
				};
	
	Matrix Bj_minusHalf_m;				matrixInit(Bj_minusHalf_m,J+1,M+1,0.0);
	Matrix Bj_plusHalf_m;				matrixInit(Bj_plusHalf_m,J+1,M+1,0.0);
	Matrix Djm_minusHalf;				matrixInit(Djm_minusHalf,J+1,M+1,0.0);
	Matrix Djm_plusHalf;				matrixInit(Djm_plusHalf,J+1,M+1,0.0);
	Matrix Wminus_j_minusHalf;			matrixInit(Wminus_j_minusHalf,J+1,M+1,0.0);
	Matrix Wminus_j_plusHalf;			matrixInit(Wminus_j_plusHalf,J+1,M+1,0.0);
	Matrix Wplus_j_minusHalf;			matrixInit(Wplus_j_minusHalf,J+1,M+1,0.0);
	Matrix Wplus_j_plusHalf;			matrixInit(Wplus_j_plusHalf,J+1,M+1,0.0);
	Matrix Vminus_m_minusHalf;			matrixInit(Vminus_m_minusHalf,J+1,M+1,0.0);
	Matrix Vminus_m_plusHalf;			matrixInit(Vminus_m_plusHalf,J+1,M+1,0.0);
	Matrix Vplus_m_minusHalf;			matrixInit(Vplus_m_minusHalf,J+1,M+1,0.0);
	Matrix Vplus_m_plusHalf;			matrixInit(Vplus_m_plusHalf,J+1,M+1,0.0);
	Matrix Tjm;							matrixInit(Tjm,J+1,M+1,0.0);
	Matrix Qjm;							matrixInit(Qjm,J+1,M+1,0.0);
	
	cout << "TRANSPORT EQUATION: Starting matrix element calculation" << endl << endl;
	
	for (size_t j=0;j<J+1;j++) {
		for (size_t m=0;m<M+1;m++) {
			double Ajm = Afun(r_local[j],g_local[m]);
			double Aj_minusHalf_m = 0.5 * (Afun(r_local_extended[j],g_local[m]) + Ajm);
			double Aj_plusHalf_m = 0.5 * (Afun(r_local_extended[j+2],g_local[m]) + Ajm);
					
			double Bjm = Bfun(r_local[j],g_local[m]);
			Bj_minusHalf_m[j][m] = 0.5 * (Bfun(r_local_extended[j],g_local[m]) + Bjm);
			Bj_plusHalf_m[j][m] = 0.5 * (Bfun(r_local_extended[j+2],g_local[m]) + Bjm);
			
			double Cjm = Cfun(r_local[j],g_local[m]);
			double Cjm_minusHalf = 0.5 * (Cfun(r_local[j],g_local_extended[m]) + Cjm);
			double Cjm_plusHalf = 0.5 * (Cfun(r_local[j],g_local_extended[m+2]) + Cjm);
			
			double Djm = Dfun(r_local[j],g_local[m]);
			Djm_minusHalf[j][m] = 0.5 * (Dfun(r_local[j],g_local_extended[m]) + Djm);
			Djm_plusHalf[j][m] = 0.5 * (Dfun(r_local[j],g_local_extended[m+2]) + Djm);
			
			Qjm[j][m] = Qfun(r_local[j],g_local[m]);
			Tjm[j][m] = Tfun(r_local[j],g_local[m]);
					
			double wj_minusHalf = Aj_minusHalf_m/Bj_minusHalf_m[j][m] * delta_r_j_minushalf[j];
			double wj_plusHalf = Aj_plusHalf_m/Bj_plusHalf_m[j][m] * delta_r_j_plushalf[j];
			double Wj_minusHalf(0.0), Wj_plusHalf(0.0);
			
			if ( abs(wj_minusHalf) < 1.0e-3 ) {
				Wj_minusHalf = pow(1.0+gsl_pow_2(wj_minusHalf)/24+gsl_pow_4(wj_minusHalf)/1920,-1);
				Wplus_j_minusHalf[j][m] = Wj_minusHalf * exp(0.5*wj_minusHalf);
				Wminus_j_minusHalf[j][m] = Wj_minusHalf * exp(-0.5*wj_minusHalf);
			} else {
				Wj_minusHalf = abs(wj_minusHalf)*exp(-0.5*abs(wj_minusHalf)) /
								(1.0-exp(-abs(wj_minusHalf)));
				if (wj_minusHalf > 0.0) {
					Wplus_j_minusHalf[j][m] = abs(wj_minusHalf) / (1.0-exp(-abs(wj_minusHalf)));
					Wminus_j_minusHalf[j][m] = Wj_minusHalf * exp(-0.5*wj_minusHalf);
				} else {
					Wplus_j_minusHalf[j][m] = Wj_minusHalf * exp(0.5*wj_minusHalf);
					Wminus_j_minusHalf[j][m] = abs(wj_minusHalf) / (1.0-exp(-abs(wj_minusHalf)));
				}
			}
			if ( abs(wj_plusHalf) < 1.0e-3 ) {
				Wj_plusHalf = pow(1.0+gsl_pow_2(wj_plusHalf)/24+gsl_pow_4(wj_plusHalf)/1920,-1);
				Wplus_j_plusHalf[j][m] = Wj_plusHalf * exp(0.5*wj_plusHalf);
				Wminus_j_plusHalf[j][m] = Wj_plusHalf * exp(-0.5*wj_plusHalf);
			} else {
				Wj_plusHalf = abs(wj_plusHalf)*exp(-0.5*abs(wj_plusHalf)) /
								(1.0-exp(-abs(wj_plusHalf)));
				if (wj_plusHalf > 0.0) {
					Wplus_j_plusHalf[j][m] = abs(wj_plusHalf) / (1.0-exp(-abs(wj_plusHalf)));
					Wminus_j_plusHalf[j][m] = Wj_plusHalf * exp(-0.5*wj_plusHalf);
				} else {
					Wplus_j_plusHalf[j][m] = Wj_plusHalf * exp(0.5*wj_plusHalf);
					Wminus_j_plusHalf[j][m] = abs(wj_plusHalf) / (1.0-exp(-abs(wj_plusHalf)));
				}
			}
			
			double vm_minusHalf = Cjm_minusHalf/Djm_minusHalf[j][m] * delta_g_m_minushalf[m];
			double vm_plusHalf = Cjm_plusHalf/Djm_plusHalf[j][m] * delta_g_m_plushalf[m];
			double Vm_minusHalf(0.0), Vm_plusHalf(0.0);
			
			if ( abs(vm_minusHalf) < 1.0e-3 ) {
				Vm_minusHalf = pow(1.0+gsl_pow_2(vm_minusHalf)/24.0+gsl_pow_4(vm_minusHalf)/1920.0,-1);
				Vplus_m_minusHalf[j][m] = Vm_minusHalf * exp(0.5*vm_minusHalf);
				Vminus_m_minusHalf[j][m] = Vm_minusHalf * exp(-0.5*vm_minusHalf);
			} else {
				Vm_minusHalf = abs(vm_minusHalf)*exp(-0.5*abs(vm_minusHalf)) /
								(1.0-exp(-abs(vm_minusHalf)));
				if (vm_minusHalf > 0.0) {
					Vplus_m_minusHalf[j][m] = abs(vm_minusHalf) / (1.0-exp(-abs(vm_minusHalf)));
					Vminus_m_minusHalf[j][m] = Vm_minusHalf * exp(-0.5*vm_minusHalf);
				} else {
					Vplus_m_minusHalf[j][m] = Vm_minusHalf * exp(0.5*vm_minusHalf);
					Vminus_m_minusHalf[j][m] = abs(vm_minusHalf) / (1.0-exp(-abs(vm_minusHalf)));
				}
			}
			if ( abs(vm_plusHalf) < 1.0e-3 ) {
				Vm_plusHalf = pow(1.0+gsl_pow_2(vm_plusHalf)/24+gsl_pow_4(vm_plusHalf)/1920,-1);
				Vplus_m_plusHalf[j][m] = Vm_plusHalf * exp(0.5*vm_plusHalf);
				Vminus_m_plusHalf[j][m] = Vm_plusHalf * exp(-0.5*vm_plusHalf);
			} else {
				Vm_plusHalf = abs(vm_plusHalf)*exp(-0.5*abs(vm_plusHalf)) /
								(1.0-exp(-abs(vm_plusHalf)));
				if (vm_plusHalf > 0.0) {
					Vplus_m_plusHalf[j][m] = abs(vm_plusHalf) / (1.0-exp(-abs(vm_plusHalf)));
					Vminus_m_plusHalf[j][m] = Vm_plusHalf * exp(-0.5*vm_plusHalf);
				} else {
					Vplus_m_plusHalf[j][m] = Vm_plusHalf * exp(0.5*vm_plusHalf);
					Vminus_m_plusHalf[j][m] = abs(vm_plusHalf) / (1.0-exp(-abs(vm_plusHalf)));
				}
			}
		}
	}
	
	ofstream fileMatrixDistEne, fileMatrixDistPos, fileMatrixDistQT;
	fileMatrixDistEne.open("matrixCoeffTransportEq_Ene.txt");
	fileMatrixDistPos.open("matrixCoeffTransportEq_Pos.txt");
	fileMatrixDistQT.open("matrixCoeffTransportEq_QT.txt");
	fileMatrixDistEne 	<< "j \t m \t Dj- \t Dj+ \t V-- \t V-+ \t V+- \t V++" << endl;
	fileMatrixDistPos 	<< "m \t j \t Bj- \t Bj+ \t W-- \t W-+ \t W+- \t W++" << endl;
	fileMatrixDistQT 	<< "j \t m \t T \t Q \t" << endl;
	
	for (size_t j=0;j<J+1;j++) {
		for (size_t m=0;m<M+1;m++) {
			fileMatrixDistEne 	<< j << "\t" << m << "\t"
								<< Djm_minusHalf[j][m] << "\t" << Djm_plusHalf[j][m] << "\t"
								<< Vminus_m_minusHalf[j][m] << "\t" << Vminus_m_plusHalf[j][m] << "\t"
								<< Vplus_m_minusHalf[j][m] << "\t" << Vplus_m_plusHalf[j][m] << endl;
			fileMatrixDistQT	<< j << "\t" << m << "\t"
								<< Tjm[j][m] << "\t" << Qjm[j][m] << endl;
		}
	}
	for (size_t m=0;m<M+1;m++) {
		for (size_t j=0;j<J+1;j++) {
			fileMatrixDistPos	<< m << "\t" << j << "\t"
								<< Bj_minusHalf_m[j][m] << "\t" << Bj_plusHalf_m[j][m] << "\t"
								<< Wminus_j_minusHalf[j][m] << "\t" << Wminus_j_plusHalf[j][m] << "\t"
								<< Wplus_j_minusHalf[j][m] << "\t" << Wplus_j_plusHalf[j][m] << endl;
		}
	}
 
	fileMatrixDistEne.close();
	fileMatrixDistPos.close();
	fileMatrixDistQT.close();
	
	cout << "TRANSPORT EQUATION: Finished matrix calculation" << endl << endl;
	
	Vector x((J+1)*(M+1),0.0);
	Matrix dist;
	
	cout << "TRANSPORT EQUATION: Starting solving" << endl << endl;
	
	Vector a((J+1)*(M+1),0.0);
	Vector ba((J+1)*(M+1),0.0);
	Vector bb((J+1)*(M+1),0.0);
	Vector bc((J+1)*(M+1),0.0);
	Vector c((J+1)*(M+1),0.0);
	Vector d((J+1)*(M+1),0.0);
	
	for (size_t j=0;j<J+1;j++) {
		for (size_t m=0;m<M+1;m++) {
			size_t l = j*(M+1)+m;
			if (j == 0) {
				bb[l] = 1.0;
			} else {
				if (j > 0)
					a[l] = - 1.0/delta_r[j] * Bj_minusHalf_m[j][m] * Wminus_j_minusHalf[j][m] / delta_r_j_minushalf[j];
				if (j < J)
					c[l] = - 1.0/delta_r[j] * Bj_plusHalf_m[j][m] * Wplus_j_plusHalf[j][m] / delta_r_j_plushalf[j];
				
				if (m > 0)
					ba[l] = - 1.0/delta_g[m] * Djm_minusHalf[j][m] * Vminus_m_minusHalf[j][m] / delta_g_m_minushalf[m];
				
				bb[l] = 1.0/delta_r[j] * ( Bj_plusHalf_m[j][m] * Wminus_j_plusHalf[j][m] / delta_r_j_plushalf[j]
								+ Bj_minusHalf_m[j][m] * Wplus_j_minusHalf[j][m] / delta_r_j_minushalf[j] )
						+ 1.0/delta_g[m] * ( Djm_plusHalf[j][m] * Vminus_m_plusHalf[j][m] / delta_g_m_plushalf[m] 
								+ Djm_minusHalf[j][m] * Vplus_m_minusHalf[j][m] / delta_g_m_minushalf[m] ) 
						+ 1.0/Tjm[j][m];
				if (m < M)
					bc[l] = - 1.0/delta_g[m] * Djm_plusHalf[j][m] * Vplus_m_plusHalf[j][m] / delta_g_m_plushalf[m];
				d[l] = Qjm[j][m];
			}
		}
	}
	
	ofstream fileElements;
	fileElements.open("elements.txt");
	fileElements << "j \t m \t a \t ba \t bb \t bc \t c \t d" << endl;
	for (size_t j=0;j<J+1;j++) {
		for (size_t m=0;m<M+1;m++) {
			size_t l = j*(M+1)+m;
			fileElements << j << "\t" << m << "\t" << a[l] << "\t" <<
			ba[l] << "\t" << bb[l] << "\t" << bc[l] << "\t" << c[l] << "\t" << d[l] << endl;
		}
	}
	fileElements.close();
	
	TriBlockDiagSys(a,ba,bb,bc,c,d,J+1,M+1);
	matrixInit(dist,J+1,M+1,0.0);
	ofstream fileDist;
	fileDist.open("dist.txt");
	for (size_t j=0;j<J+1;j++) {
		for (size_t m=0;m<M+1;m++) {
			dist[j][m] = d[j*(M+1)+m] / P2(r_local[j]*schwRadius) / factor;
			fileDist << r_local[j] << "\t" << g_local[m] << "\t" << dist[j][m] << endl;
		}
	}
	fileDist.close();*/
	/*
	int nz_num = (M+1)*(J+1) + 2*J*M + J*(M+1) + (J-1)*(M+1);
	int ia[nz_num], ja[nz_num];
	double aa[nz_num], xx[(M+1)*(J+1)];
	
	cout << "TRANSPORT EQUATION: Assigning elements to the matrix" << endl;
	
	
	
	cout << endl;
	cout << "bb" << endl << endl;
	for (size_t l=0;l<(M+1)*(J+1);l++) {
		ia[l] = l;
		ja[l] = l;
		aa[l] = bb[l];
		cout << l << "\t" << ia[l] << "\t" << ja[l] << "\t" << aa[l] << "\t" << bb[l] << endl;
	}
	cout << endl;
	cout << "bc" << endl << endl;
	for (size_t l=(M+1)*(J+1);l<(M+1)*(J+1)+J*M;l++) {
		size_t lp = l - (M+1)*(J+1);
		ia[l] = M+1 + lp + lp/M;
		ja[l] = ia[l] + 1;
		aa[l] = bc[M+1 + lp + lp/M];
		cout << l << "\t" << ia[l] << "\t" << ja[l] << "\t" << aa[l] << "\t" << bc[M+1 + lp + lp/M] << endl;
	}
	cout << endl;
	cout << "ba" << endl << endl;
	for (size_t l=(M+1)*(J+1)+J*M;l<(M+1)*(J+1)+2*J*M;l++) {
		size_t lp = l - ( (M+1)*(J+1) + J*M );
		ia[l] = M+1 + 1 + lp + lp/M;
		ja[l] = ia[l] - 1;
		aa[l] = ba[M+1 + 1 + lp + lp/M];
		cout << l << "\t" << ia[l] << "\t" << ja[l] << "\t" << aa[l] << "\t" << ba[M+1+1+lp+lp/M] << endl;
	}
	cout << "a" << endl << endl;
	for (size_t l=(M+1)*(J+1)+2*J*M;l<(M+1)*(J+1)+2*J*M+J*(M+1);l++) {
		size_t lp = l - ( (M+1)*(J+1) + 2*J*M );
		ia[l] = M+1 + lp;
		ja[l] = lp;
		aa[l] = a[M+1 + lp];
		cout << l << "\t" << ia[l] << "\t" << ja[l] << "\t" << aa[l] << "\t" << a[M+1+lp] << endl;
	}
	cout << "c" << endl << endl;
	for (size_t l=(M+1)*(J+1)+2*J*M+J*(M+1);l<(M+1)*(J+1)+2*J*M+J*(M+1)+(J-1)*(M+1);l++) {
		size_t lp = l - ( (M+1)*(J+1) + 2*J*M + J*(M+1) );
		ia[l] = M+1 + lp;
		ja[l] = lp + 2*(M+1);
		aa[l] = c[M+1 + lp];
		cout << l << "\t" << ia[l] << "\t" << ja[l] << "\t" << aa[l] << "\t" << c[M+1+lp] << endl;
	}
	
	cout << "TRANSPORT EQUATION: Start Gauss-Seidel" << endl;
	//matrixWrite("matrixElements.txt",AA,(J+1)*(M+1),(J+1)*(M+1));
	
	Matrix AA;
	matrixInit(AA,(J+1)*(M+1),(J+1)*(M+1),0.0);
	
	for (size_t l=0;l<nz_num;l++) AA[ia[l]][ja[l]] = aa[l];
	
	mgmres_st((M+1)*(J+1),nz_num,ia,ja,aa,xx,xx,10,(M+1)*(J+1)-1,1e-1,1e-3);
	double err = 0.0;
	
	//Jacobi(a,ba,bb,bc,c,d,J+1,M+1);
	//GaussSeidel(AA,d,d,(J+1)*(M+1),1,err);

	matrixInit(dist,J+1,M+1,0.0);
	ofstream fileDist;
	fileDist.open("dist.txt");
	for (size_t j=0;j<J+1;j++) {
		for (size_t m=0;m<M+1;m++) {
			dist[j][m] = x[j*(M+1)+m] / P2(r_local[j]*schwRadius) / factor;
			fileDist << r_local[j] << "\t" << g_local[m] << "\t" << dist[j][m] << endl;
		}
	}
	fileDist.close();
	
	cout << endl << "TRANSPORT EQUATION: Finished" << endl << endl;
	
	particle.ps.iterate([&] (const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		size_t j = 1;
		while (r_local[j]*schwRadius < r) j++;
		particle.ps.iterate([&] (const SpaceIterator& iRE) {
			double g = iRE.val(DIM_E) / (particle.mass*cLight2);
			size_t m = 0;
			while (g_local[m] < g) m++;
			
			double N11 = dist[j-1][m-1];
			double N12 = dist[j-1][m];
			double N1 = 0.0;
			if (N11 > 0.0 && N12 > 0.0) {
				double s1 = safeLog10(N12/N11)/safeLog10(g_local[m]/g_local[m-1]);
				N1 = N11 * pow(g / g_local[m-1], s1);
			} else {
				if (N11 > 0.0)
					N1 = - N11 * (g - g_local[m-1])/(g_local[m]-g_local[m-1]) + N11;
				else if (N12 > 0.0)
					N1 = N12 * (g - g_local[m-1])/(g_local[m]-g_local[m-1]);
				else
					N1 = 0.0;
			}
			double N21 = dist[j][m-1];
			double N22 = dist[j][m];
			double N2 = 0.0;
			if (N21 > 0.0 && N22 > 0.0) {
				double s2 = safeLog10(N22/N21)/safeLog10(g_local[m]/g_local[m-1]);
				N2 = N21 * pow(g / g_local[m-1], s2);
			} else {
				if (N21 > 0.0)
					N2 = - N21 * (g - g_local[m-1])/(g_local[m]-g_local[m-1]) + N21;
				else if (N22 > 0.0)
					N2 = N22 * (g - g_local[m-1])/(g_local[m]-g_local[m-1]);
				else
					N2 = 0.0;
			}
			
			double N = 0.0;
			if (N1 > 0.0 && N2 > 0.0) {
				double s = safeLog10(N2/N1)/safeLog10(r_local[j]/r_local[j-1]);
				N = N1 * pow(r / (r_local[j-1]*schwRadius), s);
			} else {
				if (N1 > 0.0)
					N = - N1 * (r/schwRadius - r_local[j-1])/(r_local[j]-r_local[j-1]) + N1;
				else if (N2 > 0.0)
					N = N2 * (r/schwRadius - r_local[j-1])/(r_local[j]-r_local[j-1]);
				else
					N = 0.0;
			}
			particle.distribution.set(iRE, N / (particle.mass*cLight2));
		},{-1,iR.coord[DIM_R],0});
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
*/





void distributionSpatialDiffusion(Particle& particle, State& st)
// This routine solves the complete time-dependent transport equation considering turbulence
// acceleration and spatial diffusion.
{
	double factor = 1.0;
	
	size_t J = 100;
	double r_min = sqrt(st.denf_e.ps[DIM_R][0] / schwRadius);
	double r_max = (st.denf_e.ps[DIM_R][nR-1] / schwRadius) * paso_r;
	double paso_r_local = pow(r_max/r_min,1.0/J);
	
	// We define a mesh of points for the radial coordinate
	
	Vector r_local(J+1,r_min);
	Vector r_local_extended(J+3,r_min/paso_r_local);
	Vector delta_r_j_plushalf(J+1,0.0);
	Vector delta_r_j_minushalf(J+1,0.0);
	Vector delta_r(J+1,0.0);
		
	for (size_t j=1;j<J+1;j++)	r_local[j] = r_local[j-1]*paso_r_local;
	for (size_t j=1;j<J+3;j++)	r_local_extended[j] = r_local_extended[j-1]*paso_r_local;
	for (size_t j=0;j<J+1;j++)	delta_r_j_plushalf[j] = r_local[j]*(paso_r_local-1.0);
	for (size_t j=0;j<J+1;j++)	delta_r_j_minushalf[j] = r_local[j]*(1.0-1.0/paso_r_local);
	for (size_t j=0;j<J+1;j++)	delta_r[j] = 0.5*r_local[j]*(paso_r_local-1.0/paso_r_local);
	
	// We define a mesh of points for the Lorentz factor coordinate
	size_t M = 120;
	double g_min = 0.9 * (particle.emin() / (particle.mass*cLight2));
	double g_max = 1.1 * (particle.emax() / (particle.mass*cLight2));
	
	double paso_g_local = pow(g_max/g_min,1.0/M);
	Vector g_local(M+1,g_min);
	Vector g_local_extended(M+3,g_min/paso_g_local);
	Vector delta_g_m_plushalf(M+1,0.0);
	Vector delta_g_m_minushalf(M+1,0.0);
	Vector delta_g(M+1,0.0);
		
	for (size_t m=1;m<M+1;m++)	g_local[m] = g_local[m-1]*paso_g_local;
	for (size_t m=1;m<M+3;m++)	g_local_extended[m] = g_local_extended[m-1]*paso_g_local;
	for (size_t m=0;m<M+1;m++)	delta_g_m_plushalf[m] = g_local[m]*(paso_g_local-1.0);
	for (size_t m=0;m<M+1;m++)	delta_g_m_minushalf[m] = g_local[m]*(1.0-1.0/paso_g_local);
	for (size_t m=0;m<M+1;m++)	delta_g[m] = 0.5*g_local[m]*(paso_g_local-1.0/paso_g_local);
	
	// We define a vector for the time evolution
	
	double naturalTimescale = schwRadius/cLight;
	cout << "Natural timescale = " << naturalTimescale << endl;
	
	size_t Ntime = 5e2;
	size_t NdeltaTime = 1e2;
	
	Vector time(Ntime+1,1.0);
	double energyTimescale = particle.emax()/lossesSyn(particle.emax(),magneticField(10*schwRadius),particle);
	double dt = naturalTimescale;
	//double paso_time = pow(maxTime/naturalTimescale,1.0/(Ntime-1));
	//for (size_t t=1;t<Ntime+1;t++) time[t] = time[t-1]*paso_time;
	for (size_t t=1;t<Ntime+1;t++) time[t] = time[t-1] + dt;
	Vector deltat(Ntime,1.0);
	for (size_t t=0;t<Ntime;t++) deltat[t] = (time[t+1]-time[t])/NdeltaTime;
	
	fun2 Afun = [&particle,naturalTimescale] (double r, double g) 
				{
					double rTrue = r * schwRadius;
					double height = height_fun(rTrue);
					double B = magneticField(rTrue);
					double vR = radialVel(rTrue);
					double k_rr = diffCoeff_r(g,particle,height,B);
					return - (vR/schwRadius) * naturalTimescale * ( 2.0*k_rr/(rTrue*vR) + 1.0 );
				};

	fun2 Bfun = [&particle,naturalTimescale] (double r, double g)
				{
					double rTrue = r * schwRadius;
					double height = height_fun(rTrue);
					double B = magneticField(rTrue);
					double k_rr = diffCoeff_r(g,particle,height,B);
					return k_rr / (schwRadius*schwRadius) * naturalTimescale;
				};

	fun2 Cfun = [&particle,&st,paso_r_local,naturalTimescale] (double r, double g) 
				{
					double rTrue = r * schwRadius;
					double E = g*particle.mass*cLight2;
					
					size_t jjR=0;
					if (rTrue > particle.ps[DIM_R][0]) {
						if (rTrue < particle.ps[DIM_R][nR-1]) {
							while (particle.ps[DIM_R][jjR] < rTrue) jjR++;
						} else
							jjR = nR-1;
					}
					
					double dgdt = 0.0;
					if (jjR > 0 && jjR < nR-1) {
						SpaceCoord jR1 = {0,jjR-1,0};
						SpaceCoord jR2 = {0,jjR,0};
						double r1 = st.denf_e.ps[DIM_R][jjR-1];
						double r2 = st.denf_e.ps[DIM_R][jjR];
						double dEdt_1 = losses(E,particle,st,jR1);
						double dEdt_2 = losses(E,particle,st,jR2);
						double slope = safeLog10(dEdt_2/dEdt_1)/safeLog10(r2/r1);
						dgdt = - dEdt_1 * pow(rTrue/r1, slope) / (particle.mass*cLight2);
						
					} else if (jjR == 0) {
						SpaceCoord jR = {0,0,0};
						dgdt = - losses(E,particle,st,jR) / (particle.mass*cLight2);
					} else if (jjR == nR-1) {
						SpaceCoord jR = {0,nR-1,0};
						dgdt = - losses(E,particle,st,jR) / (particle.mass*cLight2);
					}
					double vR = radialVel(rTrue);
					double dvRdr = (radialVel(rTrue*paso_r_local) - radialVel(rTrue/paso_r_local)) / 
									(rTrue*(paso_r_local-1.0/paso_r_local));
					
					dvRdr = (dvRdr > 0.0) ? dvRdr : 0.0;
					dvRdr = 0.0;
					return naturalTimescale * ( - dgdt );
					//return naturalTimescale * (-dgdt*g + 1.0/3.0*(2.0*vR/rTrue+dvRdr)*g);
				};
	
	fun2 Tfun = [&particle,naturalTimescale] (double r, double g) 
				{
					double rTrue = r * schwRadius;
					double height = height_fun(rTrue);
					double B = magneticField(rTrue);
					double vR = radialVel(rTrue);
					double rateDiff = diffCoeff_r(g,particle,height,B) / (height*height);
					double rateWinds = 2.0*s*abs(vR)/rTrue;
					return pow(rateDiff+rateWinds,-1) / naturalTimescale;
				};
		
	fun2 Qfun = [&st,&particle,paso_g_local,factor,naturalTimescale,r_local,J] (double r, double g)
				{
					SpaceCoord i = {0,0,0};
					double rTrue = r * schwRadius;
					double E = g * particle.mass*cLight2;
					double Q = (E > particle.emin() && E < particle.emax() &&
								rTrue > particle.ps[DIM_R].first() && rTrue < particle.ps[DIM_R].last()) ?
								particle.injection.interpolate({{DIM_E,E},{DIM_R,rTrue}},&i) : 0.0;
					return factor*Q*r*r*naturalTimescale * particle.mass*cLight2;
				};
	
	Matrix Bj_minusHalf_m;				matrixInit(Bj_minusHalf_m,M+1,J+1,0.0);
	Matrix Bj_plusHalf_m;				matrixInit(Bj_plusHalf_m,M+1,J+1,0.0);
	Matrix Cjm;							matrixInit(Cjm,J+1,M+1,0.0);
	Matrix Cjm_plusHalf;				matrixInit(Cjm_plusHalf,J+1,M+1,0.0);
	Matrix Cjm_minusHalf;				matrixInit(Cjm_minusHalf,J+1,M+1,0.0);
	Matrix Wminus_j_minusHalf;			matrixInit(Wminus_j_minusHalf,M+1,J+1,0.0);
	Matrix Wminus_j_plusHalf;			matrixInit(Wminus_j_plusHalf,M+1,J+1,0.0);
	Matrix Wplus_j_minusHalf;			matrixInit(Wplus_j_minusHalf,M+1,J+1,0.0);
	Matrix Wplus_j_plusHalf;			matrixInit(Wplus_j_plusHalf,M+1,J+1,0.0);
	Matrix Tjm;							matrixInit(Tjm,J+1,M+1,0.0);
	Matrix Qjm;							matrixInit(Qjm,J+1,M+1,0.0);
	
	cout << "TRANSPORT EQUATION: Starting matrix element calculation" << endl << endl;
	
	#pragma omp parallel for
	for (size_t j=0;j<J+1;j++) {
		for (size_t m=0;m<M+1;m++) {
			size_t l = j*(M+1)+m;
			double Ajm = Afun(r_local[j],g_local[m]);
			double Aj_minusHalf_m = 0.5 * (Afun(r_local_extended[j],g_local[m]) + Ajm);
			double Aj_plusHalf_m = 0.5 * (Afun(r_local_extended[j+2],g_local[m]) + Ajm);
					
			double Bjm = Bfun(r_local[j],g_local[m]);
			Bj_minusHalf_m[m][j] = 0.5 * (Bfun(r_local_extended[j],g_local[m]) + Bjm);
			Bj_plusHalf_m[m][j] = 0.5 * (Bfun(r_local_extended[j+2],g_local[m]) + Bjm);
			
			Cjm[j][m] = Cfun(r_local[j],g_local[m]);
			Cjm_minusHalf[j][m] = 0.5 * (Cfun(r_local[j],g_local_extended[m]) + Cjm[j][m]);
			Cjm_plusHalf[j][m] = 0.5 * (Cfun(r_local[j],g_local_extended[m+2]) + Cjm[j][m]);
			
			Qjm[j][m] = Qfun(r_local[j],g_local[m]);
			Tjm[j][m] = Tfun(r_local[j],g_local[m]);
					
			double wj_minusHalf = Aj_minusHalf_m/Bj_minusHalf_m[m][j] * delta_r_j_minushalf[j];
			double wj_plusHalf = Aj_plusHalf_m/Bj_plusHalf_m[m][j] * delta_r_j_plushalf[j];
			double Wj_minusHalf(0.0), Wj_plusHalf(0.0);
			
			if ( abs(wj_minusHalf) < 1.0e-3 ) {
				Wj_minusHalf = pow(1.0+gsl_pow_2(wj_minusHalf)/24+gsl_pow_4(wj_minusHalf)/1920,-1);
				Wplus_j_minusHalf[m][j] = Wj_minusHalf * exp(0.5*wj_minusHalf);
				Wminus_j_minusHalf[m][j] = Wj_minusHalf * exp(-0.5*wj_minusHalf);
			} else {
				Wj_minusHalf = abs(wj_minusHalf)*exp(-0.5*abs(wj_minusHalf)) /
								(1.0-exp(-abs(wj_minusHalf)));
				if (wj_minusHalf > 0.0) {
					Wplus_j_minusHalf[m][j] = abs(wj_minusHalf) / (1.0-exp(-abs(wj_minusHalf)));
					Wminus_j_minusHalf[m][j] = Wj_minusHalf * exp(-0.5*wj_minusHalf);
				} else {
					Wplus_j_minusHalf[m][j] = Wj_minusHalf * exp(0.5*wj_minusHalf);
					Wminus_j_minusHalf[m][j] = abs(wj_minusHalf) / (1.0-exp(-abs(wj_minusHalf)));
				}
			}
			if ( abs(wj_plusHalf) < 1.0e-3 ) {
				Wj_plusHalf = pow(1.0+gsl_pow_2(wj_plusHalf)/24+gsl_pow_4(wj_plusHalf)/1920,-1);
				Wplus_j_plusHalf[m][j] = Wj_plusHalf * exp(0.5*wj_plusHalf);
				Wminus_j_plusHalf[m][j] = Wj_plusHalf * exp(-0.5*wj_plusHalf);
			} else {
				Wj_plusHalf = abs(wj_plusHalf)*exp(-0.5*abs(wj_plusHalf)) /
								(1.0-exp(-abs(wj_plusHalf)));
				if (wj_plusHalf > 0.0) {
					Wplus_j_plusHalf[m][j] = abs(wj_plusHalf) / (1.0-exp(-abs(wj_plusHalf)));
					Wminus_j_plusHalf[m][j] = Wj_plusHalf * exp(-0.5*wj_plusHalf);
				} else {
					Wplus_j_plusHalf[m][j] = Wj_plusHalf * exp(0.5*wj_plusHalf);
					Wminus_j_plusHalf[m][j] = abs(wj_plusHalf) / (1.0-exp(-abs(wj_plusHalf)));
				}
			}
		}
	}
	
	ofstream fileMatrixDistEne, fileMatrixDistPos, fileMatrixDistQT;
	fileMatrixDistEne.open("matrixCoeffTransportEq_Ene.txt");
	fileMatrixDistPos.open("matrixCoeffTransportEq_Pos.txt");
	fileMatrixDistQT.open("matrixCoeffTransportEq_QT.txt");
	fileMatrixDistEne 	<< "j \t m \t Cj- \t Cj+" << endl;
	fileMatrixDistPos 	<< "m \t j \t Bj- \t Bj+ \t W-- \t W-+ \t W+- \t W++" << endl;
	fileMatrixDistQT 	<< "j \t m \t T \t Q \t" << endl;
	
	for (size_t j=0;j<J+1;j++) {
		for (size_t m=0;m<M+1;m++) {
			fileMatrixDistEne 	<< j << "\t" << m << "\t"
								<< Cjm_minusHalf[j][m] << "\t" << Cjm_plusHalf[j][m] << endl;
			fileMatrixDistQT	<< j << "\t" << m << "\t"
								<< Tjm[j][m] << "\t" << Qjm[j][m] << endl;
		}
	}
	for (size_t m=0;m<M+1;m++) {
		for (size_t j=0;j<J+1;j++) {
			fileMatrixDistPos	<< m << "\t" << j << "\t"
								<< Bj_minusHalf_m[m][j] << "\t" << Bj_plusHalf_m[m][j] << "\t"
								<< Wminus_j_minusHalf[m][j] << "\t" << Wminus_j_plusHalf[m][j] << "\t"
								<< Wplus_j_minusHalf[m][j] << "\t" << Wplus_j_plusHalf[m][j] << endl;
		}
	}
 
	fileMatrixDistEne.close();
	fileMatrixDistPos.close();
	fileMatrixDistQT.close();
	
	cout << "TRANSPORT EQUATION: Finished matrix calculation" << endl << endl;
	
	Vector rr((J+1)*(M+1),0.0);
	Matrix rrOld;
	matrixInit(rrOld,J+1,M+1,0.0);
	Matrix dist;
	
	cout << "TRANSPORT EQUATION: Starting time evolution" << endl << endl;
	
	for (size_t t=0;t<Ntime;t++) {
		cout << "t = " << t << endl;
		for (size_t jt=0;jt<NdeltaTime;jt++) {
			
			
			Vector d((J+1)*(M+1),0.0);
			Vector e((J+1)*(M+1),0.0);
			Vector f((J+1)*(M+1),0.0);
			
			double delta_t = deltat[t];
			
			#pragma omp parallel for
			for (size_t j=0;j<J+1;j++) {
				for (size_t m=0;m<M+1;m++) {
					if (delta_t > delta_g[m]/abs(Cjm[j][m]) && abs(Cjm[j][m]) > 0.0) delta_t = delta_g[m]/abs(Cjm[j][m]); 
				}
				delta_t = deltat[t];
				for (size_t m=0;m<M+1;m++) {
					size_t lp = j*(M+1)+m;
					
					if (m > 0 && Cjm_minusHalf[j][m] < 0.0)
						d[lp] = deltat[t]/delta_g[m] * Cjm_minusHalf[j][m];
					if (m < M && Cjm_plusHalf[j][m] > 0.0)
						f[lp] = - delta_t/delta_g[m] * Cjm_plusHalf[j][m];
					e[lp] = 1.0 + delta_t/delta_g[m] * Cjm_minusHalf[j][m];
					
					rr[lp] = rrOld[j][m] + (j > 0 && j < J ? Qjm[j][m] * delta_t : 0.0);
				}
			}
			TriDiagSys(d,e,f,rr,(J+1)*(M+1));
			
			
			#pragma omp parallel for
			for (size_t j=0;j<J+1;j++) {
				for (size_t m=0;m<M+1;m++) {
					size_t lp = j*(M+1)+m;
					rrOld[j][m] = rr[lp];
				}
			}
			
			
			/*
			double delta_t = deltat[t];
			
			for (size_t j=0;j<J+1;j++) {
				for (size_t m=0;m<M+1;m++) {
					if (delta_t > delta_g[m]/abs(Cjm[j][m]) && abs(Cjm[j][m]) > 0.0) delta_t = delta_g[m]/abs(Cjm[j][m]); 
				}
				delta_t = max(delta_t, 1e-5);
				for (size_t m=0;m<M+1;m++) {
					if (m > 1 && m < M-1) {
						
						double b_m_minusHalf = - Cjm_minusHalf[j][m];
						double b_m_plusHalf = - Cjm_plusHalf[j][m];
						double theta_m_minusHalf = (b_m_minusHalf >= 0.0 ? 1.0 : -1.0);
						double theta_m_plusHalf = (b_m_plusHalf >= 0.0 ? 1.0 : -1.0);
						
						double rm_minusHalf = (rrOld[j][m]-rrOld[j][m-1] != 0.0 ?
												( b_m_minusHalf >= 0.0 ?
												(rrOld[j][m-1]-rrOld[j][m-2])/(rrOld[j][m]-rrOld[j][m-1]) :
												(rrOld[j][m+1]-rrOld[j][m])/(rrOld[j][m]-rrOld[j][m-1]) ) : 0.0);
						double rm_plusHalf = (rrOld[j][m+1]-rrOld[j][m] != 0.0 ?
												( b_m_plusHalf >= 0.0 ?
												(rrOld[j][m]-rrOld[j][m-1])/(rrOld[j][m+1]-rrOld[j][m]) :
												(rrOld[j][m+2]-rrOld[j][m+1])/(rrOld[j][m+1]-rrOld[j][m]) ) : 0.0);
						double r = rm_minusHalf;
						double phiMC_m_minusHalf = max(0.0, min(2.0, min(2.0*r,(1.0+r)/2.0)));
						phiMC_m_minusHalf = max(0.0,max(min(1.0,2.0*r),min(2.0,r)));
						r = rm_plusHalf;
						double phiMC_m_plusHalf = max(0.0, min(2.0, min(2.0*r,(1.0+r)/2.0)));
						phiMC_m_plusHalf = max(0.0,max(min(1.0,2.0*r),min(2.0,r)));
						
						double flux_m_minusHalf = 0.5 * b_m_minusHalf * 
								( (1.0+theta_m_minusHalf)*(m > 0 ? rrOld[j][m-1] : 0.0) + (1.0-theta_m_minusHalf)*rrOld[j][m] )
								+ 0.5 * abs(b_m_minusHalf) * (1.0 - abs(b_m_minusHalf*delta_t/delta_g_m_minushalf[m]))
								* phiMC_m_minusHalf * (rrOld[j][m] - (m > 0 ? rrOld[j][m-1] : 0.0));
						double flux_m_plusHalf = 0.5 * b_m_plusHalf * 
								( (1.0+theta_m_plusHalf)*rrOld[j][m] + (1.0-theta_m_plusHalf)*(m < M ? rrOld[j][m+1] : 0.0) )
								+ 0.5 * abs(b_m_plusHalf) * (1.0 - abs(b_m_plusHalf*delta_t/delta_g_m_plushalf[m]))
								* phiMC_m_plusHalf * ((m < M ? rrOld[j][m+1] : 0.0) - rrOld[j][m]);
						// SI NO CAMBIAR EL DELTA G EN LOS FLUJOS
						rrOld[j][m] = rrOld[j][m] + delta_t/delta_g[m] * (flux_m_minusHalf-flux_m_plusHalf)
							+ (j > 0 && j < J ? 0.5*Qjm[j][m] * delta_t : 0.0);
							
						rrOld[j][m] = rrOld[j][m] - delta_t/delta_g[m] * (m < M ? 
										(Cjm[j][m+1]*rrOld[j][m+1]-Cjm[j][m]*rrOld[j][m]) : 0.0)
										+ 0.5*delta_t*Qjm[j][m];
				
					}
				}
			}
			*/
			
			Vector a((J+1)*(M+1),0.0);
			Vector b((J+1)*(M+1),0.0);
			Vector c((J+1)*(M+1),0.0);
			
			// BLOCK 1: Advances in half a time step the energy evolution operator.
			
			#pragma omp parallel for
			for (size_t m=0;m<M+1;m++) {
				for (size_t j=0;j<J+1;j++) {
					size_t l = m*(J+1)+j;
					if (j > 0 && j < J) {
						a[l] = - delta_t/delta_r[j] * Bj_minusHalf_m[m][j] * Wminus_j_minusHalf[m][j] / delta_r_j_minushalf[j];
						b[l] = 1.0 + delta_t/delta_r[j] * ( Bj_minusHalf_m[m][j] * 	Wplus_j_minusHalf[m][j] / delta_r_j_minushalf[j] +
									Bj_plusHalf_m[m][j] * Wminus_j_plusHalf[m][j] / delta_r_j_plushalf[j] ) + delta_t/Tjm[j][m];
						c[l] = - delta_t/delta_r[j] * Bj_plusHalf_m[m][j] * Wplus_j_plusHalf[m][j] / delta_r_j_plushalf[j];
			} else {
						b[l] = 1.0;
					}
					rr[l] = rrOld[j][m];
				}
			}
			TriDiagSys(a,b,c,rr,(J+1)*(M+1));
			
			#pragma omp parallel for
			for (size_t m=0;m<M+1;m++) {
				for (size_t j=0;j<J+1;j++) {
					size_t l = m*(J+1)+j;
					rrOld[j][m] = rr[l];
				}
			}
		}
		if (t % 10 == 0 || t == Ntime-2) {
			matrixInit(dist,J+1,M+1,0.0);
			ofstream fileDist;
			fileDist.open("dist"+to_string(t)+".txt");
			for (size_t j=0;j<J+1;j++) {
				for (size_t m=0;m<M+1;m++) {
					dist[j][m] = rrOld[j][m] / P2(r_local[j]*schwRadius) / factor;
					fileDist << r_local[j] << "\t" << -radialVel(r_local[j]*schwRadius)/cLight
								   << "\t" << g_local[m] << "\t" << dist[j][m] << endl;
				}
			}
			fileDist.close();
		}
	}
	
	matrixInit(dist,J+1,M+1,0.0);
	ofstream fileDist;
	fileDist.open("dist_last.txt");
	for (size_t m=0;m<M+1;m++) {
		for (size_t j=0;j<J+1;j++) {
			size_t l = m*(J+1)+j;
			dist[j][m] = rrOld[j][m] / P2(r_local[j]) / factor;
			fileDist << r_local[j] << "\t" << -radialVel(r_local[j]*schwRadius)/cLight
								   << "\t" << g_local[m] << "\t" << dist[j][m] << endl;
		}
	}
	fileDist.close();

	cout << endl << "TRANSPORT EQUATION: Finished time evolution" << endl << endl;
	
	particle.ps.iterate([&] (const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		if (r/schwRadius < r_local.back()) {
			size_t j = 1;
			while (r_local[j]*schwRadius < r) j++;
			particle.ps.iterate([&] (const SpaceIterator& iRE) {
				double g = iRE.val(DIM_E) / (particle.mass*cLight2);
				size_t m = 0;
				while (g_local[m] < g) m++;
				
				double N11 = dist[j-1][m-1];
				double N12 = dist[j-1][m];
				double N1 = 0.0;
				if (N11 > 0.0 && N12 > 0.0) {
					double s1 = safeLog10(N12/N11)/safeLog10(g_local[m]/g_local[m-1]);
					N1 = N11 * pow(g / g_local[m-1], s1);
				} else {
					if (N11 > 0.0)
						N1 = - N11 * (g - g_local[m-1])/(g_local[m]-g_local[m-1]) + N11;
					else if (N12 > 0.0)
						N1 = N12 * (g - g_local[m-1])/(g_local[m]-g_local[m-1]);
					else
						N1 = 0.0;
				}
				double N21 = dist[j][m-1];
				double N22 = dist[j][m];
				double N2 = 0.0;
				if (N21 > 0.0 && N22 > 0.0) {
					double s2 = safeLog10(N22/N21)/safeLog10(g_local[m]/g_local[m-1]);
					N2 = N21 * pow(g / g_local[m-1], s2);
				} else {
					if (N21 > 0.0)
						N2 = - N21 * (g - g_local[m-1])/(g_local[m]-g_local[m-1]) + N21;
					else if (N22 > 0.0)
						N2 = N22 * (g - g_local[m-1])/(g_local[m]-g_local[m-1]);
					else
						N2 = 0.0;
				}
				
				double N = 0.0;
				if (N1 > 0.0 && N2 > 0.0) {
					double s = safeLog10(N2/N1)/safeLog10(r_local[j]/r_local[j-1]);
					N = N1 * pow(r / (r_local[j-1]*schwRadius), s);
				} else {
					if (N1 > 0.0)
						N = - N1 * (r/schwRadius - r_local[j-1])/(r_local[j]-r_local[j-1]) + N1;
					else if (N2 > 0.0)
						N = N2 * (r/schwRadius - r_local[j-1])/(r_local[j]-r_local[j-1]);
					else
						N = 0.0;
				}
				particle.distribution.set(iRE, N / (particle.mass*cLight2));
			},{-1,iR.coord[DIM_R],0});
		}
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





void distributionSpatialDiffusionSteady(Particle& particle, State& st)
// This routine solves the complete time-dependent transport equation considering turbulence
// acceleration and spatial diffusion.
{
	size_t J = 120;
	double r_min = sqrt(st.denf_e.ps[DIM_R][0] / schwRadius);
	double r_max = (st.denf_e.ps[DIM_R][nR-1] / schwRadius) * paso_r;
	double paso_r_local = pow(r_max/r_min,1.0/J);
	
	// We define a mesh of points for the radial coordinate
	
	Vector r_local(J+1,r_min);
	Vector r_local_extended(J+3,r_min/paso_r_local);
	Vector delta_r_j_plushalf(J+1,0.0);
	Vector delta_r_j_minushalf(J+1,0.0);
	Vector delta_r(J+1,0.0);
		
	for (size_t j=1;j<J+1;j++)	r_local[j] = r_local[j-1]*paso_r_local;
	for (size_t j=1;j<J+3;j++)	r_local_extended[j] = r_local_extended[j-1]*paso_r_local;
	for (size_t j=0;j<J+1;j++)	delta_r_j_plushalf[j] = r_local[j]*(paso_r_local-1.0);
	for (size_t j=0;j<J+1;j++)	delta_r_j_minushalf[j] = r_local[j]*(1.0-1.0/paso_r_local);
	for (size_t j=0;j<J+1;j++)	delta_r[j] = 0.5*r_local[j]*(paso_r_local-1.0/paso_r_local);
	
	// We define a mesh of points for the Lorentz factor coordinate
	size_t M = 150;
	double g_min = 0.9 * (particle.emin() / (particle.mass*cLight2));
	double g_max = 1.1 * (particle.emax() / (particle.mass*cLight2));
	
	double paso_g_local = pow(g_max/g_min,1.0/M);
	Vector g_local(M+1,g_min);
	Vector g_local_extended(M+3,g_min/paso_g_local);
	Vector delta_g_m_plushalf(M+1,0.0);
	Vector delta_g_m_minushalf(M+1,0.0);
	Vector delta_g(M+1,0.0);
		
	for (size_t m=1;m<M+1;m++)	g_local[m] = g_local[m-1]*paso_g_local;
	for (size_t m=1;m<M+3;m++)	g_local_extended[m] = g_local_extended[m-1]*paso_g_local;
	for (size_t m=0;m<M+1;m++)	delta_g_m_plushalf[m] = g_local[m]*(paso_g_local-1.0);
	for (size_t m=0;m<M+1;m++)	delta_g_m_minushalf[m] = g_local[m]*(1.0-1.0/paso_g_local);
	for (size_t m=0;m<M+1;m++)	delta_g[m] = 0.5*g_local[m]*(paso_g_local-1.0/paso_g_local);
	
	// We define a vector for the time evolution
	
	double naturalTimescale = schwRadius/cLight;
	
	fun2 Afun = [&particle,naturalTimescale] (double r, double g) 
				{
					double rTrue = r * schwRadius;
					double height = height_fun(rTrue);
					double B = magneticField(rTrue);
					double vR = radialVel(rTrue);
					double k_rr = diffCoeff_r(g,particle,height,B);
					return - (vR/schwRadius) * naturalTimescale * ( 2.0*k_rr/(rTrue*vR) + 1.0 );
				};

	fun2 Bfun = [&particle,naturalTimescale] (double r, double g)
				{
					double rTrue = r * schwRadius;
					double height = height_fun(rTrue);
					double B = magneticField(rTrue);
					double k_rr = diffCoeff_r(g,particle,height,B);
					return k_rr / (schwRadius*schwRadius) * naturalTimescale;
				};

	fun2 Cfun = [&particle,&st,paso_r_local,naturalTimescale] (double r, double g) 
				{
					double rTrue = r * schwRadius;
					double E = g*particle.mass*cLight2;
					
					size_t jjR=0;
					if (rTrue > particle.ps[DIM_R][0]) {
						if (rTrue < particle.ps[DIM_R][nR-1]) {
							while (particle.ps[DIM_R][jjR] < rTrue) jjR++;
						} else
							jjR = nR-1;
					}
					
					double dgdt = 0.0;
					if (jjR > 0 && jjR < nR-1) {
						SpaceCoord jR1 = {0,jjR-1,0};
						SpaceCoord jR2 = {0,jjR,0};
						double r1 = st.denf_e.ps[DIM_R][jjR-1];
						double r2 = st.denf_e.ps[DIM_R][jjR];
						double dEdt_1 = losses(E,particle,st,jR1);
						double dEdt_2 = losses(E,particle,st,jR2);
						double slope = safeLog10(dEdt_2/dEdt_1)/safeLog10(r2/r1);
						dgdt = - dEdt_1 * pow(rTrue/r1, slope) / (particle.mass*cLight2);
						
					} else if (jjR == 0) {
						SpaceCoord jR = {0,0,0};
						dgdt = - losses(E,particle,st,jR) / (particle.mass*cLight2);
					} else if (jjR == nR-1) {
						SpaceCoord jR = {0,nR-1,0};
						dgdt = - losses(E,particle,st,jR) / (particle.mass*cLight2);
					}
					double vR = radialVel(rTrue);
					double dvRdr = (radialVel(rTrue*paso_r_local) - radialVel(rTrue/paso_r_local)) / 
									(rTrue*(paso_r_local-1.0/paso_r_local));
					
					dvRdr = (dvRdr > 0.0) ? dvRdr : 0.0;
					return - dgdt + 1.0/3.0 * (2.0*vR/rTrue + dvRdr) * g * naturalTimescale;
				};
	
	fun2 Tfun = [&particle,naturalTimescale,&st] (double r, double g) 
				{
					double E = g * particle.mass*cLight2;
					double rTrue = r * schwRadius;
					double height = height_fun(rTrue);
					double B = magneticField(rTrue);
					double vR = radialVel(rTrue);
					//double rateDiff = diffCoeff_r(g,particle,height,B) / (height*height);
					double rateDiff = 1.0/diffusionTimeTurbulence(E,height,particle,B);
					double rateWinds = s*abs(vR)/rTrue;
					
					size_t jjR=0;
					if (rTrue > particle.ps[DIM_R][0]) {
						if (rTrue < particle.ps[DIM_R][nR-1]) {
							while (particle.ps[DIM_R][jjR] < rTrue) jjR++;
						} else
							jjR = nR-1;
					}
					
					double dgdt = 0.0;
					if (jjR > 0 && jjR < nR-1) {
						SpaceCoord jR1 = {0,jjR-1,0};
						SpaceCoord jR2 = {0,jjR,0};
						double r1 = st.denf_e.ps[DIM_R][jjR-1];
						double r2 = st.denf_e.ps[DIM_R][jjR];
						double dEdt_1 = losses(E,particle,st,jR1);
						double dEdt_2 = losses(E,particle,st,jR2);
						double slope = safeLog10(dEdt_2/dEdt_1)/safeLog10(r2/r1);
						dgdt = - dEdt_1 * pow(rTrue/r1, slope) / (particle.mass*cLight2);
						
					} else if (jjR == 0) {
						SpaceCoord jR = {0,0,0};
						dgdt = - losses(E,particle,st,jR) / (particle.mass*cLight2);
					} else if (jjR == nR-1) {
						SpaceCoord jR = {0,nR-1,0};
						dgdt = - losses(E,particle,st,jR) / (particle.mass*cLight2);
					}
					double rateCool = abs(dgdt)/g;
					return pow(rateDiff+rateWinds+rateCool,-1) / naturalTimescale;
				};
		
	fun2 Qfun = [&st,&particle,paso_g_local,naturalTimescale,r_local,J] (double r, double g)
				{
					SpaceCoord i = {0,0,0};
					double rTrue = r * schwRadius;
					double E = g * particle.mass*cLight2;
					double Q = (E > particle.emin() && E < particle.emax() &&
								rTrue > particle.ps[DIM_R].first() && rTrue < particle.ps[DIM_R].last()) ?
								particle.injection.interpolate({{DIM_E,E},{DIM_R,rTrue}},&i) : 0.0;
					double smoother = 0.25 * (1.0 + gsl_sf_erf((r_local[J-5]-r))) 
											* (1.0 + gsl_sf_erf(r-r_local[0]));
					return Q * r * r * naturalTimescale * (particle.mass*cLight2) * smoother;
				};
	
	Matrix Bj_minusHalf_m;				matrixInit(Bj_minusHalf_m,M+1,J+1,0.0);
	Matrix Bj_plusHalf_m;				matrixInit(Bj_plusHalf_m,M+1,J+1,0.0);
	Matrix Cjm_plusHalf;				matrixInit(Cjm_plusHalf,J+1,M+1,0.0);
	Matrix Cjm_minusHalf;				matrixInit(Cjm_minusHalf,J+1,M+1,0.0);
	Matrix Djm_minusHalf;				matrixInit(Djm_minusHalf,J+1,M+1,0.0);
	Matrix Djm_plusHalf;				matrixInit(Djm_plusHalf,J+1,M+1,0.0);
	Matrix Wminus_j_minusHalf;			matrixInit(Wminus_j_minusHalf,M+1,J+1,0.0);
	Matrix Wminus_j_plusHalf;			matrixInit(Wminus_j_plusHalf,M+1,J+1,0.0);
	Matrix Wplus_j_minusHalf;			matrixInit(Wplus_j_minusHalf,M+1,J+1,0.0);
	Matrix Wplus_j_plusHalf;			matrixInit(Wplus_j_plusHalf,M+1,J+1,0.0);
	Matrix Vminus_m_minusHalf;			matrixInit(Vminus_m_minusHalf,J+1,M+1,0.0);
	Matrix Vminus_m_plusHalf;			matrixInit(Vminus_m_plusHalf,J+1,M+1,0.0);
	Matrix Vplus_m_minusHalf;			matrixInit(Vplus_m_minusHalf,J+1,M+1,0.0);
	Matrix Vplus_m_plusHalf;			matrixInit(Vplus_m_plusHalf,J+1,M+1,0.0);
	Matrix Tjm;							matrixInit(Tjm,J+1,M+1,0.0);
	Matrix Qjm;							matrixInit(Qjm,J+1,M+1,0.0);
	
	cout << "TRANSPORT EQUATION: Starting matrix element calculation" << endl << endl;
	
	#pragma omp parallel for
	for (size_t j=0;j<J+1;j++) {
		for (size_t m=0;m<M+1;m++) {
			size_t l = j*(M+1)+m;
			double Ajm = Afun(r_local[j],g_local[m]);
			double Aj_minusHalf_m = 0.5 * (Afun(r_local_extended[j],g_local[m]) + Ajm);
			double Aj_plusHalf_m = 0.5 * (Afun(r_local_extended[j+2],g_local[m]) + Ajm);
					
			double Bjm = Bfun(r_local[j],g_local[m]);
			Bj_minusHalf_m[m][j] = 0.5 * (Bfun(r_local_extended[j],g_local[m]) + Bjm);
			Bj_plusHalf_m[m][j] = 0.5 * (Bfun(r_local_extended[j+2],g_local[m]) + Bjm);
			
			double Cjm = Cfun(r_local[j],g_local[m]);
			Cjm_minusHalf[j][m] = 0.5 * (Cfun(r_local[j],g_local_extended[m]) + Cjm);
			Cjm_plusHalf[j][m] = 0.5 * (Cfun(r_local[j],g_local_extended[m+2]) + Cjm);
			
			Qjm[j][m] = Qfun(r_local[j],g_local[m]);
			Tjm[j][m] = Tfun(r_local[j],g_local[m]);
					
			double wj_minusHalf = Aj_minusHalf_m/Bj_minusHalf_m[m][j] * delta_r_j_minushalf[j];
			double wj_plusHalf = Aj_plusHalf_m/Bj_plusHalf_m[m][j] * delta_r_j_plushalf[j];
			double Wj_minusHalf(0.0), Wj_plusHalf(0.0);
			
			if ( abs(wj_minusHalf) < 1.0e-3 ) {
				Wj_minusHalf = pow(1.0+gsl_pow_2(wj_minusHalf)/24+gsl_pow_4(wj_minusHalf)/1920,-1);
				Wplus_j_minusHalf[m][j] = Wj_minusHalf * exp(0.5*wj_minusHalf);
				Wminus_j_minusHalf[m][j] = Wj_minusHalf * exp(-0.5*wj_minusHalf);
			} else {
				Wj_minusHalf = abs(wj_minusHalf)*exp(-0.5*abs(wj_minusHalf)) /
								(1.0-exp(-abs(wj_minusHalf)));
				if (wj_minusHalf > 0.0) {
					Wplus_j_minusHalf[m][j] = abs(wj_minusHalf) / (1.0-exp(-abs(wj_minusHalf)));
					Wminus_j_minusHalf[m][j] = Wj_minusHalf * exp(-0.5*wj_minusHalf);
				} else {
					Wplus_j_minusHalf[m][j] = Wj_minusHalf * exp(0.5*wj_minusHalf);
					Wminus_j_minusHalf[m][j] = abs(wj_minusHalf) / (1.0-exp(-abs(wj_minusHalf)));
				}
			}
			if ( abs(wj_plusHalf) < 1.0e-3 ) {
				Wj_plusHalf = pow(1.0+gsl_pow_2(wj_plusHalf)/24+gsl_pow_4(wj_plusHalf)/1920,-1);
				Wplus_j_plusHalf[m][j] = Wj_plusHalf * exp(0.5*wj_plusHalf);
				Wminus_j_plusHalf[m][j] = Wj_plusHalf * exp(-0.5*wj_plusHalf);
			} else {
				Wj_plusHalf = abs(wj_plusHalf)*exp(-0.5*abs(wj_plusHalf)) /
								(1.0-exp(-abs(wj_plusHalf)));
				if (wj_plusHalf > 0.0) {
					Wplus_j_plusHalf[m][j] = abs(wj_plusHalf) / (1.0-exp(-abs(wj_plusHalf)));
					Wminus_j_plusHalf[m][j] = Wj_plusHalf * exp(-0.5*wj_plusHalf);
				} else {
					Wplus_j_plusHalf[m][j] = Wj_plusHalf * exp(0.5*wj_plusHalf);
					Wminus_j_plusHalf[m][j] = abs(wj_plusHalf) / (1.0-exp(-abs(wj_plusHalf)));
				}
			}
		}
	}
	
	ofstream fileMatrixDistEne, fileMatrixDistPos, fileMatrixDistQT;
	fileMatrixDistEne.open("matrixCoeffTransportEq_Ene.txt");
	fileMatrixDistPos.open("matrixCoeffTransportEq_Pos.txt");
	fileMatrixDistQT.open("matrixCoeffTransportEq_QT.txt");
	fileMatrixDistEne 	<< "j \t m \t Cj- \t Cj+" << endl;
	fileMatrixDistPos 	<< "m \t j \t Bj- \t Bj+ \t W-- \t W-+ \t W+- \t W++" << endl;
	fileMatrixDistQT 	<< "j \t m \t T \t Q \t" << endl;
	
	for (size_t j=0;j<J+1;j++) {
		for (size_t m=0;m<M+1;m++) {
			fileMatrixDistEne 	<< j << "\t" << m << "\t"
								<< Cjm_minusHalf[j][m] << "\t" << Cjm_plusHalf[j][m] << endl;
			fileMatrixDistQT	<< j << "\t" << m << "\t"
								<< Tjm[j][m] << "\t" << Qjm[j][m] << endl;
		}
	}
	for (size_t m=0;m<M+1;m++) {
		for (size_t j=0;j<J+1;j++) {
			fileMatrixDistPos	<< m << "\t" << j << "\t"
								<< Bj_minusHalf_m[m][j] << "\t" << Bj_plusHalf_m[m][j] << "\t"
								<< Wminus_j_minusHalf[m][j] << "\t" << Wminus_j_plusHalf[m][j] << "\t"
								<< Wplus_j_minusHalf[m][j] << "\t" << Wplus_j_plusHalf[m][j] << endl;
		}
	}
 
	fileMatrixDistEne.close();
	fileMatrixDistPos.close();
	fileMatrixDistQT.close();
	
	cout << "TRANSPORT EQUATION: Finished matrix calculation" << endl << endl;
	
	Vector rr((J+1)*(M+1),0.0);
	Matrix rrOld;
	matrixInit(rrOld,J+1,M+1,0.0);
	Matrix dist;
	
	Vector a((J+1)*(M+1),0.0);
	Vector b((J+1)*(M+1),0.0);
	Vector c((J+1)*(M+1),0.0);
	Vector d((J+1)*(M+1),0.0);
	Vector e((J+1)*(M+1),0.0);
	Vector f((J+1)*(M+1),0.0);
	
	
	#pragma omp parallel for
	for (size_t m=0;m<M+1;m++) {
		for (size_t j=0;j<J+1;j++) {
			size_t l = m*(J+1)+j;
			if (j > 0 && j < J) {
				if (m > 0)
					a[l] = - 1.0/delta_r[j] * Bj_minusHalf_m[m][j] * Wminus_j_minusHalf[m][j] / delta_r_j_minushalf[j];
				b[l] = 1.0/delta_r[j] * ( Bj_minusHalf_m[m][j] * Wplus_j_minusHalf[m][j] / delta_r_j_minushalf[j] +
							Bj_plusHalf_m[m][j] * Wminus_j_plusHalf[m][j] / delta_r_j_plushalf[j] ) + 1.0/Tjm[j][m];
				if (m < M)
					c[l] = - 1.0/delta_r[j] * Bj_plusHalf_m[m][j] * Wplus_j_plusHalf[m][j] / delta_r_j_plushalf[j];
			} else {
				b[l] = 1.0;
			}
			
			rr[l] = (j > 0 && j < J ? Qjm[j][m] : 0.0);
		}
	}
	TriDiagSys(a,b,c,rr,(J+1)*(M+1));
	
	#pragma omp parallel for
	for (size_t m=0;m<M+1;m++) {
		for (size_t j=0;j<J+1;j++) {
			size_t l = m*(J+1)+j;
			rrOld[j][m] = rr[l];
		}
	}
	
	matrixInit(dist,J+1,M+1,0.0);
	ofstream fileDist;
	fileDist.open("dist_last.txt");
	for (size_t j=0;j<J+1;j++) {
		for (size_t m=0;m<M+1;m++) {
			size_t l = j*(M+1)+m;
			dist[j][m] = rrOld[j][m] / P2(r_local[j]);
			fileDist << r_local[j] << "\t" << -radialVel(r_local[j]*schwRadius)/cLight
								   << "\t" << g_local[m] << "\t" << dist[j][m] << endl;
		}
	}
	fileDist.close();

	cout << endl << "TRANSPORT EQUATION: Finished time evolution" << endl << endl;
	
	particle.ps.iterate([&] (const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		if (r/schwRadius < r_local.back()) {
			size_t j = 1;
			while (r_local[j]*schwRadius < r) j++;
			particle.ps.iterate([&] (const SpaceIterator& iRE) {
				double g = iRE.val(DIM_E) / (particle.mass*cLight2);
				size_t m = 0;
				while (g_local[m] < g) m++;
				
				double N11 = dist[j-1][m-1];
				double N12 = dist[j-1][m];
				double N1 = 0.0;
				if (N11 > 0.0 && N12 > 0.0) {
					double s1 = safeLog10(N12/N11)/safeLog10(g_local[m]/g_local[m-1]);
					N1 = N11 * pow(g / g_local[m-1], s1);
				} else {
					if (N11 > 0.0)
						N1 = - N11 * (g - g_local[m-1])/(g_local[m]-g_local[m-1]) + N11;
					else if (N12 > 0.0)
						N1 = N12 * (g - g_local[m-1])/(g_local[m]-g_local[m-1]);
					else
						N1 = 0.0;
				}
				double N21 = dist[j][m-1];
				double N22 = dist[j][m];
				double N2 = 0.0;
				if (N21 > 0.0 && N22 > 0.0) {
					double s2 = safeLog10(N22/N21)/safeLog10(g_local[m]/g_local[m-1]);
					N2 = N21 * pow(g / g_local[m-1], s2);
				} else {
					if (N21 > 0.0)
						N2 = - N21 * (g - g_local[m-1])/(g_local[m]-g_local[m-1]) + N21;
					else if (N22 > 0.0)
						N2 = N22 * (g - g_local[m-1])/(g_local[m]-g_local[m-1]);
					else
						N2 = 0.0;
				}
				
				double N = 0.0;
				if (N1 > 0.0 && N2 > 0.0) {
					double s = safeLog10(N2/N1)/safeLog10(r_local[j]/r_local[j-1]);
					N = N1 * pow(r / (r_local[j-1]*schwRadius), s);
				} else {
					if (N1 > 0.0)
						N = - N1 * (r/schwRadius - r_local[j-1])/(r_local[j]-r_local[j-1]) + N1;
					else if (N2 > 0.0)
						N = N2 * (r/schwRadius - r_local[j-1])/(r_local[j]-r_local[j-1]);
					else
						N = 0.0;
				}
				particle.distribution.set(iRE, N / (particle.mass*cLight2));
			},{-1,iR.coord[DIM_R],0});
		}
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





void distributionSpatialDiffusionSteady2(Particle& particle, State& st)
// This routine solves the complete time-dependent transport equation considering turbulence
// acceleration and spatial diffusion.
{
	double paso_E = pow(particle.emax()/particle.emin(),1.0/(particle.ps[DIM_E].size()-1));
	double restEnergy = particle.mass * cLight2;
	particle.ps.iterate([&](const SpaceIterator& iE) {
		
		SpaceCoord jE = {particle.ps[DIM_E].size()-1-iE.coord[DIM_E],0,0};
		double E = particle.ps[DIM_E][jE[DIM_E]];
		double g = E / restEnergy;
	
		size_t J = 120;
		double r_min = sqrt(st.denf_e.ps[DIM_R][0] / schwRadius);
		double r_max = (st.denf_e.ps[DIM_R][nR-1] / schwRadius) * paso_r;
		double paso_r_local = pow(r_max/r_min,1.0/J);
		
		// We define a mesh of points for the radial coordinate
		
		Vector r_local(J+1,r_min);
		Vector r_local_extended(J+3,r_min/paso_r_local);
		Vector delta_r_j_plushalf(J+1,0.0);
		Vector delta_r_j_minushalf(J+1,0.0);
		Vector delta_r(J+1,0.0);
			
		for (size_t j=1;j<J+1;j++)	r_local[j] = r_local[j-1]*paso_r_local;
		for (size_t j=1;j<J+3;j++)	r_local_extended[j] = r_local_extended[j-1]*paso_r_local;
		for (size_t j=0;j<J+1;j++)	delta_r_j_plushalf[j] = r_local[j]*(paso_r_local-1.0);
		for (size_t j=0;j<J+1;j++)	delta_r_j_minushalf[j] = r_local[j]*(1.0-1.0/paso_r_local);
		for (size_t j=0;j<J+1;j++)	delta_r[j] = 0.5*r_local[j]*(paso_r_local-1.0/paso_r_local);
	
	
		double naturalTimescale = schwRadius/cLight;
	
		fun1 Afun = [&particle,naturalTimescale,g] (double r) 
					{
						double rTrue = r * schwRadius;
						double height = height_fun(rTrue);
						double B = magneticField(rTrue);
						double vR = radialVel(rTrue);
						double k_rr = diffCoeff_r(g,particle,height,B);
						return - (vR/schwRadius) * naturalTimescale * ( 2.0*k_rr/(rTrue*vR) + 1.0 );
					};

		fun1 Bfun = [&particle,naturalTimescale,g] (double r)
					{
						double rTrue = r * schwRadius;
						double height = height_fun(rTrue);
						double B = magneticField(rTrue);
						double k_rr = diffCoeff_r(g,particle,height,B);
						return k_rr / (schwRadius*schwRadius) * naturalTimescale;
					};
		
		fun1 Tfun = [&particle,naturalTimescale,&st,g,paso_E] (double r) 
					{
						double E = g * particle.mass*cLight2;
						double rTrue = r * schwRadius;
						double height = height_fun(rTrue);
						double B = magneticField(rTrue);
						double vR = radialVel(rTrue);
						double rateDiff = diffCoeff_r(g,particle,height,B) / (height*height);
						double rateWinds = s*abs(vR)/rTrue;
						
						size_t jjR=0;
						if (rTrue > particle.ps[DIM_R][0]) {
							if (rTrue < particle.ps[DIM_R][nR-1]) {
								while (particle.ps[DIM_R][jjR] < rTrue) jjR++;
							} else
								jjR = nR-1;
						}
						
						SpaceCoord jjEaux = {0,jjR,0};
						double loss = losses(E,particle,st,jjEaux);
						double rateCool = loss/(E*(sqrt(paso_E)-1.0/sqrt(paso_E)));
						return pow(rateDiff+rateWinds,-1) / naturalTimescale;
					};
			
		fun1 Qfun = [&st,&particle,naturalTimescale,g,jE,paso_E] (double r)
					{
						double rTrue = r * schwRadius;
						double E = g*particle.mass*cLight2;
			
						size_t jjR=0;
						particle.ps.iterate([&](const SpaceIterator& iR) {
							if (iR.val(DIM_R) < rTrue) jjR++;
						},{0,-1,0});
						SpaceCoord jEaux = {jE[DIM_E],jjR,0};
						
						double Qplus = 0.0;
						double Nplus = 0.0;
						double loss = 0.0;
						if (jE[DIM_E] < particle.ps[DIM_E].size()-1) {
							SpaceCoord jEplus = {jE[DIM_E]+1,jjR,0};
							Nplus = (rTrue > particle.ps[DIM_R].first() &&
											rTrue < particle.ps[DIM_R].last()) ?
								particle.distribution.interpolate({{DIM_R,rTrue}},&jEplus)*r*r*particle.mass*cLight2 : 0.0;
							loss = losses(st.denf_e.ps[DIM_E][jE[DIM_E]],particle,st,jEplus);
							Qplus = Nplus * abs(loss) / (E*(sqrt(paso_E)-1.0/sqrt(paso_E)));
						}
		
						double Qlocal = (rTrue > particle.ps[DIM_R].first() && rTrue < particle.ps[DIM_R].last()) ?
									particle.injection.interpolate({{DIM_R,rTrue}},&jEaux)*r*r*particle.mass*cLight2 : 0.0;
						cout << Nplus << "\t" << abs(loss) / (E*(sqrt(paso_E)-1.0/sqrt(paso_E))) << "\t" << Qlocal << endl;
						return naturalTimescale*(Qlocal+Qplus);
					};
	
		Vector Bj_minusHalf(J+1,0.0);
		Vector Bj_plusHalf(J+1,0.0);
		Vector Wminus_j_minusHalf(J+1,0.0);
		Vector Wminus_j_plusHalf(J+1,0.0);
		Vector Wplus_j_minusHalf(J+1,0.0);
		Vector Wplus_j_plusHalf(J+1,0.0);
		Vector Tj(J+1,0.0);
		Vector Qj(J+1,0.0);
		
		cout << "TRANSPORT EQUATION: Starting matrix element calculation" << endl << endl;
		
		#pragma omp parallel for
		for (size_t j=0;j<J+1;j++) {
			double Aj = Afun(r_local[j]);
			double Aj_minusHalf = 0.5 * (Afun(r_local_extended[j]) + Aj);
			double Aj_plusHalf = 0.5 * (Afun(r_local_extended[j+2]) + Aj);
			
			double Bj = Bfun(r_local[j]);
			Bj_minusHalf[j] = 0.5 * (Bfun(r_local_extended[j]) + Bj);
			Bj_plusHalf[j] = 0.5 * (Bfun(r_local_extended[j+2]) + Bj);
			
			Qj[j] = Qfun(r_local[j]);
			Tj[j] = Tfun(r_local[j]);
			
			double wj_minusHalf = Aj_minusHalf/Bj_minusHalf[j] * delta_r_j_minushalf[j];
			double wj_plusHalf = Aj_plusHalf/Bj_plusHalf[j] * delta_r_j_plushalf[j];
			double Wj_minusHalf(0.0), Wj_plusHalf(0.0);
			
			if ( abs(wj_minusHalf) < 1.0e-3 ) {
				Wj_minusHalf = pow(1.0+gsl_pow_2(wj_minusHalf)/24+gsl_pow_4(wj_minusHalf)/1920,-1);
				Wplus_j_minusHalf[j] = Wj_minusHalf * exp(0.5*wj_minusHalf);
				Wminus_j_minusHalf[j] = Wj_minusHalf * exp(-0.5*wj_minusHalf);
			} else {
				Wj_minusHalf = abs(wj_minusHalf)*exp(-0.5*abs(wj_minusHalf)) /
								(1.0-exp(-abs(wj_minusHalf)));
				if (wj_minusHalf > 0.0) {
					Wplus_j_minusHalf[j] = abs(wj_minusHalf) / (1.0-exp(-abs(wj_minusHalf)));
					Wminus_j_minusHalf[j] = Wj_minusHalf * exp(-0.5*wj_minusHalf);
				} else {
					Wplus_j_minusHalf[j] = Wj_minusHalf * exp(0.5*wj_minusHalf);
					Wminus_j_minusHalf[j] = abs(wj_minusHalf) / (1.0-exp(-abs(wj_minusHalf)));
				}
			}
			if ( abs(wj_plusHalf) < 1.0e-3 ) {
				Wj_plusHalf = pow(1.0+gsl_pow_2(wj_plusHalf)/24+gsl_pow_4(wj_plusHalf)/1920,-1);
				Wplus_j_plusHalf[j] = Wj_plusHalf * exp(0.5*wj_plusHalf);
				Wminus_j_plusHalf[j] = Wj_plusHalf * exp(-0.5*wj_plusHalf);
			} else {
				Wj_plusHalf = abs(wj_plusHalf)*exp(-0.5*abs(wj_plusHalf)) /
								(1.0-exp(-abs(wj_plusHalf)));
				if (wj_plusHalf > 0.0) {
					Wplus_j_plusHalf[j] = abs(wj_plusHalf) / (1.0-exp(-abs(wj_plusHalf)));
					Wminus_j_plusHalf[j] = Wj_plusHalf * exp(-0.5*wj_plusHalf);
				} else {
					Wplus_j_plusHalf[j] = Wj_plusHalf * exp(0.5*wj_plusHalf);
					Wminus_j_plusHalf[j] = abs(wj_plusHalf) / (1.0-exp(-abs(wj_plusHalf)));
				}
			}
		}
		cout << "TRANSPORT EQUATION: Finished matrix calculation" << endl << endl;
		
		Vector rr(J+1,0.0);
		
		Vector a((J+1),0.0);
		Vector b((J+1),0.0);
		Vector c((J+1),0.0);
		
		#pragma omp parallel for
		for (size_t j=0;j<J+1;j++) {
			if (j > 0 && j < J) {
				if (jE[DIM_E] > 0)
					a[j] = - 1.0/delta_r[j] * Bj_minusHalf[j] * Wminus_j_minusHalf[j] / delta_r_j_minushalf[j];
				b[j] = 1.0/delta_r[j] * ( Bj_minusHalf[j] * 	Wplus_j_minusHalf[j] / delta_r_j_minushalf[j] +
							Bj_plusHalf[j] * Wminus_j_plusHalf[j] / delta_r_j_plushalf[j] ) + 1.0/Tj[j];
				if (jE[DIM_E] < particle.ps[DIM_E].size()-1)
					c[j] = - 1.0/delta_r[j] * Bj_plusHalf[j] * Wplus_j_plusHalf[j] / delta_r_j_plushalf[j];
			} else
				b[j] = 1.0;
			rr[j] = (j > 0 && j < J ? Qj[j] : 0.0);
		}
		TriDiagSys(a,b,c,rr,J+1);
		
		ofstream rrFile;
		rrFile.open("rrAux.txt", std::ofstream::app);
		for (int j=0;j<J+1;j++) {
			rrFile << r_local[j] << "\t" << rr[j] << endl;
		}
		rrFile.close();
		
		particle.ps.iterate([&] (const SpaceIterator& jER) {
			double r = jER.val(DIM_R);
			if (r/schwRadius < r_local.back()) {
				size_t j = 1;
				while (r_local[j]*schwRadius < r) j++;
				
				double N1 = rr[j-1] / P2(r_local[j-1]);
				double N2 = rr[j] / P2(r_local[j]);
				double N = 0.0;
				if (N1 > 0.0 && N2 > 0.0) {
					double ss = safeLog10(N2/N1)/safeLog10(r_local[j]/r_local[j-1]);
					N = N1 * pow(r / (r_local[j-1]*schwRadius), ss);
				} else {
					if (N1 > 0.0)
						N = - N1 * (r/schwRadius - r_local[j-1])/(r_local[j]-r_local[j-1]) + N1;
					else if (N2 > 0.0)
						N = N2 * (r/schwRadius - r_local[j-1])/(r_local[j]-r_local[j-1]);
					else
						N = 0.0;
				}
				particle.distribution.set(jER, N / (particle.mass*cLight2));
				
				double x=3;
				
				x = x+1;
				
			}
		},{jE[DIM_E],-1,0});
	},{-1,0,0});
	
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