#include "NTdistribution.h"
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
					else if (p.id == "ntElectron") {
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
			if (p.id == "ntProton")	writeEandRParamSpace("linearEmitter_p",Nle,0);
			if (p.id == "ntElectron")	writeEandRParamSpace("linearEmitter_e",Nle,0);
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
}


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
		size_t nPoints = 30;
		const double pasoRc = pow(p.ps[DIM_R].last()/r,1.0/(nPoints-1));
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
		cout << "uNTh/uTh = " << sumNT << endl;
	},{0,-1,0});
}

/*
#include "Particles.h"
#include "Radiation.h"
void distributionGAMERA(Particle& p, State& st)
{
	ofstream file;
	file.open("distributionGAMERA.dat");
	p.ps.iterate([&](const SpaceIterator& iR) {
		Vector e(nE,0.0);
		Vector q(nE,0.0);
		Particles *electron = new Particles();
		double tacc = accretionTime(iR.val(DIM_R));
		p.ps.iterate([&](const SpaceIterator& iRE) {
			e[iRE.coord[DIM_E]] = iRE.val(DIM_E);
			q[iRE.coord[DIM_E]] = p.injection.get(iRE);
		},{-1,iR.coord[DIM_R],0});
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