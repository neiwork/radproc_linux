#include "distribution.h"
#include "injection.h"
#include "messages.h"
#include "write.h"
#include "adafFunctions.h"
#include "globalVariables.h"
#include "losses.h"
#include <fmath/RungeKutta.h>
//#include "timeDistribution.h"

#include <iostream>

#include <fparameters/Dimension.h>
#include <fparameters/SpaceIterator.h>
//#include <fparameters/parameters.h>
//#include <fmath/physics.h>
//#include <boost/property_tree/ptree.hpp>


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

int effectiveE_2(Vector& Eeff, Vector& rCurrent, double rEnd, size_t nRR, double Emax, Particle& p,
						State& st, const SpaceCoord& i)
{
	double pasoR = pow(rEnd/rCurrent[0],1.0/nRR);
	size_t j=0;
	while (rCurrent[j] < rEnd && Eeff[j] < Emax) {
		double rAux = rCurrent[j]*sqrt(pasoR);
		double be = -b(Eeff[j],rAux,p,st,i);
		double v = -radialVel(rAux);
		double dr = rCurrent[j]*(pasoR-1.0);
		Eeff[j+1] = Eeff[j] + be/v * dr;
		rCurrent[j+1] = rCurrent[j]*pasoR;
		j++;
	}
	return --j;
}

double nLocal(double dist, Vector Eeff, Vector rCurrent, size_t nRR, Particle& p, State& st,
					const SpaceCoord& i)
{
	for (size_t j=0;j<nRR;j++) {
		size_t k = nRR-1-j;
		double dE,dbde;
		double dr = (rCurrent[nRR-1]/rCurrent[0],1.0/nRR);
		if (k-1 >= 0) {
			dE = Eeff[k]-Eeff[k-1];
			dbde = (b(Eeff[k],rCurrent[k],p,st,i)-b(Eeff[k-1],rCurrent[k-1],p,st,i))/dE;
		} else {
			dE = Eeff[k+1]-Eeff[k];
			dbde = (b(Eeff[k+1],rCurrent[k+1],p,st,i)-b(Eeff[k],rCurrent[k],p,st,i))/dE;
		}
		dist *= (1.0+dbde/(-radialVel(rCurrent[k]))*dr);
	}
	return (dist > 0.0) ? dist : 0.0;
}


void distributionWorking(Particle& p, State& st)
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
				//double ratioLosses = Eeff/E;
				double dist2 = dist*ratioLosses;
				Nle.set(itRRE,dist2);
			},{-1,itRR,0});
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
			if (p.id == "ntProton") writeEandRParamSpace("linearEmitter_p",Nle,0);
			if (p.id == "ntElectron") writeEandRParamSpace("linearEmitter_e",Nle,0);
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


void distribution2(Particle& p, State& st)
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
							return p.injection.interpolate({{DIM_E,e}},&itRR.coord);});
						Nle.set(itRRE,integ/losses(energy,p,st,itRRE)*vol);
					}
				} else
					Nle.set(itRRE,0.0);
			},{-1,itRR.coord[DIM_R],0});
		},{0,-1,0});

		if (p.id == "ntProtton") {
			for (int itRR=itR.coord[DIM_R]-1;itRR>=0;itRR--) {
				double rprim = p.ps[DIM_R][itRR];
				double delta_r = (rprim/sqrt(paso_r))*(paso_r-1.0);
				double vel = -radialVel(rprim);
				double tcell = delta_r/vel;
				p.ps.iterate([&](const SpaceIterator& itRRE) {  // para cada energía
					
					/*
					const double E = itRRE.val(DIM_E);
					double Emax = p.emax();
					double Eeff = Emax;
					if (itRRE.its[0].canPeek(+1)) { 
						Eeff = effectiveE(E,Emax,tcell,p,st,itRRE);
					}
					SpaceCoord itRRE_plus_1 = itRRE.moved({0,+1,0});
					double dist = (Eeff < Emax) ? 
									Nle.interpolate({{DIM_E,Eeff}},&itRRE_plus_1) : 0.0;
					//double ratioLosses = losses(Eeff,p,st,itRRE_plus_1)/losses(Eeff,p,st,itRRE);
					//double ratioLosses = 2.0-losses(E,p,st,itRRE_plus_1)/losses(Eeff,p,st,itRRE_plus_1);
					//double ratioLosses = Eeff/E;
					
					double deltaR = itRRE_plus_1.val(DIM_R) - itRRE.val(DIM_R);
					double dr = deltaR/100.0;
					double rcurrent = itRRE_plus_1.val(DIM_R);
					double ecurrent = Eeff;
					double dist2 = dist;
					while (rcurrent > itRRE.val(DIM_R)){
						double dt = dr/(-radialVel(rcurrent));
						double dbde = losses(ecurrent,p,st,itRRE_plus_1)-losses(ecurrent-de,p,st,itRRE_plus_1);
						dbde /= de;
						dist2 = dist2 * (1.0+dbde*dt);
						ecurrent = ecurrent - losses(ecurrent,p,st,itRRE_plus_1)*dt;
						rcurrent -= dr;
					}
					
					double logb2 = 0.5*(log10(losses(Eeff,p,st,itRRE_plus_1))+log10(losses(Eeff,p,st,itRRE)));
					double b2 = pow(10.0,logb2);
					double logb1 = 0.5*(log10(losses(E,p,st,itRRE_plus_1))+log10(losses(E,p,st,itRRE)));
					double b1 = pow(10.0,logb1);
					double dbde = (b2-b1)/(Eeff-E);
					double ratioLosses = 1.0 + dbde * tcell;
					double dist2 = dist*ratioLosses;
					*/
					double Emax = p.emax();
					double E = itRRE.val(DIM_E);
					size_t nRR = 10;
					Vector Eeff(nRR+2,E);
					Vector rCurrent(nRR+2,rprim);
					SpaceCoord itRRE_plus_1 = itRRE.moved({0,+1,0});
					double rEnd = p.ps[DIM_R][itRR+1];
					effectiveE_2(Eeff,rCurrent,rEnd,nRR,Emax,p,st,itRRE);
					double dist = (Eeff[nRR-1] < Emax) ? 
									Nle.interpolate({{DIM_E,Eeff[nRR-1]}},&itRRE_plus_1) : 0.0;
					double dist2 = nLocal(dist,Eeff,rCurrent,nRR,p,st,itRRE);
					
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
                double totalNumberOfParticles = 0.0;
                double pasoE = pow(p.ps[DIM_E][nE-1]/p.ps[DIM_E][0],1.0/nE);
				
				totalNumberOfParticles = RungeKuttaSimple(p.emin(),p.emax(),[&](double e) {
					return e*Nle.interpolate({{DIM_E,e}},&itRR.coord);});
				p.ps.iterate([&](const SpaceIterator& itRRE) {
                    double dE = itRRE.val(DIM_E)*(pasoE-1.0);
                    //totalNumberOfParticles += Nle.get(itRRE)*dE;
					Nle.set(itRRE,Nle.get(itRRE)/vol);
				},{-1,itRR.coord[DIM_R],0});
                if (p.id == "ntProton")
                    cout << "Total Number of Particles = " << totalNumberOfParticles << endl;
			},{0,-1,0});
			if (p.id == "ntProton") writeEandRParamSpace("linearEmitter_p",Nle,0);
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
	
	if (0) {
        double rMax = p.ps[DIM_R].last();
		double rMin = p.ps[DIM_R].first();
		p.ps.iterate([&](const SpaceIterator& it) {
			double r = it.val(DIM_R);
			double e = it.val(DIM_E);
			size_t nPoints = 50;
			Vector Rc(nPoints,r);
			Vector Ec(nPoints,e);
			double int_Rc = pow(p.ps[DIM_R].last()/r,1.0/nPoints);
			for (size_t j=0;j<nPoints-1;j++) {
				double vr = 0.0;
				double be = 1.0;
				if (r > rMin && r < rMax) {
					vr = radialVel(Rc[j]);
					be = b(Ec[j],Rc[j],p,st,it);
				} else {
					vr = radialVel(Rc[j]*1.1);
					be = b(Ec[j],Rc[j],p,st,it);
				}
				double deriv = be/vr;
				double dRc = Rc[j]*(sqrt(int_Rc)-1.0/sqrt(int_Rc));
				Ec[j+1] = Ec[j] + deriv*dRc;
				Rc[j+1] = Rc[j]*int_Rc;
			}
			
			double N = 0.0;
			
			for (size_t j=0;j<nPoints-1;j++) {
				size_t k = nPoints-1-j;
				if (((Ec[k] > p.emin() && Ec[k+1] < p.emax()) && (Ec[k+1] > p.emin() && Ec[k] < p.emax())) &&
				(Rc[k] > p.ps[DIM_R][0] && Rc[k] < p.ps[DIM_R][nR-1])) {
					double Q = p.injection.interpolate({{DIM_E,Ec[k]},{DIM_R,Rc[k]}},&it.coord);
					double dRc = Rc[k]*(sqrt(int_Rc)-1.0/sqrt(int_Rc));
					double dEc = (Ec[k+1] > Ec[k]) ? Ec[k+1]-Ec[k] : 0.0;
					double dbde = (dEc > 0.0) ? (b(Ec[k+1],Rc[k],p,st,it.coord) - 
									b(Ec[k],Rc[k],p,st,it))/dEc : 0.0;
                                    
					double deriv = (Q-N*dbde) / radialVel(Rc[k]);
					N += deriv * (-dRc);
				}
			}
			p.distribution.set(it,N);
			

			/*for (size_t jE=0;jE<nPoints;jE++) {
				if ((Ec[jE] > p.emin() && Ec[jE] < p.emax()) &&
				(Rc[jE] > p.ps[DIM_R][0] && Rc[jE] < p.ps[DIM_R][nR-1])) {
					double mu = 0.0;
					for (size_t jjE=0;jjE<=jE;jjE++) {
						if ((Ec[jjE] > p.emin() && Ec[jjE] < p.emax()) &&
						(Rc[jjE]/sqrt(int_Rc) > p.ps[DIM_R][0] && 
						Rc[jjE]*sqrt(int_Rc) < p.ps[DIM_R][nR-1])) {
                            double dEE = 0.0;
                            if (jjE+1 < nPoints)
                                dEE = (Ec[jjE+1]-Ec[jjE]) ? Ec[jjE+1]-Ec[jjE] : 0.0;
							double vr = abs(radialVel(Rc[jjE]));
							double dRc = Rc[jjE]*(sqrt(int_Rc)-1.0/sqrt(int_Rc));
							double dbdr = (b(Ec[jjE],Rc[jjE]*sqrt(int_Rc),p,st,it) - 
											b(Ec[jjE],Rc[jjE]/sqrt(int_Rc),p,st,it))/dRc;
							double be = b(Ec[jjE],Rc[jjE],p,st,it);
							mu += vr*dbdr/(be*be)*dEE;
							jjE++;
						}
					}
					mu = exp(mu);
                    double dEc = 0.0;
                    if (jE+1 < nPoints)
                        dEc = (Ec[jE+1] > Ec[jE]) ? Ec[jE+1]-Ec[jE] : 0.0;
					N += p.injection.interpolate({{DIM_E,Ec[jE]},{DIM_R,Rc[jE]}},&it.coord)*mu*dEc;
				}
			}
			N = (e > p.emin() && e < p.emax()) &&
				(r > p.ps[DIM_R][0] && r) ? N/abs(b(e,r,p,st,it)) : 0.0;
			p.distribution.set(it,N);*/
            
		},{-1,-1,0});
	}

	if (p.id == "ntElectron") show_message(msgEnd,Module_electronDistribution);
	else if (p.id == "ntProton") show_message(msgEnd,Module_protonDistribution);
}

int rangeE(double e, Particle& p)
{
	return (e > p.emin() && e < p.emax());
}

int rangeR(double r, Particle& p)
{
	return (r > p.ps[DIM_R].first() && r < p.ps[DIM_R].last());
}

double fAux(double r, double pasoR)
{
	double cos1 = costhetaH(r);
	double cos2 = costhetaH(r*pasoR);
	double cos0 = costhetaH(r*sqrt(pasoR));
	double dcos = cos2-cos1;
	double d2cos = (cos2+cos1-2.0*cos0)/(paso_r-1.0);
	return -(2.0*(pasoR-1.0) + 4.0/3.0 * dcos) / (cos0 + 1.0/3.0 * dcos/ (paso_r-1.0)) - 
				d2cos/(cos1+1.0/3.0 * dcos);
}


void distribution3(Particle& p, State& st)
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
			/*
			for (size_t j=0;j<nPoints;j++) {
				double Qinj = (rangeE(Ec[j],p) && rangeR(Rc[j],p)) ? 
						p.injection.interpolate({{DIM_E,Ec[j]},{DIM_R,Rc[j]}},&itRE.coord) : 0.0;
				double mu = 0.0;
				double pasoEaux = 1.000001;
				for (size_t jj=0;jj<j;jj++) {
					double dEe = Ec[jj]*(sqrt(pasoEaux)-1.0/sqrt(pasoEaux));
					double dbde = ( (rangeE(Ec[jj]*sqrt(pasoEaux),p) &&
							rangeE(Ec[jj]/sqrt(pasoEaux),p)) && rangeR(Rc[jj],p) ) ?
							(b(Ec[jj]*sqrt(pasoEaux),Rc[j],p,st,itRE)-
									b(Ec[jj]/sqrt(pasoEaux),Rc[j],p,st,itRE))/dEe : 0.0;
					double vR = abs(radialVel(Rc[jj]));
					double dRR = Rc[jj]*(pasoRc-1.0);
					mu += dbde * (dRR/vR);
				}
				mu = exp(-mu);
				double f1 = fAux(Rc[j]);
				double f2 = fAux(r);
				cout << f1 << "\t" << f2 << endl;
				double dRc = Rc[j]*(pasoRc-1.0);
				N += Qinj * P2(Rc[j]/r) * (f1/f2) * mu * dRc;
			}
			*/
			
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
					double f = fAux(Rc[jj],pasoRc);
					
					mu += f + dbde*(dRR/vR);
				}
				mu = exp(-mu);
				double dRc = Rc[j]*(pasoRc-1.0);
				N += Qinj * mu * dRc;
			}

			N /= v;

			/*double tCell = itR.val(DIM_R)*(sqrt(paso_r)-1.0/sqrt(paso_r)) / (-radialVel(r));
			cout << itR.coord[DIM_R] << "\t Ninj = "
									 << p.injection.get(itRE)*tCell << "\t N = " << N << endl;*/
			p.distribution.set(itRE,N);
		},{-1,itR.coord[DIM_R],0});
	},{0,-1,0});
}