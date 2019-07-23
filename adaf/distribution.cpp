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

/*void distribution(Particle& p, State& st)
{
	show_message(msgStart, Module_electronDistribution);

	ParamSpaceValues N1(p.ps), N2(p.ps);

	N1.fill([&](const SpaceIterator& i){
		return 0.0;
	},{-1,-1,0});

	p.ps.iterate([&](const SpaceIterator& r_glob){
		
		const double E = r_glob.val(DIM_E);
		const double r = r_glob.val(DIM_R);
		double rB1 = r/sqrt(paso_r);
		double delta_r = rB1*(paso_r-1.0);
		double vel = -radialVel(r); 
		
		double tcell = delta_r/vel;
		double tcool = E/losses(E,p,st,r_glob);
		
		double dist1 = p.injection.get(r_glob)*min(tcell,tcool);
		N1.set(r_glob,dist1);   //hay injeccion solo en la primera celda de cada emisor lineal
		//N1.set(r_glob,safeLog10(dist1));
		
	},{-1,-1,0});

	for (int r_glob=0;r_glob<p.ps[1].size();r_glob++){

		N2.fill([&](const SpaceIterator& h){
			return 0.0;
		}, {-1,-1,0});  
		
		for (int r_i=r_glob-1;r_i>=0;r_i--) {
		
			double r = p.ps[DIM_R][r_i];
			double rB1 = r/sqrt(paso_r);
			double delta_r = rB1*(paso_r-1.0);
			double vel = -radialVel(r);
			double tcell = delta_r/vel;

			p.ps.iterate([&](const SpaceIterator& i){
				const double E = i.val(DIM_E);
				const double r = i.val(DIM_R);
			
				const double magf = st.magf.get(i);
				//double Emax = p.emax(); //VER eEmax(r, magf);
				double Emax = eEmax(p,r,magf,vel);
				
				double dist2 = 0.0;
				if (i.its[1].canPeek(+1)) { //if (r_ix != 0)
					double Eeff=Emax;

					if (i.its[0].canPeek(+1)) { 
						Eeff = effectiveE(E,Emax,tcell,p,st,i);
					}
				
					SpaceCoord i_plus_1 = i.moved({i.its[0].peek(0),+1,0}); //muevo el r uno para adentro
					//double r_plus_1 = i.its[1].peek(+1);
					double dist = N1.interpolate({{DIM_E,Eeff}},&i_plus_1);

					double ratioLosses = losses(Eeff,p,st,i_plus_1)/losses(E,p,st,i);
					dist2 = dist*ratioLosses;
				}
				N2.set(i, dist2);
			},{-1,r_i,0});
		}
		
		for (int E_i = 0; E_i < p.ps[0].size(); E_i++){
			
			double sum2 = 0.0;
		
			p.ps.iterate([&](const SpaceIterator& i2){
				sum2 += N2.get(i2);
				p.distribution.set(i2,sum2);
			}, { E_i, -1, 0 });
		}
	}

	for (int E_i = 0; E_i < p.ps[0].size(); E_i++){
		double sum1 = 0.0;
		p.ps.iterate([&](const SpaceIterator& i1){
			sum1 += N1.get(i1);
			double tot = sum1 + sum1;
			p.distribution.set(i1,tot);
		}, { E_i, -1, 0 });
	}
	
	show_message(msgEnd, Module_electronDistribution);
}*/

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
					SpaceCoord itRcoord = itRR;
					if (p.id == "ntElectron") {
						double integ = RungeKuttaSimple(energy,p.emax()*0.99,[&](double e) {
							return p.injection.interpolate({{0,e}},&itRcoord);});
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
	
	if (1) {
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
			SpaceCoord itAux = it;
			double int_Ec = pow(min(Ec[nPoints-1],p.emax())/max(min(Ec[0],p.emax()),p.emin()),1.0/nPoints);
			/*
			for (size_t j=0;j<nPoints;j++) {
				size_t k = nPoints-1-j;
				if ((Ec[k]/sqrt(int_Ec) > p.emin() && Ec[k]*sqrt(int_Ec) < p.emax()) &&
				(Rc[k] > p.ps[DIM_R][0] && Rc[k] < p.ps[DIM_R][nR-1])) {
					double Q = p.injection.interpolate({{DIM_E,Ec[k]},{DIM_R,Rc[k]}},&itAux);
					double dRc = Rc[k]*(sqrt(int_Rc)-1.0/sqrt(int_Rc));
					double dEc = Ec[k]*(sqrt(int_Ec)-1.0/sqrt(int_Ec));
					double dbde = (b(Ec[k]*sqrt(int_Ec),Rc[k],p,st,it) - 
									b(Ec[k]/sqrt(int_Ec),Rc[k],p,st,it))/dEc;
					double deriv = (Q-N*dbde) / radialVel(Rc[k]);
					N += deriv * (-dRc);
				}
			}
			p.distribution.set(it,N);
			*/

			for (size_t jE=0;jE<nPoints;jE++) {
				if ((Ec[jE] > p.emin() && Ec[jE] < p.emax()) &&
				(Rc[jE] > p.ps[DIM_R][0] && Rc[jE] < p.ps[DIM_R][nR-1])) {
					double mu = 0.0;
					for (size_t jjE=0;jjE<=jE;jjE++) {
						if ((Ec[jjE] > p.emin() && Ec[jjE] < p.emax()) &&
						(Rc[jjE]/sqrt(int_Rc) > p.ps[DIM_R][0] && 
						Rc[jjE]*sqrt(int_Rc) < p.ps[DIM_R][nR-1])) {
							double dEE = Ec[jjE]*(sqrt(int_Ec)-1.0/sqrt(int_Ec));
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
					double dEc = Ec[jE]*(sqrt(int_Ec)-1.0/sqrt(int_Ec));
					N += p.injection.interpolate({{DIM_E,Ec[jE]},{DIM_R,Rc[jE]}},&itAux)*mu*dEc;
				}
			}
			N = (e > p.emin() && e < p.emax()) &&
				(r > p.ps[DIM_R][0] && r) ? N/abs(b(e,r,p,st,it)) : 0.0;
			p.distribution.set(it,N);
		},{-1,-1,0});
	}

	if (p.id == "ntElectron") show_message(msgEnd,Module_electronDistribution);
	else if (p.id == "ntProton") show_message(msgEnd,Module_protonDistribution);
}





			/*for (int t_ix = 0; t_ix < p.ps[2].size(); t_ix++) {

		if (t_ix % 5 == 0) {
			std::cout << "time: " << t_ix << std::endl;
		}
		 
		//for (int z_ix = 0; z_ix < p.ps[1].size(); z_ix++) { //emisores para todo z
		//for (int z_ix = 0; z_ix <= t_ix; z_ix++) {	//unico emisor en z=0

		if (single) { superior = t_ix + 1; }
		//el +1 es para que el for incluya el t_ix

		for (int z_ix = 0; z_ix < superior; z_ix++) { */	


/*  paso 2
	
				double tp = t / Gamma; //time in the FF

				if (E < 10.0*Emax){
					double Eeff = effectiveE(E, Emax, tp, r, p, st, i);  //este tp se usa
					double dist1(0.0), dist2(0.0);

					dist1 = timeDistribution(E, r, tp, p, st, Eeff, i); //este tp NO se usa

					if (t_ix != 0)
					{	//estos son los puntos donde Q=0, y las particulas vienen de ti-1
						//if (i.its[2].canPeek(-1)) 

						double dist = N2.interpolate({ { DIM_E, Eeff }, { DIM_R, r }, { DIM_T, i.its[2].peek(-1) } });
						double ratioLosses = losses(Eeff, r, p, st, i) / losses(E, r, p, st, i);
						dist2 = dist*ratioLosses;
					}

					N12.set(i, dist1 + dist2); //lo cargo en N12 mientras interpolo de N2
					//N12.set(i, dist1); //lo cargo en N12 mientras interpolo de N2
				}
				else
				{
					N12.set(i, 0.0);  //si E>Emax entonces anulo la N
				}

			}, { -1, z_ix, t_ix });*/
			
			
	/*		paso 3
			
			
			p.ps.iterate([&](const SpaceIterator& i){

				double ni = N12.get(i); //el ni es que que obtengo con N12

				double dist3;

				double ri = i.its[1].peek(0);
				double rip1;

				if (i.its[1].canPeek(-1)) //if (z_ix != 0)
				{
					SpaceCoord coord = i.moved({ 0, -1, 0 }); //N(ri-1)
					double ni_1 = p.distribution.get(coord); //este lo calculo con el p.dist porque ya esta en r-1
				
					double ri_1 = i.its[1].peek(-1);
					
					if (i.its[1].canPeek(1)) { rip1 = i.its[1].peek(+1); }
					else{
						double r_int = pow((rmax / rmin), (1.0 / nR));
						rip1 = ri*r_int; 
					}

					dist3 = ni*(1.0 - ri / rip1) + ni_1*(ri_1 / ri);
				}
				else //z_ix=0
				{
					rip1 = i.its[1].peek(+1);
					dist3 = ni*(1.0 - ri / rip1);  
				}

				p.distribution.set(i,dist3); // lleno p.distribution e interpolo en N12

				//p.distribution.set(i, ni); // prueba
				
			} , { -1, z_ix, t_ix });*/