#include "distribution.h"

#include "messages.h"

#include "adafFunctions.h"
#include "globalVariables.h"
#include "losses.h"
//#include "timeDistribution.h"

#include <iostream>

#include <fparameters/Dimension.h>
#include <fparameters/SpaceIterator.h>
//#include <fparameters/parameters.h>
//#include <fmath/physics.h>
//#include <boost/property_tree/ptree.hpp>


double effectiveE(double Ee, double Emax, double t, Particle& p, State& st, const SpaceCoord& i)
{
	double Eeffmin, Eeffmax, Eeff_int, sum_Eeff, dEeff, dtau_eff;

	Eeffmin = Ee;
	Eeffmax = Emax;
	double Eeff = Eeffmin;

	double nEeff = 50;

	Eeff_int = pow((Eeffmax / Eeffmin), 1.0/nEeff);  //0.00001
	sum_Eeff = 0.0;

	while ((sum_Eeff < t) && ( (Eeffmax-Eeff) < 1.0e-2)){//  (Eeff < Eeffmax))	{// for (int k=0; k<=100000; ++k)	{

		dEeff = Eeff*(Eeff_int - 1.0);
		dtau_eff = dEeff / losses(Eeff, p, st, i); //losses(Eeff, r, p, st, i);
		sum_Eeff = sum_Eeff + dtau_eff;

		//	if (sum_Eeff > times[j]) goto

		Eeff = Eeff*Eeff_int;

	}

	return Eeff;

}


void distribution(Particle& p, State& st)
{

	show_message(msgStart, Module_electronDistribution);


	ParamSpaceValues N1(p.ps);
	N1.fill([&](const SpaceIterator& i){
		return 0.0;
	}, { -1, -1, 0 });

	ParamSpaceValues N2(p.ps);

	p.ps.iterate([&](const SpaceIterator& r_glob){
		
		const double E = r_glob.val(DIM_E);
		const double r = r_glob.val(DIM_R);
		
		double rB1   = r/sqrt(paso_r);
		double rB2   = r*sqrt(paso_r);
		double delta_r = rB2-rB1;
		double vel	 = -radialVel(r); 
		
		double tcell = delta_r/vel;
		double tcool = E/losses(E, p, st, r_glob);
		
		double dist1 = p.injection.get(r_glob)*min(tcell,tcool); //Q(E)×min(tcell,tcool), where tcell=scell(θ)/v‖(θ)

		N1.set(r_glob,dist1);   //hay injeccion solo en la primera celda de cada emisor lineal 
		
	}, { -1, -1, 0 });

	

	for (int r_glob = 0; r_glob < p.ps[1].size(); r_glob++){
						
		N2.fill([&](const SpaceIterator& h){
			return 0.0;
		}, { -1, -1, 0 });  
		
		
		for (int r_i = r_glob-1; r_i >= 0; r_i--){  //(int z_ix = 0; z_ix < r_i; z_ix++){ //p.ps[1].size(); z_ix++) { //inicio de los emisores lineales
		
		
			double r = p.ps[DIM_R][r_i];
			double rB1   = r/sqrt(paso_r);
			double rB2   = r*sqrt(paso_r);
			double delta_r = rB2-rB1;
			double vel	 = -radialVel(r); 
			
			double tcell = delta_r/vel;
		
		
			p.ps.iterate([&](const SpaceIterator& i){
				const double E = i.val(DIM_E);
				const double r = i.val(DIM_R);
			
				const double magf = st.magf.get(i);
				
				double Emax = p.emax(); //VER eEmax(r, magf);
				
				double dist2 = 0.0;
				
				if (i.its[1].canPeek(+1))  //if (r_ix != 0)
					{	
						
						double Eeff;

						if (i.its[0].canPeek(+1)){ Eeff = effectiveE(E, Emax, tcell, p, st, i);}
						else { Eeff = Emax; }
						
						SpaceCoord i_plus_1 = i.moved({ i.its[0].peek(0), +1, 0 }); //muevo el r uno para adentro
						
						double r_plus_1 = i.its[1].peek(+1);
						double dist =0.0; // N1.interpolate({ { DIM_E, Eeff }}, &i_plus_1); //{ DIM_R, r_plus_1 }}, &i);
						
						double ratioLosses = losses(Eeff, p, st, i_plus_1) / losses(E, p, st, i);
						
						dist2 = dist*ratioLosses;
					}

					N2.set(i, dist2);
				
				}, { -1, r_i, 0 });

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