#include "injection.h"

#include "messages.h"
#include "globalVariables.h"
//#include "modelParameters.h"

#include <fparameters/parameters.h>
#include <fparameters/SpaceIterator.h>
#include <fparameters/Dimension.h>
#include <fmath/fbisection.h>

#include <fmath/RungeKutta.h>
#include <fmath/physics.h>

#include <boost/property_tree/ptree.hpp>

#include <iostream>

double powerLaw(double E, double Emin, double Emax)
{
	static const double primaryIndex = GlobalConfig.get<double>("primaryIndex");
	
	double result = pow(E, (-primaryIndex))*exp(-E / Emax)*exp(-5 * Emin / E);
	return result;
}



void injection(Particle& p, State& st)
{
	
	static const double etaInj = GlobalConfig.get<double>("etaInj");



	double Emin1 = p.emin();  //esta es la primera que uso de prueba
	double Emax = p.emax();

	double int_E = RungeKuttaSimple(Emin1, Emax, [&Emax, &Emin1](double E) {
		return E*powerLaw(E, Emin1, Emax);
	});  //integra E*Q(E)  entre Emin y Emax


	//p.injection.fill([&](const SpaceIterator& i) {
	p.ps.iterate([&](const SpaceIterator& i) {

		const double r = i.val(DIM_R);
		
		double rB1   = r/sqrt(paso_r);
		double rB2   = r*sqrt(paso_r);
		double vol_i = (rB2*rB2*rB2-rB1*rB1*rB1)*(4.0/3.0)*pi; //*cos(thetaH);
		double area  = 4.0*pi*rB2*rB2; //*(2.0*cos(thetaH)+sin(thetaH)*sin(thetaH));
		
		double uth,Emin;
		
		if(p.id == "ntElectron"){
			double norm_temp = boltzmann*st.tempElectrons.get(i)/(p.mass*cLight2);
			double dens = st.denf_e.get(i); 
			
			uth = dens*boltzmann*norm_temp*(p.mass*cLight2)*3.0/2.0; //cambiar 3/2 por factor a(theta_e)   erg cm^-3
		
			Emin =  fbisection ([&](double E)    //({ { 0, eval } }, &psc)
			{return st.electron.distribution.interpolate({ {DIM_E, E },{ DIM_R, r },{ 3, 0 } }) - powerLaw(E, Emin, Emax)*uth/int_E; },
			//aca normalizo al Q' con unidades de erg^-1 cm^-3, para compararlo con la termica
									st.electron.emin(), st.electron.emax(), 1.0e-2);
								
			
		}
		else if(p.id == "ntProton"){
			double norm_temp = boltzmann*st.tempIons.get(i)/(p.mass*cLight2);
			double dens = st.denf_i.get(i); 
			
			uth = dens*boltzmann*norm_temp*(p.mass*cLight2)*3.0/2.0; //cambiar 3/2 por factor a(theta_e)   erg cm^-3
		
			//fbisection(fun1 func,double x1,double x2,double xacc)
			Emin =  fbisection ([&](double E) 
			{return st.proton.distribution.interpolate({ {DIM_E, E },{ DIM_R, r },{ 3, 0 } }) - powerLaw(E, Emin, Emax)*uth/int_E; },
			//aca normalizo al Q' con unidades de erg^-1 cm^-3, para compararlo con la termica
									st.proton.emin(), st.proton.emax(), 1.0e-2);
			
		}
		
		double Q0 = etaInj*uth*vol_i/(area*cLight);   //power injected in nt particles
	
		double Q0p = Q0/ (int_E);


		p.ps.iterate([&](const SpaceIterator& jE) {
				
			const double E = jE.val(DIM_E);
			double total = powerLaw(E, Emin, Emax)*Q0p; 

			p.injection.set(jE,total); //en unidades de erg^-1 s^-1
		
		}, { -1, i.coord[DIM_R],0} );
		
	}, { 0, -1, 0} );

}
