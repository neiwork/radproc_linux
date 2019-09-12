#include "neutronDecay.h"


#include <fmath/RungeKutta.h>
#include <fmath/physics.h>


double f_NeuDec(double Epi, double E, double mass, const Particle& neutron, const SpaceCoord& distCoord)         //funcion a integrar variable Epi
{
	
	double distNeutron;
	if (Epi < neutron.emin() || Epi > neutron.emax()){
		distNeutron = 0.0;
	}
	else{
		distNeutron = neutron.distribution.interpolate({ { 0, Epi } }, &distCoord); 
	}
	//double distNeutron = interpol(Epi,Ecreator,Ncreator,Ncreator.size()-1);

	double r = P2(mass/neutronMass);

	double x = E/Epi;
	double decayTime = neutronMeanLife*Epi/(neutronMass*cLight2);

	double Q = distNeutron/(Epi*decayTime); //*r*(1-x)/(Epi*x*P2(1-r)*decayTime);   

	return Q;   
}


double injNeutronDecay(double E, Particle& p, Particle& neutron, const SpaceCoord& distCoord)  
{

	double Emax = neutron.emax();  
	
	double mass = p.mass;
	
	//.Ecreator = neutron.eDim()->values; //el creator es neutron

	double injection = 0.0;

	double rpi = P2(mass/neutronMass);

	double sup = std::min(Emax,E/rpi);  // transformo la condicion de la heaviside en un limite superior

	injection = RungeKuttaSimple(E, sup, 
				[E,mass,&neutron,&distCoord](double Epi){
						return f_NeuDec(Epi,E,mass,neutron,distCoord);
					});
		
	return injection;

}




/* la de arriba era protonNeutron, lo de abajo es de electrones, son lo mismo, cambia la masa

double fQelectron(double Epi, double E, const Particle& p, const SpaceCoord& distCoord)         //funcion a integrar variable Epi
{
	double mass = p.mass;
	
	double distNeutron = interpol(Epi,Ecreator,Ncreator,Ncreator.size()-1);

	double r = P2(mass/neutronMass);

	double x = E/Epi;
	double decayTime = neutronMeanLife*Epi/(neutronMass*cLight2);

	double Q = distNeutron/(Epi*decayTime); //*r*(1-x)/(Epi*x*P2(1-r)*decayTime);   

	return Q;   
}


double electronNeutron(double E, Vector Ncreator, Particle& particle, Particle& neutron)  
{
	ParticleType particleName = particle.type; 

	double Emax = 1.6e-12*pow(10.0,neutron.logEmax);  

	double mass = p.mass;

	DataInjection data;

	data.E        = E;
	data.mass     = particle.mass;
	data.Ncreator = Ncreator;
	data.Ecreator = neutron.energyPoints;

	double injection = 0.0;

	double rpi = P2(electronMass/neutronMass);

	double sup = std::min(Emax,E/rpi);  // transformo la condicion de la heaviside en un limite superior

	injection = RungeKuttaSimple(E, sup, fQelectron, &data);
		
	return injection;

}*/