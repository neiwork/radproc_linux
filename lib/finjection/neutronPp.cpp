#include "neutronPp.h"

#include <fparameters/parameters.h>
#include <flosses/lossesHadronics.h>
//#include <fmath/interpolation.h>


//double neutronPp(double E, Vector Nproton, Particle& particle, Particle& proton)  

double neutronPp(double E, Particle& p ,
	const double density, const SpaceCoord& psc)

{
	
	double protonDist = p.distribution.interpolate({ { 0, E } }, &psc);  // proton.dist(E);// interpol(E, proton.energyPoints, Nproton, Nproton.size() - 1);
		
	double t_1   = lossesHadronics(E, density, p)/E;
	
	double emissivity = 0.5*t_1*protonDist;
	
	return emissivity;
}

