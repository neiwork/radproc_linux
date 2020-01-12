#include "lossesHadronics.h"


#include "crossSectionInel.h"
#include <fparameters/parameters.h>
#include <fmath/physics.h>

double lossesHadronics(double E, double density, Particle& particle)
{
	double inelasticity = 0.5; 
	if(particle.id == "ntPion")
		return (2.0/3.0)*cLight*density*inelasticity*crossSectionHadronic(E)*E;
	else	//particle == proton or neutron
		return cLight*density*inelasticity*crossSectionHadronic(E)*E;
}