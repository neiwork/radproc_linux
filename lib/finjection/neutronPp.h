#pragma once

#include <fparticle/Particle.h>

//double neutronPp(double E, Vector Nproton, Particle& particle, Particle& proton);

double neutronPp(double E, Particle& p,
	const double density, const SpaceCoord& psc);