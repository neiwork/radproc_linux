#pragma once

#include <fparticle/Particle.h>


/*erg/s/cm^3 */

//double luminosityIC(double E, const Particle& creator, const SpaceCoord& distCoord, fun1 tpf, double phEmin);
double luminosityIC(double E, const Particle& creator, const SpaceCoord& distCoord, const ParamSpaceValues& tpf,
					double tpEmin, double tpEmax);