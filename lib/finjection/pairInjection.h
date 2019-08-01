#pragma once

#include <fparticle/Particle.h>


/* pairInjection calculates the pair injection due to photon-photon annihilation 
Aharonian, F. A.; Atoian, A. M.; Nagapetian, A. M. (1983); Vieyro & Romero 2012 */ 

double pairInjection(double E, const ParamSpaceValues& ntPh, const ParamSpaceValues& tpf, const SpaceCoord& distCoord, double tpEmin, double tpEmax);