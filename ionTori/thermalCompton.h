#pragma once

#include <fparticle/Particle.h>
//#include <fparameters/SpaceIterator.h>
//#include <fparameters/Dimension.h>
//#include <fparameters/parameters.h>

double thCompton(double E, const Particle& creator, const SpaceCoord& distCoord, const ParamSpaceValues& tpf, double phEmin);