#pragma once

#include <fmath/mathematics.h>
#include <fmath/physics.h>
#include <fparticle/Particle.h>
#include <boost/property_tree/ptree_fwd.hpp>

const DimensionCoord
	DIM_E = 0,
	DIM_R = 1,
	DIM_Rcd = 2;
	
/* define the inital values of the global parameters*/
void prepareGlobalCfg();
void initGridLogarithmically(Vector& v, double min, double max);
void initGridLinearly(Vector& v, double min, double max);
void initEnergyPoints(Vector& v, double logEmin, double logEmax);





