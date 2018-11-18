#pragma once

#include <fmath/mathematics.h>
#include <fparticle/Particle.h>
#include <boost/property_tree/ptree_fwd.hpp>

//class Electron : public ParticleCfg<Electron> {};
//class Photon : public ParticleCfg<Photon> {};
const DimensionCoord
	DIM_E = 0,
	DIM_R = 1,
	DIM_THETA = 2;
/* define the inital values of the global parameters*/
void prepareGlobalCfg();
void initThetaPoints(Vector& v, double min, double max);
void initRadiiPoints(Vector& v, double min, double max);
void initEnergyPoints(Vector& v, double logEmin, double logEmax);
void initGridLinearly(Vector& v, double min, double max);
void initGridLogarithmically(Vector& v, double min, double max);





