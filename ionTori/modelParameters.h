#pragma once

#include <fmath/mathematics.h>
#include <fparticle/Particle.h>
#include <boost/property_tree/ptree_fwd.hpp>

//
//class Electron : public ParticleCfg<Electron> {};
//class Photon : public ParticleCfg<Photon> {};

const DimensionCoord
	DIM_E = 0,
	DIM_R = 1,
	DIM_THETA = 2;


/* define the inital values of the global parameters*/
void prepareGlobalCfg();

double eEmax(double z, double B);

//las fill no las estoy usando, lo calculo todo dentro de State
//void fillPhotonField(State& st);

//void fillMagnetic(State& st);

void initializePoints(Vector& v, double min, double max);

void torusParameters(double *l_0, double *rCusp, double *rCenter);

void initializeEnergyPoints(Vector& v, double logEmin, double logEmax);





