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
	DIM_T = 2;

/* define the inital values of the global parameters*/
void prepareGlobalCfg();


double computeMagField(double z);

double eEmax(double z, double B);

double blackBody(double E, double z);


//las fill no las estoy usando, lo calculo todo dentro de State
//void fillPhotonField(State& st);

//void fillMagnetic(State& st);


void initializeRPoints(Vector& v, double Rmin, double Rmax);

void initializeCrossingTimePoints(Vector& time, double rMin, double rMax);

void initializeEnergyPoints(Vector& v, double logEmin, double logEmax);





