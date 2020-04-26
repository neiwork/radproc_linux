#pragma once

//#include <fparticle\Particle.h>

#include "State.h"

/* Returns Q(E,z) in units of 1/erg/s  (it is multiplied by the volume of each celd)*/

double findGammaMin(double,double);
void injection(Particle& p, State& st);
void injectionBurst(Particle& p, State& st);

double eEmax(Particle&,double,double,double,double,int);
double eEmax_numerical(Particle& p, double r, double B, double v, double dens, int jR, State& st,
						SpaceCoord i);

void injectionChargedPion(Particle&, State&);
void injectionMuon(Particle&, State&);
void injectionNeutrino(Particle&, State&);
void injectionPair(Particle& p, State& st, int);
