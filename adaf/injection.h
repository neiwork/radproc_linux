#pragma once

//#include <fparticle\Particle.h>

#include "State.h"

/* Returns Q(E,z) in units of 1/erg/s  (it is multiplied by the volume of each celd)*/

void injection(Particle& p, State& st);
void injectionBurst(Particle& p, State& st);

double eEmax(Particle&,double,double,double);

void injectionChargedPion(Particle&, State&);
void injectionMuon(Particle&, State&);
void injectionPair(Particle& p, State& st);
