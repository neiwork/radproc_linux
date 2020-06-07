#pragma once



#include "State.h"



/* returns dE/dt in units of erg/s*/

double losses(double E, Particle& p, State& st, const SpaceCoord& i);
double b(double E, double r, Particle& p, State& st, const SpaceCoord& i);
double t_cool(double E, double r, Particle& p);
double annihilationRate(double Ep, Particle& p, const SpaceCoord& psc);