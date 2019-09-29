#pragma once

#include "State.h"

void opticalDepth(Vector&,int,State&,int);
double ssaAbsorptionCoeff(double,double,Particle&,SpaceCoord&);

//double internalAbs(int E_ix, State& st, double r_current) ;