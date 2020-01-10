#pragma once


#include <fparticle/Particle.h>




//double lossesIC(double E, Particle& particle, fun1 tpf, double phEmin, double phEmax);

double lossesIC(double E, Particle& particle, const ParamSpaceValues& tpf, 
				const SpaceCoord& psc, double phEmin, double phEmax);
double lossesIC_Th(double E, Particle& particle, const ParamSpaceValues& tpf, 
				const SpaceCoord& psc, double phEmin, double phEmax);
double lossesICnew(double E, Particle& particle, const ParamSpaceValues& tpf, 
				const SpaceCoord& psc, double phEmin, double phEmax);