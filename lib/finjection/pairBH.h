#pragma once

#include <fparticle/Particle.h>



/* pairBH estimates the pair injection due to the channel 
for pair production of photohadronic interaction; Bethe Heitler*/ 
//double pairBH(double E, Particle& particle, Particle& proton, fun1 tpf);  

double pairBH(double E, const Particle& creator, const ParamSpaceValues& tpf, 
				const SpaceCoord& distCoord, double tpEmin, double tpEmax);