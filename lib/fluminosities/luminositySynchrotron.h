#pragma once

#include <fparticle/Particle.h>


/* this is the Synchrotron luminosity without Synchrotron self absorption; erg/s/cm^3*/ 
double luminositySynchrotron(double E, const Particle& c, const SpaceCoord& distCoord, const ParamSpaceValues& magf);
double luminositySynchrotron2(double E, const Particle& c, const SpaceCoord& distCoord, double magf);
double luminositySynchrotronExact(double E, const Particle& c, const SpaceCoord& distCoord, double magf);
double luminositySynchrotronExactSec(double E, const Particle& c, const SpaceCoord& distCoord, double magf);

/* this is the Synchrotron luminosity with Synchrotron self absorption*/ 
//double luminositySynchrotron_conSSA(double E, const Particle& creator);


/* this is the Synchrotron luminosity produced by secondary pairs */ 
double luminositySynchrotronSec(double E, const Particle& creator, const SpaceCoord& psc, double magf);