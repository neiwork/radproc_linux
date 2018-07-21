#pragma once
#include <fmath/mathematics.h>
#include <boost/property_tree/ptree_fwd.hpp>

// Global variables
double massBH;                // Black hole mass
double RS;                    // Schwarzschild Radius
double alphapar;              // Viscosity parameter
double betapar;               // Magnetic parameter
double gammapar;              // Polytropic index
double rmin;                  // Inner edge of the ADAF
double rmax;                  // Outer edge of the ADAF
double f;                     // Advection parameter
double r;                     // Radius

const DimensionCoord
	DIM_E = 0,
	DIM_R = 1,
	DIM_THETA = 2;
	
/* define the inital values of the global parameters*/
void prepareGlobalCfg();
void initializeThetaPoints(Vector& v, double min, double max);
void initializeRadiiPoints(Vector& v, double min, double max);
void initializeEnergyPoints(Vector& v, double logEmin, double logEmax);





