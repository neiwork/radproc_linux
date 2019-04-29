#pragma once
#include <cstddef>
#include <fmath/physics.h>

// PARAMETERS (from parameter.json)
double blackHoleMass;               	// Black hole mass [M_Sol]
double accRateOut;						// Accretion rate at rOut
double s;								// Wind index
double iMeanMolecularWeight;			// Ion mean molecular weight
double eMeanMolecularWeight;			// Electron mean molecular weight
double magFieldPar;						// Magnetic field parameter (beta)
double j;
double alpha;
double paso_r;
double paso_rCD;
double rTr;
double inclination;

double logMinEnergy;					// Logarithm of the minimum energy for photons
double logMaxEnergy;					// Logarithm of the maximum energy for photons
std::size_t nR,nE,nRcd;					// Number of points in each dimension

int numProcesses;
int calculateScatt;

size_t nRaux;
Vector logr;
Vector logTi;
Vector logTe;
Vector logv;

// PARAMETERS
double schwRadius;						// Gravitationl Radius [cm]

void adafParameters();