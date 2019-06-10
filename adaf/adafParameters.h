#pragma once
#include <cstddef>
#include <fmath/physics.h>

// PARAMETERS (from parameter.json)
double blackHoleMass;               	// Black hole mass [g].
double accRateOut;						// Accretion rate at rOut.
double s;								// ADAF wind index.
double iMeanMolecularWeight;			// Ion mean molecular weight.
double eMeanMolecularWeight;			// Electron mean molecular weight.
double magFieldPar;						// Magnetic field parameter (beta).
double j;								// ADAF angular momentum at rS.
double alpha;							// Viscosity parameter.
double paso_r;							// logarithmic step in rADAF.
double paso_rCD;						// logarithmic step in rCD.
double rTr;								// Transition radius.
double inclination;						// Inclination of the line of sight.

double logMinEnergy;					// Logarithm of the minimum energy for photons.
double logMaxEnergy;					// Logarithm of the maximum energy for photons.
std::size_t nR,nE,nRcd;					// Number of points in each dimension.

int numProcesses;						// Number of thermal radiative processes.
int calculateScatt;						// If 1 calculate the Compton scattering 
										// probability matrices.
int calculateComptonRedMatrix;			// If 1 calculate the Comptonization energy redistribution
										// matrices. It only works if comptonMethod=0 or 1.
int comptonMethod;						// Method used to calculate the Comptonized spectrum.
										// 0: Interpolation using precomputed P(v,v',Te).
										// 1: Interpolation using precomputed P(v,v',Y).
										// 2: Local Comptonization.
										// 3: Compute the P(v,v',Te) at the energies required.

size_t nGammaCompton;					// Number of points in Y  for the precomputed matrices
										// (only if comptonMethod = 1).
size_t nTempCompton;					// Number of points in Temps for the precomputed matrices
										// (only if comptonMethod = 0).
size_t nNuPrimCompton;					// Number of points in v' for the precomputed matrices.
size_t nNuCompton;						// Number of points in v for the precomputed matrices.
double gammaMinCompton;					// Minimum Y for the precomputed matrices (only if 
										// comptonMethod = 1).
double gammaMaxCompton;					// Maximum Y for the precomputed matrices (only if
										// comptonMethod = 1).
double tempMinCompton;					// Minimum Te (norm.) for the precomputed matrices (only 
										// if comptonMethod = 0).
double tempMaxCompton;					// Maximum Te (norm.) for the precomputed matrices (only 
										// if comptonMethod = 0).
double nuPrimMinCompton;				// Minimum v' for the precomputed matrices.
double nuPrimMaxCompton;				// Maximum v' for the precomputed matrices.
double nuMinCompton;					// Minimum v  for the precomputed matrices.
double nuMaxCompton;					// Maximum v for the precomputed matrices.

size_t nRaux;							// Number of radius points in the precomputed ADAF.
Vector logr;							// Precomputed ADAF radii.
Vector logTi;							// Precomputed ADAF Ti.
Vector logTe;							// Precomputed ADAF Te.
Vector logv;							// Precomputed ADAF	radial velocity.

// PARAMETERS
double schwRadius;						// Schwarzschild Radius [cm]

void adafParameters();