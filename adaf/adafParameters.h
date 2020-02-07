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
double delta;							// Fraction of energy released by viscosity that directly
										// heats electrons.
double paso_r;							// logarithmic step in rADAF.
double paso_rCD;						// logarithmic step in rCD.
double rTr;								// Transition radius.
double rOutCD;							// Outer radius of the cold disk.
double inclination;						// Inclination of the line of sight.
double Ainjection;						// Ainj;

double logMinEnergy;					// Logarithm of the minimum energy for photons.
double logMaxEnergy;					// Logarithm of the maximum energy for photons.
std::size_t nR,nE,nRcd;					// Number of points in each dimension.

int numProcesses;						// Number of thermal radiative processes.
int calculateComptonScatt;				// If 1 calculate the Compton scattering 
										// probability matrices.
int calculateThermal;   				// If 1 calculate the thermal processes.
int calculatePhotonDensityGap;          // If 1 calculate the photon density at zGap.
int calculateNonThermal;   				// If 1 calculate the non-thermal processes.
int calculateNTprotons;					// If 1 calculate the non-thermal processes for protons.
int calculateNTelectrons;				// If 1 calculate the non-thermal processes for electrons.
int calculateLosses;                    // If 1 calculate radiative losses for electrons and protons.
int calculateNTdistributions;           // If 1 calculate non-thermal particle distributions.
int calculateNonThermalLum;             // If 1 calculate non-thermal luminosities.
int calculateSecondaries;				// If 1 calculate electromagnetic cascades.
int calculateFlare;						// If 1 calculate non-thermal flare distributions.
int calculateNeutronInj;                // If 1 calculate neutron injection.
int calculateNeutronDis;				// If 1 calculate neutron propagation.
int calculateJetDecay;					// If 1 calculate power injection in jet by neutron decay.
int calculateJetEmission;				// If 1 calculate emission of a jet.
int calculateComptonRedMatrix;			// If 1 calculate the Comptonization energy redistribution
										// matrices. It only works if comptonMethod=0 or 1.
int comptonMethod;						// Method used to calculate the Comptonized spectrum.
										// 0: Interpolation using precomputed P(v,v',Te).
										// 1: Interpolation using precomputed P(v,v',Y).
										// 2: Local Comptonization.
										// 3: Compute the P(v,v',Te) at the energies required.
int height_method;						// 0: Use r*costhetaH(r)
										// 1: Use sqrt(cs/omK)
int calculateNewTemp;					// If 1 calculate again the electron temperatures.

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

Matrix scattAA;							// Probab. for a photon emitted in an ADAF shell
										// to be Compton-scattered in another shell.
Matrix scattDA;							// Probab. for a photon emitted in a thin disk shell to
										// be Compton-scattered in an ADAF shell.
Matrix reachAD;							// Probab. for a photon emitted in an ADAF shell to
										// reach a thin disk shell.
Matrix reachAA;							// Probab. for a photon emitted in an ADAF shell
										// to reach another shell.
Matrix reachDA;							// Probab. for a photon emitted in a thin disk shell to
										// reach an ADAF shell.
Vector escapeAi;						// Probab. for a photon emitted in an ADAF shell to escape
										// to infinity at an angle i.
Vector escapeDi;						// Probab. for a photon emitted in a thin disk shell to
										// escape to infinity at an angle i.

// PARAMETERS
double schwRadius;						// Schwarzschild Radius [cm]

// BLOB PARAMETERS
double rBlob;
double timeAfterFlare;
double factorDensity;

double maxRadius;
double minRadius;
double etaInj;
double pIndex;

void adafParameters();