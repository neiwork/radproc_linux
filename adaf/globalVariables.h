#pragma once
#include <cstddef>
#include <fmath/physics.h>
#include <fparticle/Particle.h>
#include <boost/property_tree/ptree_fwd.hpp>

// PARAMETERS (from parameter.json)
extern const double blackHoleMass;              // Black hole mass
extern const double iMeanMolecularWeight;		// Ion mean molecular weight
extern const double eMeanMolecularWeight;		// Electron mean molecular weight
extern const double magFieldPar;				// Magnetic field parameter (beta)
extern const double accRateOut;					// Accretion rate at rOut
extern const double s;							// Wind index
extern const double alpha;
extern const double delta;
extern const double j;
extern const double paso_r;
extern const double paso_rCD;
extern const double rTr;
extern const double powerIndex;
extern const double rOutCD;
extern const double inclination;
extern double Ainjection;

extern const double logMinEnergy;				// Logarithm of the minimum energy for photons
extern const double logMaxEnergy;				// Logarithm of the maximum energy for photons
extern const size_t nR,nE,nRcd;					// Number of points in each dimension

extern const size_t nRaux;
extern const Vector logr;
extern const Vector logTi;
extern const Vector logTe;
extern const Vector logv;

extern Matrix scattAA;
extern Matrix scattDA;
extern Matrix reachAA;
extern Matrix reachAD;
extern Matrix reachDA;
extern Vector escapeAi;
extern Vector escapeDi;

extern Vector redshift_to_inf;
extern Vector redshift_CD_to_inf;
extern Matrix redshift;
extern Matrix redshift_RIAF_to_CD;
extern Matrix redshift_CD_to_RIAF;

extern const int calculateComptonScatt;
extern const int calculateThermal;
extern const int calculatePhotonDensityGap;
extern const int calculateNonThermal;
extern const int calculateLosses;
extern const int calculateNTdistributions;
extern const int calculateFlare;
extern const int calculateNonThermalLum;
extern const int calculateNonThermalHE;
extern const int calculateSecondaries;
extern const int calculateNTprotons;
extern const int calculateNTelectrons;
extern const int calculateNeutronInj;
extern const int calculateNeutronDis;
extern const int calculateJetDecay;
extern const int calculateJetEmission;
extern const int calculateComptonProbs;
extern const int calculateComptonRedMatrix;
extern const int comptonMethod;
extern const int calculateNewTemp;
extern const int height_method;
extern const int accMethod;

// PARAMETERS
extern const double schwRadius;					// Gravitationl Radius [cm]

// BLOB PARAMETERS
extern double rBlob;
extern double timeAfterFlare;
extern double factorDensity;
extern const double uth_power;

extern const int numProcesses;                  // Number of thermal processes

extern const size_t nGammaCompton;
extern const size_t nTempCompton;
extern const size_t nNuPrimCompton;
extern const size_t nNuCompton;
extern const double gammaMinCompton;
extern const double gammaMaxCompton;
extern const double tempMinCompton;
extern const double tempMaxCompton;
extern const double nuPrimMinCompton;
extern const double nuPrimMaxCompton;
extern const double nuMinCompton;
extern const double nuMaxCompton;

//FLARE

extern double maxRadius;
extern double minRadius;
extern double etaInj;
extern double pIndex;

extern const DimensionCoord
	DIM_E,
	DIM_R,
	DIM_Rcd;

void adafParameters();