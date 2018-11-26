#pragma once
#include <cstddef>

// PARAMETERS (from parameter.json)
extern const double blackHoleMass;              // Black hole mass [M_Sol]
extern const double blackHoleSpinPar;           // Specific black hole spin
extern const double specificAngMomPar;			// Torus specific angular momentum parameter (lambda)
extern const double iMeanMolecularWeight;		// Ion mean molecular weight
extern const double eMeanMolecularWeight;		// Electron mean molecular weight
extern const double magFieldPar;				// Magnetic field parameter (beta)
extern const double eiTempRatio;				// Electron-Ion temperature ratio
extern const double eCentralTemp;				// Electron temperature at the torus center [K]
extern const double centralMassDens;			// Gas energy density at the torus center [g cm^-3]
extern const double polytropIndex;				// Polytropic Index

extern const double minPolarAngle;				// Minimum polar angle
extern const double maxPolarAngle;				// Maximum polar angle
extern const double logMinEnergy;				// Logarithm of the minimum energy for photons
extern const double logMaxEnergy;				// Logarithm of the maximum energy for photons
extern const size_t nR,nE,nTheta;			// Number of points in each dimension
extern const double auxM0,auxM1;				// Auxiliary parameters

// PARAMETERS
extern const double specificAngMom;				// Specific Angular Momentum
extern const double gravRadius;					// Gravitationl Radius [cm]
extern const double polytropConst;				// Polytropic Constant [P / e^{1+1/n}]
extern const double torusCenterRadius;			// Radius of the torus center [rg]
extern const double cuspRadius;					// Radius of the internal edge of the torus (cusp) [rg]
extern const double edgeRadius;					// Radius of the external edge of the torus [rg]

extern const int numProcesses;                  // Number of thermal processes

void torusParameters();