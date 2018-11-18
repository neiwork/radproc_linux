#pragma once
#include <cstddef>

// PARAMETERS (from parameter.json)
double blackHoleMass;               	// Black hole mass [M_Sol]
double blackHoleSpinPar;            	// Specific black hole spin
double specificAngMomPar;				// Torus specific angular momentum parameter (lambda)
double iMeanMolecularWeight;			// Ion mean molecular weight
double eMeanMolecularWeight;			// Electron mean molecular weight
double magFieldPar;						// Magnetic field parameter (beta)
double eiTempRatio;						// Electron-Ion temperature ratio
double eCentralTemp;					// Electron temperature at the torus center [K]
double centralMassDens;					// Gas energy density at the torus center [g cm^-3]
double polytropIndex;					// Polytropic Index

double minPolarAngle;					// Minimum polar angle
double maxPolarAngle;					// Maximum polar angle
double logMinEnergy;					// Logarithm of the minimum energy for photons
double logMaxEnergy;					// Logarithm of the maximum energy for photons
std::size_t nR,nE,nTheta;					// Number of points in each dimension
double auxM0,auxM1;						// Auxiliary parameters

int numProcesses;

// PARAMETERS
double specificAngMom;					// Specific Angular Momentum
double gravRadius;						// Gravitationl Radius [cm]
double polytropConst;					// Polytropic Constant [P / e^{1+1/n}]
double torusCenterRadius;				// Radius of the torus center [rg]
double cuspRadius;						// Radius of the internal edge of the torus (cusp) [rg]
double edgeRadius;						// Radius of the external edge of the torus [rg]