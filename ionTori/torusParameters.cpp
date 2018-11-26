#include "torusParameters.h"
#include "torusFunctions.h"
#include <math.h>
#include <iostream>
#include <fmath/physics.h>
#include <fparameters/parameters.h>
#include <boost/property_tree/ptree.hpp>

using namespace std;

///////////////////////

void specificAngularMomentum(double r_ms, double r_mb)
{
	double l_ms = keplAngularMom(r_ms);          // Keplerian specific angular momentum at r = r_ms
	double l_mb = keplAngularMom(r_mb);          // Keplerian specific angular momentum at r = r_mb
	specificAngMom = (1.0 - specificAngMomPar) * l_ms + specificAngMomPar * l_mb;
}

void criticalRadii(double r_ms, double r_mb)
{
	int maxmitr = 1000;
	double allerr = 1.0e-3;
	
	cuspRadius = bisection(r_mb, r_ms, allerr, maxmitr,modfKepl);
	torusCenterRadius = bisection(r_ms, 20.0*r_ms, allerr, maxmitr, modfKepl);
	edgeRadius = bisection(torusCenterRadius,10.0*torusCenterRadius,allerr,maxmitr,modfw);
}

void torusParameters() 
{	
	blackHoleMass = GlobalConfig.get<double>("blackHoleMass")*solarMass;
	blackHoleSpinPar = GlobalConfig.get<double>("blackHoleSpinPar");
	eMeanMolecularWeight = GlobalConfig.get<double>("eMeanMolecularWeight");
	iMeanMolecularWeight = GlobalConfig.get<double>("iMeanMolecularWeight");
	eiTempRatio = GlobalConfig.get<double>("eiTempRatio");
	eCentralTemp = GlobalConfig.get<double>("eCentralTemp");
	polytropIndex = GlobalConfig.get<double>("polytropIndex");
	magFieldPar = GlobalConfig.get<double>("magFieldPar");
	centralMassDens = GlobalConfig.get<double>("centralMassDens");
	specificAngMomPar = GlobalConfig.get<double>("specificAngMomPar");
	nR = GlobalConfig.get<int>("model.particle.default.dim.radius.samples");
	nE = GlobalConfig.get<int>("model.particle.default.dim.energy.samples");
	nTheta = GlobalConfig.get<int>("model.particle.default.dim.theta.samples");
	logMinEnergy = GlobalConfig.get<double>("model.particle.default.dim.energy.min");
	logMaxEnergy = GlobalConfig.get<double>("model.particle.default.dim.energy.max");
	minPolarAngle = GlobalConfig.get<double>("model.particle.default.dim.theta.min");
	maxPolarAngle = GlobalConfig.get<double>("model.particle.default.dim.theta.max");
	numProcesses = GlobalConfig.get<int>("numProcesses");
	
	// DERIVED CONSTANTS
	auxM0 = iMeanMolecularWeight / (eMeanMolecularWeight + iMeanMolecularWeight);
	auxM1 = iMeanMolecularWeight*eiTempRatio /
			(eMeanMolecularWeight+iMeanMolecularWeight*eiTempRatio);
	polytropConst = boltzmann*eCentralTemp/( (1.0-magFieldPar)*atomicMassUnit
			*pow(centralMassDens,2.0/3.0)*eMeanMolecularWeight*auxM1 );
	gravRadius = gravitationalConstant*blackHoleMass / cLight2;
	
	double r_ms,r_mb;
	marginalOrbits(r_ms, r_mb);
	specificAngularMomentum(r_ms,r_mb);
	criticalRadii(r_ms,r_mb);
}