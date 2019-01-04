#include "torusParameters.h"
#include "torusFunctions.h"
#include <math.h>
#include <iostream>
#include <fmath/physics.h>
#include <fparameters/parameters.h>
#include <boost/property_tree/ptree.hpp>

#include <fmath/brent.h>

using namespace std;

///////////////////////

void specificAngularMomentum(double r_ms, double r_mb)
{
	double l_ms = keplAngularMom(r_ms);          // Keplerian specific angular momentum at r = r_ms
	double l_mb = keplAngularMom(r_mb);          // Keplerian specific angular momentum at r = r_mb
	specificAngMom = specificAngMomPar * (l_mb-l_ms) + l_ms;
}

void criticalRadii(double r_ms, double r_mb)
{
	int status1,status2,count=0;
	double x_lo,x_hi;
	
	struct zero_params modfKepl_params = {};
	struct zero_params modfNormPot_params = {};
	gsl_function gsl_modfKepl, gsl_modfNormPot;
		gsl_modfKepl.function= &modfKepl;
		gsl_modfKepl.params = &modfKepl_params;
		gsl_modfNormPot.function = &modfNormPot;
		gsl_modfNormPot.params = &modfNormPot_params;
		
	x_lo = r_mb;	x_hi = r_ms;
	cuspRadius = brent(&gsl_modfKepl,x_lo,x_hi,&status1,&status2);
	
	x_lo = r_ms;	x_hi = x_lo;
	do {
		x_hi += 1.0;
	} while ((keplAngularMom(x_hi)-specificAngMom)*(keplAngularMom(x_lo)-specificAngMom) > 0.0);
	torusCenterRadius = brent(&gsl_modfKepl,x_lo,x_hi,&status1,&status2);
	
	x_lo = torusCenterRadius;	x_lo = x_hi;
	do {
		x_hi += 10.0;
	} while (normalizedPotential(x_lo,0.0)*normalizedPotential(x_hi,0.0) > 0.0 && ++count < 1000);
	edgeRadius = brent(&gsl_modfNormPot,x_lo,x_hi,&status1,&status2);
	double aux=edgeRadius;
	
	double s=0.;
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