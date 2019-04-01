#include "adafParameters.h"
#include "adafFunctions.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <fmath/physics.h>
#include <fparameters/parameters.h>
#include <boost/property_tree/ptree.hpp>

using namespace std;

///////////////////////

void adafParameters() 
{	
	ifstream adafFile, adafParams;
	
	adafFile.open("adafFile.txt"); adafParams.open("adafParameters.txt");
	adafFile >> nRaux;
	double accRateNorm;
	adafParams >> blackHoleMass >> accRateNorm >> s >> magFieldPar >> alpha >> j;
	adafParams.close();
	
	logr.resize(nRaux,0.0);
	logTi.resize(nRaux,0.0);
	logTe.resize(nRaux,0.0);
	logv.resize(nRaux,0.0);
	
	for (size_t i=0;i<nRaux;i++) {
		adafFile >> logr[i] >> logTi[i] >> logTe[i] >> logv[i];
	}
	adafFile.close();
	
	blackHoleMass *= solarMass;
	double eddAccRate = 1.39e18 * blackHoleMass/solarMass;
	accRateOut = accRateNorm*eddAccRate;
	eMeanMolecularWeight = GlobalConfig.get<double>("mu_e");
	iMeanMolecularWeight = GlobalConfig.get<double>("mu_i");
	nR = GlobalConfig.get<int>("model.particle.default.dim.radius.samples");
	nE = GlobalConfig.get<int>("model.particle.default.dim.energy.samples");
	logMinEnergy = GlobalConfig.get<double>("model.particle.photon.dim.energy.min");
	logMaxEnergy = GlobalConfig.get<double>("model.particle.photon.dim.energy.max");
	numProcesses = GlobalConfig.get<int>("numProcesses");
	
	paso_r = pow(exp(logr.back())/1.1,1.0/nR);
	// DERIVED CONSTANTS
	schwRadius = 2.0*gravitationalConstant*blackHoleMass / cLight2;
}