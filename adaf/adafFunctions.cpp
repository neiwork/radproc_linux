#include "adafFunctions.h"
#include "globalVariables.h"
#include "read.h"
#include <iostream>
#include <fmath/physics.h>
#include <fparameters/parameters.h>

// Electrons

double electronTemp(double r)
{
	double logr_actual = log(r/schwRadius);
	double aux = logr.front();
	size_t pos_r = 0;
	while (aux < logr_actual) {
		pos_r++;
		aux = logr[pos_r];
	}
	double m = (logTe[pos_r]-logTe[pos_r-1])/(logr[pos_r]-logr[pos_r-1]);
	double logT = m*(logr_actual-logr[pos_r-1])+logTe[pos_r-1];
	double temp = eMeanMolecularWeight*exp(logT);
	double r_rg = r/(schwRadius/2.0);
	return (r_rg <= 30) ? temp * 1.4/pow(r_rg,0.097) : temp;
	//return temp;
}

// Ions
double ionTemp(double r) {
	double logr_actual = log(r/schwRadius);
	double aux = logr.front();
	size_t pos_r=0;
	while (aux < logr_actual) {
		pos_r++;
		aux = logr[pos_r];
	}
	double m = (logTi[pos_r]-logTi[pos_r-1])/(logr[pos_r]-logr[pos_r-1]);
	double logT = m*(logr_actual-logr[pos_r-1])+logTi[pos_r-1];
	
	return iMeanMolecularWeight*exp(logT);
}

double radialVel(double r) {
	double logr_actual = log(r/schwRadius);
	double aux = logr.front();
	size_t pos_r=0;
	while (aux < logr_actual) {
		pos_r++;
		aux = logr[pos_r];
	}
	double m = (logv[pos_r]-logv[pos_r-1])/(logr[pos_r]-logr[pos_r-1]);
	double logrv = m*(logr_actual-logr[pos_r-1])+logv[pos_r-1];
	double r_rg = r / (schwRadius/2.0);
	double vel = -exp(logrv);
	return (r_rg <= 30 ? vel / (0.93*exp(2.13/r_rg)) : vel);
}

double keplAngVel(double r)
{
	return sqrt(gravitationalConstant*blackHoleMass/r) / (r-schwRadius);
}

double sqrdSoundVel(double r)
{
    return boltzmann / (magFieldPar*atomicMassUnit) * ( ionTemp(r) / iMeanMolecularWeight 
							+ electronTemp(r) / eMeanMolecularWeight );
}

double height_fun(double r)
{
    return sqrt(sqrdSoundVel(r))/keplAngVel(r);
}

double angularVel(double r)
{
	return alpha*sqrdSoundVel(r)/(-radialVel(r)*r) + (j*schwRadius*cLight)/(r*r);
}

double costhetaH(double r)
{
	double cs = sqrt(sqrdSoundVel(r));
	double omR = angularVel(r)*r;
	double result = sqrt(pi/2.0) * cs/omR * erf(omR / (sqrt(2.0)*cs));
	return result;
}

double accRateADAF(double r)
{
	double accRateCorona = GlobalConfig.get<double>("accretionRateCorona")*accRateOut;
	double rOut = exp(logr.back())*schwRadius;
	//double result = accRateOut * (r/schwRadius > 5 ? pow(r/rOut,s) : pow(5*schwRadius/rOut,s));
	double result = accRateOut * pow(r/rOut,s);
	int processesFlags[numProcesses];
	readThermalProcesses(processesFlags);
	if (processesFlags[3]) {
		if (r > rTr) result = accRateCorona*pow(r/rOut,s)*(rTr/r);
	}
	return result;
}

double massDensityADAF(double r)
{
	//double result = accRateADAF(r) / (4.0*pi*r*r*costhetaH(r)*(-radialVel(r)));
	double result = accRateADAF(r) / (4.0*pi*r*height_fun(r)*(-radialVel(r)));
	//if (r < rBlob) result *= factorDensity;
	return result;
}

double magneticField(double r)
{
	return sqrt(8.0*pi*(1.0-magFieldPar)*massDensityADAF(r)*sqrdSoundVel(r));
}

double accRateColdDisk(double r)
{
	double rOut = exp(logr.back())*schwRadius;
	int processesFlags[numProcesses];
	readThermalProcesses(processesFlags);
	if (processesFlags[3] && rTr < 5.0) return accRateOut * (1.0-rTr/rOut);
	return accRateOut*(1.0-rTr/r);
}

double auxCD(double r)
{
	double lj = r/sqrt(paso_r);
	double lj1 = r*sqrt(paso_r);
	return 3.0*gravitationalConstant*blackHoleMass*accRateColdDisk(r)/2.0 *
			(1.0/lj * (1.0-2.0/3.0*sqrt(rTr/lj))-1.0/lj1*(1.0-2.0/3.0*sqrt(rTr/lj1)));
}

double electronDensity(double r)
{
	return (r > schwRadius) ? massDensityADAF(r)/(atomicMassUnit*eMeanMolecularWeight) : 0.0;
}

double electronDensityTheta(double r, double theta)
{
	double thetaMinLocal = acos(costhetaH(r));
	return (theta > thetaMinLocal && theta < pi-thetaMinLocal) ?
										electronDensity(r) : 0.0;
}

double ionDensity(double r)
{
	return massDensityADAF(r)/(atomicMassUnit*iMeanMolecularWeight);
}