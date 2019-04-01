#include "adafFunctions.h"
#include "globalVariables.h"
#include <fmath/physics.h>

// Electrons

double electronTemp(double r)
{
	double logr_actual = log(r/schwRadius);
	double aux = logr.front();
	size_t pos_r=0;
	while (aux < logr_actual) {
		pos_r++;
		aux = logr[pos_r];
	}
	double m = (logTe[pos_r]-logTe[pos_r-1])/(logr[pos_r]-logr[pos_r]);
	double logT = m*(logr_actual-logr[pos_r-1])+logTe[pos_r-1];
	return eMeanMolecularWeight*exp(logT);
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
	double m = (logTi[pos_r]-logTi[pos_r-1])/(logr[pos_r]-logr[pos_r]);
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
	double m = (logv[pos_r]-logv[pos_r-1])/(logr[pos_r]-logr[pos_r]);
	double logrv = m*(logr_actual-logr[pos_r-1])+logv[pos_r-1];
	return -exp(logrv);
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
	return sqrt(pi/2.0) * cs/omR * erf(omR / (sqrt(2.0)*cs));
}

double massDensityADAF(double r)
{
	double rOut = exp(logr.back())*schwRadius;
	return accRateOut*pow(r/rOut,s) / 
		(4.0*pi*r*r*(-radialVel(r))*costhetaH(r));
}

double massDensityCorona(double r)
{
	double rtr = GlobalConfig.get<double>("rtr");
	return massDensityADAF(r)*(rt/r);
}

double electronDensity(double r)
{
	return massDensity(r)/(atomicMassUnit*eMeanMolecularWeight);
}

double ionDensity(double r)
{
	return massDensity(r)/(atomicMassUnit*iMeanMolecularWeight);
}