#include "torusFunctions.h"
#include "globalVariables.h"
#include <fmath/physics.h>
#include "metric.h"

double massDensity(double r, double theta)  // = mass density (c=1)
{
	double w;
	double result = (w = normalizedPotential(r,theta)) > 0.0 ? 
		pow((pow(polytropConst*pow(centralMassDens,1.0/polytropIndex) + 
		1.0,w)-1.0)/polytropConst,polytropIndex) : 0.0;

	return result;
}

double electronDensity(double r, double theta)
{
	double result = massDensity(r,theta)/(atomicMassUnit*eMeanMolecularWeight);
	return result;
}

double totalPressure(double r, double theta) 
{
    return polytropConst*pow(massDensity(r, theta),1.0+1.0/polytropIndex);
}

/////////////////////////////////////////////
// TEMPERATURES

// Electrons
double electronTemp(double r, double theta) {
	double w;
    double result = (w = normalizedPotential(r,theta)) > 0.0 ? 
		(1.0-w)*auxM0 + w*auxM1*eMeanMolecularWeight* 
		((1.0-magFieldPar)*atomicMassUnit*totalPressure(r,theta))
		/(boltzmann*massDensity(r,theta))+2.7 : 2.7;
	return result;
}

// Ions
double ionTemp(double r, double theta) {
	double w;    
    double result = (w = normalizedPotential(r,theta)) > 0.0 ? 
		(eMeanMolecularWeight/iMeanMolecularWeight*auxM0 + w*(auxM0-auxM1))*iMeanMolecularWeight*
		( (1.0-magFieldPar)*atomicMassUnit*totalPressure(r,theta) ) / 
		(boltzmann*massDensity(r,theta))+2.7 : 2.7;
	return result;
}


// IMPLEMENTATION OF INTERNAL FUNCTIONS

void marginalOrbits(double& r_ms, double& r_mb) //calculates r_ms and r_mb
{
	double z1 = 1.0 + pow( 1.0 - blackHoleSpinPar*blackHoleSpinPar, 1.0/3.0) * 
		( pow(1.0 + blackHoleSpinPar, 1.0/3.0) + pow(1.0 - blackHoleSpinPar, 1.0/3.0) );
	double z2 = sqrt( 3.0 * blackHoleSpinPar*blackHoleSpinPar + z1*z1);

	// marginally stable circular orbit
	r_ms = 3.0 + z2 - sqrt( (3.0 - z1) * (3.0 + z1 + 2.0*z2) );
	// marginally bound circular orbit
	r_mb = 2.0 - blackHoleSpinPar + 2.0 * sqrt(1.0-blackHoleSpinPar);
}

// Keplerian specific angular momentum
double keplAngularMom(double r)
{
    return (r*r-2.0*blackHoleSpinPar*sqrt(r)+blackHoleSpinPar*blackHoleSpinPar)/
			(pow(r, 1.5)-2.0*sqrt(r)+blackHoleSpinPar);
}

// ANGULAR VELOCITIY OF THE TORUS
double torusAngularVel(double r, double theta)
{
	return -(g_tphi(r,theta)+specificAngMom*g_tt(r,theta)) / 
			(g_phiphi(r,theta)+specificAngMom*g_tphi(r,theta));
}

// REDSHIFT FACTOR
double redshiftFactor(double r, double theta)
{
	double torusAngVel = torusAngularVel(r,theta);
	double aux = g_tt(r,0.0)+2.0*torusAngVel*g_tphi(r,0.0)+g_phiphi(r,0.0)*torusAngVel*torusAngVel;
	return sqrt(-aux);
}

// POTENTIAL FUNCTION
double gravPotential(double r, double theta) {
	double angVel = torusAngularVel(r,theta);
	double gtt=g_tt(r,theta);
	double gtphi=g_tphi(r,theta);
	double gphiphi=g_phiphi(r,theta);
	double aux = gtt + 2.0*angVel*gtphi + gphiphi*angVel*angVel;
	return (aux < 0.0) ? 0.5 * log(-aux / ((gtt+angVel*gtphi)*(gtt+angVel*gtphi))) : 0.0;
}

// NORMALIZED POTENTIAL FUNCTION
double normalizedPotential(double r, double theta)
{
	double potentialS = gravPotential(cuspRadius, 0.0);        		// potential at the torus surface
	double potentialC = gravPotential(torusCenterRadius, 0.0);      // potential at the torus center
	
	return (r > cuspRadius) ? (gravPotential(r,theta)-potentialS)/(potentialC-potentialS) : -1.0;
}

double modfKepl(double r) {
	return keplAngularMom(r) - specificAngMom; }
double modfw(double r)
	{return normalizedPotential(r,0.0);}