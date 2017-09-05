// FUNCTIONS
#include "functions.h"
#include "modelParameters.h"
#include <fmath/physics.h>
#include <fmath/bisection.h>

#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>



void readSpinMbh(double massBH, double spinBH)
{	
	massBH = GlobalConfig.get<double>("massBH");
    spinBH = GlobalConfig.get<double>("spinBH") * massBH;
	
}


// Keplerian specific angular momentum
double keplAngularMom(double r) {
	
	double massBH, spinBH;
	readSpinMbh(massBH, spinBH);
	
    return sqrt(massBH) * ( r*r - 2.0 * spinBH * sqrt(massBH*r) + spinBH*spinBH ) /
               (pow(r, 1.5) - 2.0 * massBH * sqrt(r) + spinBH* sqrt(massBH) );
}

// Torus Parameters
void torusParameters(double *l_0, double *rCusp, double *rCenter) {
    
	double massBH, spinBH;
	readSpinMbh(massBH, spinBH);
	
	// Auxiliary variables
	
    double z1 = 1.0 + pow( 1.0 - (spinBH/massBH)*(spinBH/massBH), 1.0/3.0) * 
    ( pow(1.0 + spinBH/massBH, 1.0/3.0) + pow(1.0 - spinBH/massBH, 1.0/3.0) );
    double z2 = sqrt( 3.0 * (spinBH/massBH)*(spinBH/massBH) + z1*z1);

    // marginally stable circular orbit
    double r_ms = massBH * (3.0 + z2 - sqrt( (3.0 - z1) * (3.0 + z1 + 2.0*z2) ) );

    // marginally bound circular orbit
    double r_mb = 2.0 * massBH - spinBH + 2.0 * sqrt(massBH) * sqrt(massBH-spinBH);

    double l_ms = keplAngularMom(r_ms);              // Keplerian specific angular momentum at r = r_ms
    double l_mb = keplAngularMom(r_mb);              // Keplerian specific angular momentum at r = r_mb
    
    *l_0 = (1.0 - lambda) * l_ms + lambda * l_mb;
    *rCusp = bisection(r_mb, r_ms, angular_mom(r) - (*l_0));
    *rCenter = bisection(r_ms, 1000.0, angular_mom(r) - (*l_0));
}

/////////////////////////////////////////////////////////////////////////////////////
// METRIC COMPONENTS (in Boyer-Lindquist coordinates)

double g_tt(double r, double theta) {
	
	double massBH, spinBH;
	readSpinMbh(massBH, spinBH);
	
    double sigma = r*r + spinBH*spinBH*sin(theta)*sin(theta);
    return - (1.0 - 2.0*massBH*r / sigma);
}

double g_rr(double r, double theta) {
	
	double massBH, spinBH;
	readSpinMbh(massBH, spinBH);
	
    double delta = r*r - 2.0 * massBH * r + spinBH*spinBH;
    double sigma = r*r + spinBH*spinBH*sin(theta)*sin(theta);
    return sigma / delta;
}

double g_thetatheta(double r, double theta) {    // = Sigma
	double massBH, spinBH;
	readSpinMbh(massBH, spinBH);

    return r*r + spinBH*spinBH*sin(theta)*sin(theta);
}

double g_tphi(double r, double theta) {
	
	double massBH, spinBH;
	readSpinMbh(massBH, spinBH);
	
    double sigma = r*r + spinBH*spinBH*sin(theta)*sin(theta);
    return -(2.0*massBH*r*spinBH / sigma) * cos(theta)*cos(theta);
}

double g_phiphi(double r, double theta) {

	double massBH, spinBH;
	readSpinMbh(massBH, spinBH);
	
    double sigma = r*r + spinBH*spinBH*sin(theta)*sin(theta);
    return ( r*r + spinBH*spinBH + 2.0*massBH*r*spinBH*spinBH*cos(theta)*cos(theta) / sigma )
    * sin(theta)*sin(theta);
}

/////////////////////////////////////////////////////////////////////////////////////

// ANGULAR VELOCITIY OF THE TORUS
double angularVel(double r, double theta)  {
  return - ( g_tphi(r, theta) + (*l_0) * g_tt(r,theta) ) / ( g_phiphi(r, theta) + * (*l_0) * g_phiphi(r, theta) );
}

// POTENTIAL FUNCTION
double potential(double r, double theta) {
    double aux = g_tt(r,theta) + 2.0*angularVel(r,theta)*g_tphi(r,theta) + 
    g_phiphi(r,theta) * P2(angularVel(r,theta));
    return (aux < 0.0) ? 0.5 * log(-aux / P2(g_tt(r,theta)+angularVel(r,theta)*g_tphi(r,theta))) : 0.0;
}

// NORMALIZED POTENTIAL FUNCTION
double w(double r, double theta) {
    double potentialS = potential(*rCusp, 0.0)         // potential at the torus surface
    double potentialC = potential(*rCenter, 0.0)      // potential at the torus center
    return (potential(r,theta) - potentialS) / (potentialC - potentialS);
}

// ENERGY DENSITY
double energyDensity (double r, double theta) {
    return ( w(r, theta) > 0 ) ? pow(pK, -n) * pow( pow(pK*pow(energyC, 1.0/n) + 1.0, w(r, theta)) - 1.0, n) : 0.0;
}

// PRESSURE
double pressureTot (double r, double theta) {
    return pK * pow(energyDensity(r, theta), 1.0 + 1.0/n);
}

////////////////////////////////////////////////////////////////////////////////////////////
// TEMPERATURES

// Electrons
double temp_e(double r, double theta) {
    return (w(r, theta) > 0.0) ? (1.0 - w(r,theta)) * M_0 + w(r,theta) * M_1 *
                mu_e * ( (1.0-beta)*mu*pressureTot(r,theta) ) / ( boltzmann * energyDensity(r,theta) ) : 0.0;
}

// Ions
double temp_i(double r, double theta) {
    return (w(r, theta) > 0.0) ? ( (1 - w(r,theta))*M_0 + w(r,theta)*M_1 ) *
                mu_i * ( (1.0-beta)*mu*pressureTot(r,theta) ) / ( boltzmann * energyDensity(r,theta) ) : 0.0;
}
///////////////////////////////////////////////////////////////////////////////////////////
