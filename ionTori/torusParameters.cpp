#include "torusParameters.h"
#include "metric.h"


#include <fmath/physics.h>
extern "C" {
#include <nrMath/bisection.h>
}

#include <fparameters/parameters.h>
#include <boost/property_tree/ptree.hpp>
//#include "auxFunctions.h"
//#include <math.h>

#include <math.h>
#define ALLERR 1.0e-5

/*
double rCenter, rMin, l_0;
double n = 1.5;                  // Polytropic index
double temp_ec = 1.0e11;         // Electron temperature at the torus center
double xi = 1.0e-3;              // Electron over Ion temperature
double mu_i = 1.23;              // Mean molecular weight for ions
double mu_e = 1.14;              // Mean molecular weight for electrons
double energyC = 1.e-15;         // Energy density at the torus center
double beta = 0.5;               // Magnetic over total pressure
double pK;
*/


void marginalOrbits(double &r_ms, double &r_mb) //calculates r_ms and r_mb
{
	static const double massBH=GlobalConfig.get<double>("massBH");
    static const double spinBH=GlobalConfig.get<double>("spinBH") * massBH;
	
	double z1 = 1.0 + pow( 1.0 - (spinBH/massBH)*(spinBH/massBH), 1.0/3.0) * 
	( pow(1.0 + spinBH/massBH, 1.0/3.0) + pow(1.0 - spinBH/massBH, 1.0/3.0) );
	double z2 = sqrt( 3.0 * (spinBH/massBH)*(spinBH/massBH) + z1*z1);

	// marginally stable circular orbit
	r_ms = massBH * (3.0 + z2 - sqrt( (3.0 - z1) * (3.0 + z1 + 2.0*z2) ) );
	// marginally bound circular orbit
	r_mb = 2.0 * massBH - spinBH + 2.0 * sqrt(massBH) * sqrt(massBH-spinBH);
}
//////////

//now the functions to obtain rCusp and rCenter: criticalRadii

// Keplerian specific angular momentum
double keplAngularMom(double r) {
	static const double massBH=GlobalConfig.get<double>("massBH");
    static const double spinBH=GlobalConfig.get<double>("spinBH") * massBH;
    
    return sqrt(massBH)*(r*r-2.0*spinBH*sqrt(massBH*r)+spinBH*spinBH)/
               (pow(r, 1.5)-2.0*massBH*sqrt(r)+spinBH*sqrt(massBH));
}

double modfKepl(double r)
{
  static const double l_0=GlobalConfig.get<double>("l_0");
  return keplAngularMom(r) - l_0;
}

void criticalRadii(double &rCusp, double &rCenter, double &r_ms, double &r_mb, double lambda)
{
  //double r_ms, r_mb;
  marginalOrbits(r_ms, r_mb);

  int maxmitr = 1000;
  double allerr = 1.0e-3;
  rCusp = bisection(r_mb, r_ms, allerr, maxmitr, modfKepl);
  rCenter = bisection(r_ms, 500.0, allerr, maxmitr, modfKepl);
}
///////////////////////

//the following functions are to obtain rEdge: function edge

// ANGULAR VELOCITIY OF THE TORUS
double angularVel(double r, double theta)
{  
	static const double l_0=GlobalConfig.get<double>("l_0");
	return - ( g_tphi(r, theta) + l_0 * g_tt(r,theta) ) / ( g_phiphi(r, theta) + l_0 * g_tphi(r, theta) );
}

// POTENTIAL FUNCTION
double potential(double r, double theta) {
  double aux = g_tt(r,theta) + 2.0*angularVel(r, theta)*g_tphi(r, theta) + 
    g_phiphi(r, theta) * angularVel(r, theta) * angularVel(r, theta);
  return (aux < 0.0) ? 0.5 * log(-aux / (g_tt(r, theta)+angularVel(r, theta)*g_tphi(r, theta)) / 
				 (g_tt(r, theta)+angularVel(r, theta)*g_tphi(r, theta))) : 0.0;
}

// NORMALIZED POTENTIAL FUNCTION
double w(double r, double theta)
{ 
	static const double rMin=GlobalConfig.get<double>("model.particle.default.dim.radius.min");
	static const double rCenter=GlobalConfig.get<double>("rCenter");
	
	double potentialS = potential(rMin, 0.0);         // potential at the torus surface
	double potentialC = potential(rCenter, 0.0);      // potential at the torus center

  return (potential(r, theta) - potentialS) / (potentialC - potentialS);
}

double modfw(double r)
{
  return w(r,0.0);
}

double edge(double &rCenter)
{
  int maxmitr = 100000;
  double allerr = 1.0e-5;
  return bisection(rCenter, 500.0*rCenter, allerr, maxmitr, modfw);
}




double specificAngularMom(double r_ms, double r_mb, double lambda)
{
  //double r_ms, r_mb;
  marginalOrbits(r_ms, r_mb);
  
  double l_ms = keplAngularMom(r_ms);             // Keplerian specific angular momentum at r = r_ms
  double l_mb = keplAngularMom(r_mb);             // Keplerian specific angular momentum at r = r_mb
  
  return (1.0 - lambda) * l_ms + lambda * l_mb;
}

// AUXILIARY FUNCTIONS /////////////////////////////






//////////////////////////////////////


// Torus Parameters
void torusParameters() 
{
	
	static const double massBH=GlobalConfig.get<double>("massBH");
    static const double spinBH=GlobalConfig.get<double>("spinBH")*massBH;
	static const double lambda=GlobalConfig.get<double>("lambda");
	
	double r_ms, r_mb, rCusp, rCenter;
	marginalOrbits(r_ms, r_mb); //returns r_ms and r_mb
	
	double l_0 = specificAngularMom(r_ms,r_mb,lambda);
	GlobalConfig.put("l_0",GlobalConfig.get<double>("l_0",l_0));
	
	criticalRadii(rCusp, rCenter, r_ms, r_mb, lambda); //returns  rCusp and rCenter
	GlobalConfig.put("rCusp",GlobalConfig.get<double>("rCusp",rCusp));
    GlobalConfig.put("rCenter",GlobalConfig.get<double>("rCenter",rCenter));
	
	
    double rEdge=edge(rCenter);
	GlobalConfig.put("rEdge",GlobalConfig.get<double>("rEdge",rEdge));
	
	
	
	/*
	// Auxiliary variables
    double z1=1.0+pow(1.0-(spinBH/massBH)*(spinBH/massBH),1.0/3.0)* 
    (pow(1.0+spinBH/massBH,1.0/3.0)+pow(1.0-spinBH/massBH,1.0/3.0));
    double z2=sqrt(3.0*(spinBH/massBH)*(spinBH/massBH)+z1*z1);
    // marginally stable circular orbit
    double r_ms=massBH*(3.0+z2-sqrt((3.0-z1)*(3.0+z1+2.0*z2)));
    // marginally bound circular orbit
    double r_mb=2.0*massBH-spinBH+2.0*sqrt(massBH)*sqrt(massBH-spinBH);
    
	//double l_ms=keplAngularMom(r_ms);              // Keplerian specific angular momentum at r = r_ms
    //double l_mb=keplAngularMom(r_mb);             // Keplerian specific angular momentum at r = r_mb
    
    //l_0=(1.0-lambda)*l_ms+lambda*l_mb;
 
	//rCusp=bisection([&l_0](double r){return keplAngularMom(r)-l_0;},r_mb,r_ms,ALLERR);
	//rCenter=bisection(([&l_0](double r){return keplAngularMom(r)-l_0;},r_ms,100*rCusp,ALLERR);
    rEdge=bisection([&](double r){return w(r,0.0);},rCenter,100*rCenter,ALLERR);
	*/
	
}

// energyDensity
double energyDensity(double r, double theta)  // = mass density (c=1)
{
	static const double M_1 =GlobalConfig.get<double>("M_1");
    static const double pK  = GlobalConfig.get<double>("pK");
	static const double energyC = GlobalConfig.get<double>("energyC");
	static const double n = GlobalConfig.get<double>("n");
	
 // double M_1 = mu_i * xi / (mu_e + mu_i * xi);
  //pK = boltzmann * temp_ec / ( (1.0 - beta) * atomicMassUnit * pow(energyC, 2.0/3.0) * mu_e * M_1 );
  return ( w(r, theta) > 0.0 ) ? pow( (pow(pK*pow(energyC, 1.0/n) + 1.0, w(r, theta)) - 1.0) / pK, n) : 0.0;
}

double electronDensity(double r, double theta)
{
	static const double mu_e = GlobalConfig.get<double>("mu_e");
	return energyDensity(r, theta) / (atomicMassUnit*mu_e);
}


// PRESSURE
double pressureTot (double r, double theta) {
    static const double pK = GlobalConfig.get<double>("pK");
    static const double n = GlobalConfig.get<double>("n");

    return pK * pow(energyDensity(r, theta), 1.0 + 1.0/n);
}

/////////////////////////////////////////////
// TEMPERATURES

// Electrons
double temp_e(double r, double theta) {
    
    static const double M_0 = GlobalConfig.get<double>("M_0");
    static const double M_1 = GlobalConfig.get<double>("M_1");
    static const double mu_e = GlobalConfig.get<double>("mu_e");
    static const double beta = GlobalConfig.get<double>("beta");
    
    //Falta definir mu en constantes fundamentales (tambien podría ser mu_e y mu_i)

    return (w(r, theta) > 0.0) ? (1.0 - w(r,theta)) * M_0 + w(r,theta) * M_1 *
                mu_e * ( (1.0-beta)*atomicMassUnit*pressureTot(r,theta) ) / ( boltzmann * energyDensity(r, theta) ) + 2.7 : 2.7;
}

// Ions
double temp_i(double r, double theta) {
    
    static const double M_0 = GlobalConfig.get<double>("M_0");
    static const double M_1 = GlobalConfig.get<double>("M_1");
    static const double mu_i = GlobalConfig.get<double>("mu_i");
    static const double beta = GlobalConfig.get<double>("beta");
    
    return (w(r, theta) > 0.0) ? ( (1 - w(r, theta))*M_0 + w(r, theta)*M_1 ) *
                mu_i * ( (1.0-beta)*atomicMassUnit*pressureTot(r, theta) ) / ( boltzmann * energyDensity(r, theta) ) + 2.7: 2.7;
}

