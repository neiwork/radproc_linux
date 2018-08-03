#include "modelParameters.h"
#include <fparameters/parameters.h>
#include <fmath/physics.h>
#include <fmath/configure.h>
#include <iostream>
#include <algorithm>

// Global variables
/*double massBH;                // Black hole mass
double RS;                    // Schwarzschild Radius
double alphapar;              // Viscosity parameter
double betapar;               // Magnetic parameter
double gammapar;              // Polytropic index
double rmin;                  // Inner edge of the ADAF
double rmax;                  // Outer edge of the ADAF
double f;                     // Advection parameter
double r;                     // Radius*/

void prepareGlobalCfg()
{
	static const double massBH=GlobalConfig.get<double>("massBH")*solarMass;
	double rg=gravitationalConstant*Mbh/cLight2;
    GlobalConfig.put("rg", rg);
	fmath_configure(GlobalConfig);
}

void initializeEnergyPoints(Vector& v, double logEmin, double logEmax)
{
	double Emax=1.6e-12*pow(10,logEmax);
	double Emin=1.6e-12*pow(10,logEmin);
	double E_int=pow(Emax/Emin,1.0/(v.size()-1));
    v[0]=Emin;
	for (size_t i=1;i<v.size();++i) {
		v[i]=v[i-1]*E_int;
	}
}

void initializeRadiiPoints(Vector& v,double min,double max) 
{
    double var_int=pow(max/min,1.0/(v.size()-1));
    v[0]=min;
    for (size_t i=1;i<v.size();++i) {
        v[i]=v[i-1]*var_int;
    }
}

void initializeThetaPoints(Vector& v,double min,double max)
{
    double var_int=(sin(max)-sin(min))/(v.size()-1);
    v[0]=min;
    double sin0=sin(v[0]);
    for (size_t i=1;i<v.size();++i) {
        v[i]=asin(sin0+var_int);
        sin0+=var_int;
    }
}