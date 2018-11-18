#include "metric.h"
#include <math.h>
#include <fparameters/parameters.h>
#include <boost/property_tree/ptree.hpp>

extern const double blackHoleSpinPar;

double g_tt(double r, double theta)
{  
	double a2 = blackHoleSpinPar*blackHoleSpinPar;
	double sigma = r*r + a2*sin(theta)*sin(theta);
	return - (1.0 - 2.0*r / sigma);
}

double g_rr(double r, double theta)
{  	
	double a2 = blackHoleSpinPar*blackHoleSpinPar;
	double delta = r*r - 2.0 * r + a2;
	double sigma = r*r + a2*sin(theta)*sin(theta);
	return sigma / delta;
}

double g_thetatheta(double r, double theta)
{  
	double a2 = blackHoleSpinPar*blackHoleSpinPar;
	return r*r + a2*sin(theta)*sin(theta);
}

double g_tphi(double r, double theta)
{
	double a2 = blackHoleSpinPar*blackHoleSpinPar;
	double sigma = r*r + a2*sin(theta)*sin(theta);
	return - 2.0*r*blackHoleSpinPar / sigma * cos(theta)*cos(theta);
}

double g_phiphi(double r, double theta)
{
	double a2 = blackHoleSpinPar*blackHoleSpinPar;
	double sigma = r*r + a2*sin(theta)*sin(theta);
	return (r*r + a2 + 2.0*r*a2*cos(theta)*cos(theta) / sigma )
    * cos(theta)*cos(theta);
}
