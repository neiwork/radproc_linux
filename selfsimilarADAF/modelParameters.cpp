#include "modelParameters.h"
#include <fparameters/parameters.h>
#include <fmath/physics.h>
#include <fmath/configure.h>
#include <iostream>
#include <algorithm>

void prepareGlobalCfg()
{
	adafParameters();
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

void initGridLogarithmically(Vector& v, double min, double max)
{
	double ratio = pow(max/min,1.0/(v.size()-1));
	v[0]=min;
	for (size_t i=1; i < v.size(); ++i) 
		v[i] = v[i-1]*ratio;
}

void initGridLinearly(Vector& v, double min, double max)
{
	double dv = (max-min)/(v.size()-1);
	v[0]=min;
	for (size_t i=1;i<v.size();++i)
		v[i] = v[i-1]+dv;
}