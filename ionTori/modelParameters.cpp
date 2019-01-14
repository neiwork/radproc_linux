#include "modelParameters.h"
#include "globalVariables.h"
#include <fparameters/parameters.h>
#include <fmath/physics.h>
#include <fmath/configure.h>
#include <iostream>
#include <algorithm>

void prepareGlobalCfg()
{
	torusParameters();
	fmath_configure(GlobalConfig);
}

void initEnergyPoints(Vector& v, double logEmin, double logEmax)
{
	double Emax = pow(10,logEmax)*EV_TO_ERG;
	double Emin =  pow(10,logEmin)*EV_TO_ERG;
	double E_int = pow(Emax/Emin,1.0/(v.size()-1));
    v[0]=Emin;
	for (size_t i=1;i<v.size();++i)
		v[i]=v[i-1]*E_int;
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