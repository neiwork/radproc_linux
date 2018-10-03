#include "modelParameters.h"
#include "torusParameters.h"
#include <fparameters/parameters.h>
#include <fmath/physics.h>
#include <fmath/configure.h>
#include <iostream>
#include <algorithm>

void prepareGlobalCfg()
{
    /*static const double massBH = GlobalConfig.get<double>("massBH");
    static const double spinBH = GlobalConfig.get<double>("spinBH") * massBH;
    static const double lambda = GlobalConfig.get<double>("lambda");
    static const double xi = GlobalConfig.get<double>("xi");
    static const double n = GlobalConfig.get<double>("n");
    static const double energyC = GlobalConfig.get<double>("energyC");
	static const double beta = GlobalConfig.get<double>("beta");*/
    
	static const double Mbh=GlobalConfig.get<double>("realMassBH")*solarMass;
	double rg=gravitationalConstant*Mbh/cLight2;
    GlobalConfig.put("rg", rg);
	
	//double rCusp,rCenter,rEdge;
	torusParameters();
	
	//GlobalConfig.put("l_0",GlobalConfig.get<double>("l_0",l_0));
    //GlobalConfig.put("rCusp",GlobalConfig.get<double>("rCusp",rCusp));
    //GlobalConfig.put("rCenter",GlobalConfig.get<double>("rCenter",rCenter));
    //GlobalConfig.put("rEdge",GlobalConfig.get<double>("rEdge",rEdge));
    
    static const double mu_i=GlobalConfig.get<double>("mu_i");
    static const double mu_e=GlobalConfig.get<double>("mu_e");
    static const double xi=GlobalConfig.get<double>("xi");
    double M_0=mu_i/(mu_e + mu_i);
    double M_1=mu_i*xi/(mu_e+mu_i*xi);
    GlobalConfig.put("M_0",GlobalConfig.get<double>("M_0",M_0));
    GlobalConfig.put("M_1",GlobalConfig.get<double>("M_1",M_1));
    static const double temp_ec=GlobalConfig.get<double>("temp_ec");
    static const double beta=GlobalConfig.get<double>("beta");
    static const double energyC=GlobalConfig.get<double>("energyC");
    double pK=boltzmann*temp_ec/((1.0-beta)*atomicMassUnit*pow(energyC,2.0/3.0)*mu_e*M_1);
    GlobalConfig.put("pK",GlobalConfig.get<double>("pK",pK));
    
	//GlobalConfig.put("Dlorentz", GlobalConfig.get<double>("Dlorentz", computeDlorentz(Gamma)));
	//DefOpt_IntLosses.samples_x = GlobalConfig.get<int>("integrate-losses.samples.x", DefOpt_IntLosses.samples_x);
	//DefOpt_IntLosses.samples_t = GlobalConfig.get<int>("integrate-losses.samples.t", DefOpt_IntLosses.samples_t);
	//DefOpt_IntLosses.samples_y = GlobalConfig.get<int>("integrate-losses.samples.y", DefOpt_IntLosses.samples_y);
	fmath_configure(GlobalConfig);
}

void initializeEnergyPoints(Vector& v, double logEmin, double logEmax)
{
	double Emax=1.6e-12*pow(10,logEmax);
	double Emin=1.6e-12*pow(10,logEmin);
	//double E_int=pow((50*Emax/Emin),1.0/(v.size()-1));
	double E_int=pow(Emax/Emin,1.0/(v.size()-1));
    v[0]=Emin;
	for (size_t i=1;i<v.size();++i) {
		v[i]=v[i-1]*E_int;
	}
}

void initializeRadiiPoints(Vector& v,double min,double max) 
{
    double dv=(max-min)/v.size();
    v[0]=min+dv/2.0;
    for (size_t i=1;i<v.size();++i) {
        v[i]=v[i-1]+dv;
    }
}

void initializeThetaPoints(Vector& v,double min,double max)
{
/*    double var_int=(sin(max)-sin(min))/(v.size()-1);
    v[0]=min;
    double sin0=sin(v[0]);
    for (size_t i=1;i<v.size();++i) {
        v[i]=asin(sin0+var_int);
        sin0+=var_int;
    }
*/
	double dv=(max-min)/v.size();
	v[0]=min;
	for (size_t i=1;i<v.size();++i) {
		v[i]=v[i-1]+dv;
	}
}