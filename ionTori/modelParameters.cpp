#include "modelParameters.h"

#include "State.h"

#include <fmath/interpolation.h>
#include <fmath/RungeKutta.h>

#include <fparameters/parameters.h>
#include <fmath/physics.h>
#include <fmath/configure.h>
#include <iostream>
#include <algorithm>

// void fillMagnetic(State& st)
// {
	//static const double Gj = GlobalConfig.get<double>("Gamma");

//	st.magf.fill([&](const SpaceIterator& i) {
//		double z = i.val(DIM_R);
//		int z_ix = i.coord[DIM_R];

//		return computeMagField(z);  
//	});
//}*/

void prepareGlobalCfg()
{
    /*static const double massBH = GlobalConfig.get<double>("massBH");
    static const double spinBH = GlobalConfig.get<double>("spinBH") * massBH;
    static const double lambda = GlobalConfig.get<double>("lambda");
    static const double xi = GlobalConfig.get<double>("xi");
    static const double n = GlobalConfig.get<double>("n");
    static const double energyC = GlobalConfig.get<double>("energyC");
	static const double beta = GlobalConfig.get<double>("beta");*/
    
	
	/*GlobalConfig.put("massBH", GlobalConfig.get<double>("massBH"));
    GlobalConfig.put("spinBH", GlobalConfig.get<double>("spinBH") * massBH);
    GlobalConfig.put("lambda", GlobalConfig.get<double>("lambda"));
    GlobalConfig.put("xi", GlobalConfig.get<double>("xi"));
    GlobalConfig.put("n", GlobalConfig.get<double>("n"));
    GlobalConfig.put("energyC", GlobalConfig.get<double>("energyC"));
	GlobalConfig.put("beta", GlobalConfig.get<double>("beta"));*/
	
	//OJO con los nombres beta y n!!
	
	//GlobalConfig.put("Dlorentz", GlobalConfig.get<double>("Dlorentz", computeDlorentz(Gamma)));
	//DefOpt_IntLosses.samples_x = GlobalConfig.get<int>("integrate-losses.samples.x", DefOpt_IntLosses.samples_x);
	//DefOpt_IntLosses.samples_t = GlobalConfig.get<int>("integrate-losses.samples.t", DefOpt_IntLosses.samples_t);
	//DefOpt_IntLosses.samples_y = GlobalConfig.get<int>("integrate-losses.samples.y", DefOpt_IntLosses.samples_y);
	fmath_configure(GlobalConfig);
}

void initializeEnergyPoints(Vector& v, double logEmin, double logEmax)
{
	double Emax = 1.6e-12*pow(10, logEmax);
	double Emin = 1.6e-12*pow(10, logEmin);
	double E_int = pow((10 * Emax / Emin), (1.0 / (v.size() - 1)));
	v[0] = Emin;
	for (size_t i = 1; i < v.size(); ++i){
		v[i] = v[i - 1] * E_int;
	}
}

void initializeRPoints(Vector& v, double Rmin, double Rmax)
{

	double R_int = pow((Rmax / Rmin), (1.0 / (v.size() - 1.0)));

	v[0] = Rmin;

	for (size_t i = 1; i < v.size(); ++i){
		v[i] = v[i - 1] * R_int;
	}
}

void initializeThetaPoints(Vector& v, double Thetamin, double Thetamax)
{

	double Theta_int = pow((Thetamax / Thetamin), (1.0 / (v.size() - 1.0)));

	v[0] = Thetamin;

	for (size_t i = 1; i < v.size(); ++i){
		v[i] = v[i - 1] * Theta_int;
	}
}
