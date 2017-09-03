#include "modelParameters.h"

#include "State.h"


#include <fmath/interpolation.h>
#include <fmath/RungeKutta.h>

#include <fparameters/parameters.h>
#include <fmath/physics.h>
#include <fmath/configure.h>
#include <iostream>
#include <algorithm>




double computeMagField(double z) {

	//static const double Lj = GlobalConfig.get<double>("Lt");

	return  (5.0e16/z);  

} 


double blackBody(double E, double z)
{

	//static const double starT = GlobalConfig.get<double>("starT");
	
	double Epeak = boltzmann*1.0e10; //VER starT;
	double Ephmin = Epeak / 1000.0;
	double Ephmax = Epeak*1000.0;

	//static const double int_E = RungeKuttaSimple(Ephmin, Ephmax, [&Ephmin, &Epeak](double E) {
	//	return (P3(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));
	//});

	double normalizacion = 8.0*pi/P3(planck*cLight); //VER  wph(z) / int_E;

	double nph = normalizacion*(P2(E)*exp(-Ephmin / E) / (exp(E / Epeak) - 1));

	return nph;
	
}


/*void fillPhotonField(State& st)
{

	st.tpf.fill([&](const SpaceIterator& i) {
		double z = i.val(DIM_R);
		double E = i.val(DIM_E);

		return blackBody(E, z);  
	});
}



void fillMagnetic(State& st)
{
	//static const double Gj = GlobalConfig.get<double>("Gamma");

	st.magf.fill([&](const SpaceIterator& i) {
		double z = i.val(DIM_R);
		int z_ix = i.coord[DIM_R];

		return computeMagField(z);  
	});
}*/




double eEmax(double z, double B)
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");
	static const double Gamma = GlobalConfig.get<double>("Gamma");
	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency");

	double Reff = 10.0*z;
	double vel_lat = cLight*openingAngle;

	double Emax_ad = accEfficiency*3.0*(z*openingAngle)*cLight*electronCharge*B / (vel_lat*Gamma); //
	double Emax_syn = electronMass*cLight2*sqrt(accEfficiency*6.0*pi*electronCharge / (thomson*B));
	
	double ampl = Gamma; //factor de amplificaci'on de B en la zona del choque
	double Emax_diff = electronCharge*B*Reff*sqrt(3.0*accEfficiency*ampl/2.0);
	double min1 = std::min(Emax_syn, Emax_ad);
	double min2 = std::min(min1, Emax_diff);


	return min2;
		
}




void prepareGlobalCfg()
{
	//static const double Gamma = GlobalConfig.get<double>("Gamma", 10);
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

//VER
void initializeRPoints(Vector& v, double Rmin, double Rmax)
{

	double R_int = pow((Rmax / Rmin), (1.0 / (v.size() - 1)));

	v[0] = Rmin;

	for (size_t i = 1; i < v.size(); ++i){
		v[i] = v[i - 1] * R_int;
	}

}

void initializeCrossingTimePoints(Vector& time, double rMin, double rMax)
{
	double R_int = pow((rMax / rMin), (1.0 / (time.size()-1)));

	Vector v(time.size()+1, 0.0);

	v[0] = rMin;

	for (size_t i = 0; i < time.size(); ++i){
		v[i+1] = v[i] * R_int;
		double delta_x = v[i+1]-v[i];
		time[i] = delta_x / cLight;  //construyo los t_i como los crossing time de las celdas i
	}

}

