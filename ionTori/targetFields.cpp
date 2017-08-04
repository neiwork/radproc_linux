#include "targetFields.h"

#include "modelParameters.h"

#include <fparameters/SpaceIterator.h>
#include <fmath/RungeKutta.h>
#include <fmath/physics.h>
#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>



void fillPhotonField(State& st)
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

