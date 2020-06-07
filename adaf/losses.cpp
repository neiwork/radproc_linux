#include "losses.h"
#include "globalVariables.h"
#include <fmath/physics.h>
#include <fmath/RungeKutta.h>
#include <gsl/gsl_math.h>

//#include "targetFields.h"

//#include <flosses/nonThermalLosses.h>
#include <flosses/lossesSyn.h>
#include <flosses/lossesIC.h>
#include <flosses/lossesHadronics.h>
#include <flosses/lossesPhotoHadronic.h>
#include <fparameters/parameters.h>
//#include <fparticle/Particle.h>


//#include <iostream>
//#include <map>

double losses(double E, Particle& p, State& st, const SpaceCoord& i)
{
	double B = st.magf.get(i);
	double density = st.denf_e.get(i)+st.denf_i.get(i);
	double losses = 0.0;
	if (p.id == "ntProton" || p.id == "ntChargedPion")
		losses = lossesSyn(E,B,p) + lossesHadronics(E,density,p)
			+ lossesPhotoHadronic_simple(E,p,st.photon.distribution,i,st.photon.emin(),st.photon.emax());
	else if(p.id == "ntElectron" || p.id == "ntPair" || p.id == "ntMuon") {
		double lossSy = lossesSyn(E,B,p);
		double lossIC1, lossIC2;
		lossIC1 = lossIC2 = 0.0;
		int comptonLossesImportant = GlobalConfig.get<int>("nonThermal.comptonLosses");
		if (comptonLossesImportant) {
			lossIC1 = lossesIC(E,p,st.photon.distribution,i,st.photon.emin(),st.photon.emax());
			lossIC2 = lossesIC(E,p,st.ntPhoton.distribution,i,st.ntPhoton.emin(),st.ntPhoton.emax());
		}
		losses = lossSy+lossIC1+lossIC2;
	}
	return losses;
}

#include "adafFunctions.h"
#include <fparameters/Dimension.h>

double b(double E, double r, Particle& p, State& st, const SpaceCoord& psc) {

	double B = (r < p.ps[DIM_R].last()) ? st.magf.interpolate({{1,r}},&psc) : 0.0;
	double density = (r < p.ps[DIM_R].last()) ? st.denf_e.interpolate({{1,r}},&psc) : 0.0;
	double losses = 0.0;

	if(p.id == "ntProton"){
			losses = lossesSyn(E, B, p) + lossesHadronics(E, density, p);// + lossesPhotoHadronic(E, p, st.photon.distribution, i, st.photon.emin(), st.photon.emax());
	}
	else if(p.id == "ntElectron" || p.id == "ntPair"){	
			losses = lossesSyn(E, B, p);//  + lossesIC(E, p, st.photon.distribution, i, st.photon.emin(), st.photon.emax());
	}
	return -losses;
}

double t_cool(double E, double r, Particle& p) {

	double B = (r < p.ps[DIM_R].last()) ? magneticField(r) : 0.0;
	double density = (r < p.ps[DIM_R].last()) ? ionDensity(r) : 0.0;
	double losses = 0.0;

	if(p.id == "ntProton"){
			losses = lossesSyn(E, B, p) + lossesHadronics(E, density, p);// + lossesPhotoHadronic(E, p, st.photon.distribution, i, st.photon.emin(), st.photon.emax());
	}
	else if(p.id == "ntElectron" || p.id == "ntPair"){	
			losses = lossesSyn(E, B, p);//  + lossesIC(E, p, st.photon.distribution, i, st.photon.emin(), st.photon.emax());
	}
	return E/losses;
}

double annihilationRate(double Ep, Particle& p, const SpaceCoord& psc)
{
	double gamma_p = Ep / electronRestEnergy;
	double integral = integSimpsonLog(p.emin(),p.emax(),[gamma_p,&p,&psc](double Em)
	{
		double gamma_m = Em / electronRestEnergy;
		double Nm = (Em > p.emin() && Em < p.emax()) ?
				p.distribution.interpolate({{DIM_E,Em}},&psc) : 0.0;
		double result = Nm / gamma_m * (log(4.0*gamma_m*gamma_p)-2.0);
		return (result > 0.0 ? result : 0.0);
	},80);
	return pi*gsl_pow_2(electronRadius)*cLight / gamma_p * integral;
}