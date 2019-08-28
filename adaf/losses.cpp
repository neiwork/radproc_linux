#include "losses.h"
#include "globalVariables.h"
#include <fmath/physics.h>

//#include "targetFields.h"

//#include <flosses/nonThermalLosses.h>
#include <flosses/lossesSyn.h>
#include <flosses/lossesIC.h>
#include <flosses/lossesHadronics.h>
#include <flosses/lossesPhotoHadronic.h>

//#include <fparticle/Particle.h>


//#include <iostream>
//#include <map>

double losses(double E, Particle& p, State& st, const SpaceCoord& i)
{
	double B = st.magf.get(i);
	double density = st.denf_e.get(i)+st.denf_i.get(i);
	
	double losses = 0.0;

	if(p.id == "ntProton"){
			losses = lossesSyn(E,B,p) + lossesHadronics(E,density,p);// + lossesPhotoHadronic(E, p, st.photon.distribution, i, st.photon.emin(), st.photon.emax());
	}
	else if(p.id == "ntElectron" || p.id == "ntPair"){	
			losses = lossesSyn(E, B, p);//+lossesIC(E,p,st.photon.distribution,i,st.photon.emin(),st.photon.emax());
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
