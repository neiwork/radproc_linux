#include "State.h"
#include "globalVariables.h"
#include "adafFunctions.h"
#include "modelParameters.h"
#include <fparameters/Dimension.h>
#include <fmath/physics.h>
#include <fparameters/SpaceIterator.h>
#include <fparameters/parameters.h>
#include <boost/property_tree/ptree.hpp>

using namespace std;

State::State(boost::property_tree::ptree& cfg) :
 electron{ "electron" },
 proton{ "proton" },
 photon{ "photon" },
 ntProton("ntProton"),
 ntElectron{ "ntElectron" },
 magf(photon.ps, false),
 denf_i(photon.ps, false),
 denf_e(photon.ps, false),
 tempElectrons(photon.ps, false),
 tempIons(photon.ps, false),
 thetaH(photon.ps, false),
 height(photon.ps, false)
 {
	particles.push_back(&electron);
	particles.push_back(&proton);
	particles.push_back(&ntElectron);
	particles.push_back(&ntProton);
	particles.push_back(&photon);
	for (auto p : particles) {
		initializeParticle(*p, cfg);
	}
	magf.initialize();
	magf.fill([&](const SpaceIterator& i){
	    double r=i.val(DIM_R);
		return sqrt( (1.0-magFieldPar)*8.0*pi*massDensityADAF(r)*sqrdSoundVel(r) );
    });
    denf_i.initialize();
    denf_i.fill([&](const SpaceIterator& i) {
        double r=i.val(DIM_R); 
		return massDensityADAF(r)/(atomicMassUnit*iMeanMolecularWeight);
    });
    denf_e.initialize();
    denf_e.fill([&](const SpaceIterator& i) {
        double r=i.val(DIM_R);
		return massDensityADAF(r)/(atomicMassUnit*eMeanMolecularWeight);
    });
    tempElectrons.initialize();
    tempElectrons.fill([&](const SpaceIterator& i) {
        double r=i.val(DIM_R);
		return electronTemp(r);
    });
    tempIons.initialize();
    tempIons.fill([&](const SpaceIterator& i) {
        double r=i.val(DIM_R);
		return ionTemp(r);
    });
	thetaH.initialize();
	thetaH.fill([&](const SpaceIterator& i) {
		double r=i.val(DIM_R);
		return acos(costhetaH(r));
	});
	height.initialize();
	height.fill([&](const SpaceIterator& i) {
		double r=i.val(DIM_R);
		return height_fun(r);
	});
}

Dimension* State::createDimension(Particle& p, string dimid, 
		function<void(Vector&,double,double)> initializer, boost::property_tree::ptree& cfg)
{
	int samples=p.getpar<int>(cfg,"dim."+dimid+".samples");
	double min=p.getpar<double>(cfg,"dim."+dimid+".min");
	double max=p.getpar<double>(cfg,"dim."+dimid+".max");
	return new Dimension(samples,bind(initializer,placeholders::_1,min,max));
}

void State::initializeParticle(Particle& p,boost::property_tree::ptree& cfg)
{
	using std::bind;
	p.configure(cfg.get_child("particle.default"));
	p.configure(cfg.get_child("particle."+p.id));

	// add dimension for energies
	p.ps.add(createDimension(p,"energy",initEnergyPoints,cfg));
	//p.ps.add(new Dimension(nE,bind(initEnergyPoints,placeholders::_1,
	//							logMinEnergy,logMaxEnergy)));
	
	// add dimension for r
	double innerRadius = schwRadius*1.1;
	double edgeRadius = exp(logr.back())*schwRadius;
	p.ps.add(new Dimension(nR,bind(initGridLogarithmically,placeholders::_1,
							innerRadius*sqrt(paso_r),edgeRadius/sqrt(paso_r))));
							
	// add dimension for rcd
	p.ps.add(new Dimension(nRcd,bind(initGridLogarithmically,placeholders::_1,
							rTr*sqrt(paso_rCD),edgeRadius/sqrt(paso_rCD))));
	
	p.initialize();
}
