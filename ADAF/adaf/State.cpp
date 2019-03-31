#include "State.h"
#include "globalVariables.h"
#include "torusFunctions.h"
#include "modelParameters.h"
#include <fparameters/Dimension.h>
#include <fmath/physics.h>
#include <fparameters/SpaceIterator.h>
#include <fparameters/parameters.h>
#include <boost/property_tree/ptree.hpp>
//#include <fmath/bisection.h>

using namespace std;

State::State(boost::property_tree::ptree& cfg) :
 electron{ "electron" },
 photon{ "photon" },
 proton("proton"),
 magf(photon.ps, false),
 tempElectrons(photon.ps, false),
 tempIons(photon.ps, false),
 denf_i(photon.ps, false),
 denf_e(photon.ps, false)
 {
	particles.push_back(&electron);
	particles.push_back(&photon);
	particles.push_back(&proton);
	for (auto p : particles) {
		initializeParticle(*p, cfg);
	}
	magf.initialize();
	magf.fill([&](const SpaceIterator& i){
	    double r=i.val(DIM_R);
		double theta=i.val(DIM_THETA);
		return sqrt(magFieldPar*24.0*pi*totalPressure(r,theta));
    });
    denf_i.initialize();
    denf_i.fill([&](const SpaceIterator& i) {
        double r=i.val(DIM_R);
		double theta=i.val(DIM_THETA); 
		return massDensity(r,theta)/(atomicMassUnit*iMeanMolecularWeight);
    });
    denf_e.initialize();
    denf_e.fill([&](const SpaceIterator& i) {
        double r=i.val(DIM_R);
		double theta=i.val(DIM_THETA);
		return massDensity(r,theta)/(atomicMassUnit*eMeanMolecularWeight);
    });
    tempElectrons.initialize();
    tempElectrons.fill([&](const SpaceIterator& i) {
        double r=i.val(DIM_R);
		double theta=i.val(DIM_THETA);
		return electronTemp(r,theta);
    });
    tempIons.initialize();
    tempIons.fill([&](const SpaceIterator& i) {
        double r=i.val(DIM_R);
		double theta=i.val(DIM_THETA);
		return ionTemp(r,theta);
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
	//p.ps.add(createDimension(p,"energy",initEnergyPoints,cfg));
	p.ps.add(new Dimension(nE,bind(initEnergyPoints,placeholders::_1,
								logMinEnergy,logMaxEnergy)));
	// add dimension for r
	double dr = (edgeRadius-cuspRadius)/nR;
	p.ps.add(new Dimension(nR,bind(initGridLinearly,placeholders::_1,
							cuspRadius+dr/2,edgeRadius-dr/2)));
    
    // add dimension for theta
    p.ps.add(new Dimension(nTheta,bind(initGridLinearly,std::placeholders::_1,
							minPolarAngle,maxPolarAngle)));
	p.initialize();
}
