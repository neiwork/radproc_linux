#include "State.h"
//#include "functions.h"
#include "torusParameters.h"
#include "modelParameters.h"
#include <fparameters/Dimension.h>
#include <fmath/physics.h>
#include <fparameters/SpaceIterator.h>
#include <fparameters/parameters.h>
#include <boost/property_tree/ptree.hpp>
//#include <fmath/bisection.h>

State::State(boost::property_tree::ptree& cfg) :
 electron{ "electron" },
 photon{ "photon" },
 //proton("proton"),
 magf(photon.ps, false),
 tempElectrons(photon.ps, false),
 tempIons(photon.ps, false),
 tpf1(photon.ps, false),
 tpf2(photon.ps, false),
 denf_i(photon.ps, false),
 denf_e(photon.ps, false)
 {
	particles.push_back(&electron);
	particles.push_back(&photon);
	for (auto p : particles) {
		initializeParticle(*p, cfg);
	}
	magf.initialize();
	magf.fill([&](const SpaceIterator& i){
		static const double beta = GlobalConfig.get<double>("beta");
	    double r=i.val(DIM_R);
		double theta=i.val(DIM_THETA);
		return sqrt(beta*24.0*pi*pressureTot(r,theta));
    });
    denf_i.initialize();
    denf_i.fill([&](const SpaceIterator& i) {
        static const double mu=GlobalConfig.get<double>("mu_i");
        double r=i.val(DIM_R);
		double theta=i.val(DIM_THETA); 
		return energyDensity(r,theta)/(atomicMassUnit*mu);
    });
    denf_e.initialize();
    denf_e.fill([&](const SpaceIterator& i) {
        static const double mu=GlobalConfig.get<double>("mu_e");
        double r=i.val(DIM_R);
		double theta=i.val(DIM_THETA);
		return energyDensity(r,theta)/(atomicMassUnit*mu);
    });
    tempElectrons.initialize();
    tempElectrons.fill([&](const SpaceIterator& i) {
        double r=i.val(DIM_R);
		double theta=i.val(DIM_THETA);
		return temp_e(r, theta);
    });
    tempIons.initialize();
    tempIons.fill([&](const SpaceIterator& i) {
        double r=i.val(DIM_R);
		double theta=i.val(DIM_THETA);
		return temp_i(r, theta);
    });
}

Dimension* State::createDimension(Particle& p, std::string dimid, std::function<void(Vector&,double,double)> initializer, boost::property_tree::ptree& cfg) {
	int samples=p.getpar<int>(cfg,"dim."+dimid+".samples");
	double min=p.getpar<double>(cfg,"dim."+dimid+".min");
	double max=p.getpar<double>(cfg,"dim."+dimid+".max");
	return new Dimension(samples,bind(initializer,std::placeholders::_1,min,max));
}

void State::initializeParticle(Particle& p,boost::property_tree::ptree& cfg)
{
	using std::bind;

	p.configure(cfg.get_child("particle.default"));
	p.configure(cfg.get_child("particle."+p.id));
	p.ps.add(createDimension(p,"energy",initializeEnergyPoints,cfg));
	// we can't use createDimension because we're multiplying by pc before creating them
	
    // add dimension for R
    double rmin=GlobalConfig.get<double>("rCusp");
    double rmax=GlobalConfig.get<double>("rEdge");
	int nR=p.getpar(cfg,"dim.radius.samples",20);
	p.ps.add(new Dimension(nR,bind(initializeRadiiPoints,std::placeholders::_1,rmin,rmax)));
    
    // add dimension for theta
    double thetamin=p.getpar(cfg,"dim.theta.min",0.0);
    double thetamax=pi/2.0*p.getpar(cfg,"dim.theta.max",0.8);
    int nTheta=p.getpar(cfg,"dim.theta.samples",5);
    GlobalConfig.put("thetamin", GlobalConfig.get<double>("thetamin", thetamin));
    GlobalConfig.put("thetamax", GlobalConfig.get<double>("thetamax", thetamax));
    p.ps.add(new Dimension(nTheta,bind(initializeThetaPoints,std::placeholders::_1,thetamin,thetamax)));

	p.initialize();
}




