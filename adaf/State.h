#pragma once

#include <fparticle/Particle.h>
#include <fmath/mathematics.h>
#include <fmath/physics.h>
#include <boost/property_tree/ptree_fwd.hpp>

class State {
public:
	std::vector<Particle*> particles;

	Particle photon;
	Particle ntPhoton;
    Particle ntElectron;
    Particle ntProton;
    Particle ntNeutron;
	Particle ntChargedPion;
	Particle ntMuon;
	Particle neutrino;
	Particle ntPair;
    
	ParamSpaceValues magf;
	ParamSpaceValues denf_i;
    ParamSpaceValues denf_e;
    ParamSpaceValues tempElectrons, tempIons;
	ParamSpaceValues thetaH, height;
	ParamSpaceValues tau_gg;

	State(boost::property_tree::ptree& cfg);
	
	static Dimension* createDimension(Particle& p, string dimid, 
		function<void(Vector&,double,double)> initializer, function<double(double)> to_linear, function<double(double)> from_linear, boost::property_tree::ptree& cfg);

	static void initializeParticle(Particle& p, boost::property_tree::ptree& cfg);
};