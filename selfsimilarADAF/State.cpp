#include "State.h"
#include "modelParameters.h"

#include <ssADAF/ssfunctions.h>
#include <ssADAF/ssADAF.h>
#include <ssADAF/packageData.h>

#include <fparameters/Dimension.h>
#include <fmath/physics.h>
#include <fparameters/SpaceIterator.h>
#include <fparameters/parameters.h>
#include <boost/property_tree/ptree.hpp>
#include <fmath/fbisection.h>
#include <nrMath/nrutil.h>

State::State(boost::property_tree::ptree& cfg) :
 electron{ "electron" },
 photon{ "photon" },
 magfield(photon.ps, false),
 tempElectrons(photon.ps, false),
 tempIons(photon.ps, false),
 denf_i(photon.ps, false),
 denf_e(photon.ps, false)
 {
	particles.push_back(&electron);
	particles.push_back(&photon);
	for (auto p : particles) {
		initializeParticle(*p, cfg);
	}
	magfield.initialize();
	denf_e.initialize();
	denf_i.initialize();
	tempElectrons.initialize();
	tempIons.initialize();
	
	double *x;
	x=dvector(1,3);
	
	//constants();
	dataADAF data;
	
	static const double massBH=GlobalConfig.get<double>("massBH")*solarMass;
	static const double f=GlobalConfig.get<double>("f");
	static const double alphapar=GlobalConfig.get<double>("alpha");
	static const double betapar=GlobalConfig.get<double>("beta");
	static const double rmin   = GlobalConfig.get<double>("rMin");
	static const double rmax   = GlobalConfig.get<double>("rMax");
	
	static const double mdot=GlobalConfig.get<double>("mdot");
	
	double gammapar=(32.0-24.0*betapar-3.0*betapar*betapar)/(24.0-21.0*betapar);
	
	GlobalConfig.put("gammapar", gammapar);

	data = (dataADAF){.massBH = massBH, .f = f, .alphapar = alphapar, .betapar  = betapar, .rmin = rmin, .rmax = rmax}; 
	//, .eDensity = eDensity, .iDensity = iDensity, .magField = magField, .eTemp = eTemp}; 

	
	x[1]=(boltzmann*1.0e12)/(protonMass*cLight2);
	x[2]=(boltzmann*1.0e9)/(electronMass*cLight2);
	x[3]=1.0e-3; //= mdot
	photon.ps.iterate([&](const SpaceIterator& iR) {
		double rB1=iR.val(DIM_R);
		if(iR.its[DIM_R].canPeek(1)) {
			
			double rB2=iR.its[DIM_R].peek(1);
			double r=sqrt(rB1*rB2);
			
			
			data.eDensity = ne(x[3], r, massBH, alphapar, gammapar, f);
			data.iDensity = ni(x[3], r, massBH, alphapar, gammapar, f);
			data.magField = magf(x[3], r, massBH, betapar, alphapar, gammapar, f);
			data.eTemp = x[1];

	
			adafSol(x, r, data);
			photon.ps.iterate([&](const SpaceIterator& iTh)  {  
				magfield.set(iTh,magf(x[3], r, massBH, betapar, alphapar, gammapar, f));
				denf_e.set(iTh,ne(x[3], r, massBH, alphapar, gammapar, f));
				denf_i.set(iTh,ni(x[3], r, massBH, alphapar, gammapar, f));
				tempElectrons.set(iTh,x[1]);
				tempIons.set(iTh,x[2]);
			},{0,iR.coord[DIM_R],-1});
		}
	},{0,-1,0});
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
    double rmin=GlobalConfig.get<double>("rMin");
    double rmax=GlobalConfig.get<double>("rMax");
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