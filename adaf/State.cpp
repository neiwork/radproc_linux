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

State::State(boost::property_tree::ptree& cfg):
    photon{ "photon" },
	ntPhoton{ "ntPhoton" },
    ntProton{ "ntProton" },
    ntElectron{ "ntElectron" },
    ntNeutron{ "ntNeutron" },
	ntChargedPion{ "ntChargedPion" },
	ntMuon{ "ntMuon" },
	neutrino{ "neutrino" },
	ntPair{ "ntPair" },
	tau_gg(ntPhoton.ps, false),
    magf(ntPhoton.ps, false),
    denf_i(ntPhoton.ps, false),
    denf_e(ntPhoton.ps, false),
    tempElectrons(ntPhoton.ps, false),
    tempIons(ntPhoton.ps, false),
    thetaH(ntPhoton.ps, false),
    height(ntPhoton.ps, false)
    {
		particles.push_back(&photon);
		particles.push_back(&ntPhoton);
        particles.push_back(&ntElectron);
        particles.push_back(&ntProton);
        particles.push_back(&ntNeutron);
		particles.push_back(&ntChargedPion);
		particles.push_back(&ntMuon);
		particles.push_back(&neutrino);
		particles.push_back(&ntPair);
        for (auto p : particles) {
            initializeParticle(*p, cfg);
        }
        magf.initialize();
        magf.fill([&](const SpaceIterator& i){
            double r = i.val(DIM_R);
            return magneticField(r);
        });
        denf_i.initialize();
        denf_i.fill([&](const SpaceIterator& i) {
            double r=i.val(DIM_R); 
            double dens = massDensityADAF(r)/(atomicMassUnit*iMeanMolecularWeight);
			return dens;
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
		tau_gg.initialize();
}

Dimension* State::createDimension(Particle& p, string dimid, 
		function<void(Vector&,double,double)> initializer, function<double(double)> to_linear,  function<double(double)> from_linear, boost::property_tree::ptree& cfg)
{
	int samples=p.getpar<int>(cfg,"dim."+dimid+".samples");
	double min=p.getpar<double>(cfg,"dim."+dimid+".min");
	double max=p.getpar<double>(cfg,"dim."+dimid+".max");
	return new Dimension(samples,bind(initializer,placeholders::_1,min,max),to_linear,from_linear);
}

auto l10 = [](double x) { return (x > 0.0) ? log10(x) : -300.0;};
auto e10 = [](double x) { return exp10(x); };

void State::initializeParticle(Particle& p,boost::property_tree::ptree& cfg)
{
	using std::bind;
	p.configure(cfg.get_child("particle.default"));
	p.configure(cfg.get_child("particle."+p.id));

	// add dimension for energies
	p.ps.add(
		createDimension(
			p,
			"energy",
			initEnergyPoints,
			l10, e10,
			cfg
		)
	);
	
	// add dimension for r
	double innerRadius = exp(logr.front())*schwRadius;
	double edgeRadius = exp(logr.back())*schwRadius;
	//p.ps.add(new Dimension(nR,bind(initGridLogarithmically,placeholders::_1,
	//						innerRadius*sqrt(paso_r),edgeRadius/sqrt(paso_r)),l10,e10));
    //p.ps.add(new Dimension(nR,bind(initGridLogarithmically,placeholders::_1,
	//						innerRadius*sqrt(paso_r),edgeRadius/sqrt(paso_r)),l10,e10));
	p.ps.add(new Dimension(nR,bind(initGridLogarithmically,placeholders::_1,
							innerRadius, edgeRadius) ,l10,e10));
							
	// add dimension for rcd
	//p.ps.add(new Dimension(nRcd,bind(initGridLogarithmically,placeholders::_1,
	//						rTr*sqrt(paso_rCD),rOutCD/sqrt(paso_rCD)),l10,e10));
	p.ps.add(new Dimension(nRcd,bind(initGridLogarithmically,placeholders::_1,
							rTr, rOutCD),l10,e10));
	
	p.initialize();
}
