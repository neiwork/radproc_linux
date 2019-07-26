#include "radiativeLosses.h"
#include "write.h"
#include "messages.h"
#include "globalVariables.h"
#include "adafFunctions.h"

#include <flosses/nonThermalLosses.h>
#include <flosses/lossesSyn.h>
#include <flosses/lossesIC.h>
#include <flosses/lossesBrem.h>
#include <flosses/lossesHadronics.h>
#include <flosses/lossesPhotoHadronic.h>
#include <fparameters/SpaceIterator.h>
#include <fparameters/Dimension.h>
#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>

void radiativeLosses(Particle& p, State& st, const std::string& filename)
{
	static const double accEfficiency = GlobalConfig.get<double>("nonThermal.injection.accEfficiency");
	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);

	file << "Log(E/eV)" 
		<< "\t" << "r [M]"
		<< "\t" << "Acc"
		<< "\t" << "tCell"
		<< "\t" << "Adv"
		<< "\t" << "Diff"
        << "\t" << "EmaxHillas"
		<< "\t" << "Sy"
		<< "\t" << "IC/pp"
		<< "\t" << "pg/Bremss"
		<< std::endl;

	double phEmin = st.photon.emin();
	double phEmax = st.photon.emax();
	p.ps.iterate([&](const SpaceIterator& i) {
		
		double E = i.val(DIM_E);
		double r = i.val(DIM_R);
		const double B = st.magf.get(i);
		double fmtE = log10(E/1.6e-12);
		double delta_r = (r/sqrt(paso_r))*(paso_r-1.0);
		double v = -radialVel(r);
		
		double eAdv = v/r;
		double eCell = v/delta_r;
		double eAcc = accelerationRate(E,B,accEfficiency);
		double eDiff = diffusionRate(E,r,B);
        double eSyn = lossesSyn(E,B,p)/E;
        
        double eMaxHillas = electronCharge*B*r*(pi-2.0*st.thetaH.get(i));

		file << fmtE << "\t" << r/schwRadius
					 << "\t" << safeLog10(eAcc)
					 << "\t" << safeLog10(eCell)
					 << "\t" << safeLog10(eAdv)
				 	 << "\t" << safeLog10(eDiff)
                     << "\t" << safeLog10(eMaxHillas/1.6e-12)
					 << "\t" << safeLog10(eSyn);
		
		if (p.id == "ntElectron") {
			double eIC = lossesIC(E,p,st.photon.distribution,i.coord,phEmin,phEmax)/E;
			double eBrem = lossesBremss(E,st.denf_e.get(i)+st.denf_i.get(i),p)/E;  
			file << "\t" << safeLog10(eIC)
				 << "\t" << safeLog10(eBrem)
				 << std::endl;
		} else if(p.id == "ntProton") {
			double pPP = lossesHadronics(E,st.denf_i.get(i),p)/E;
			double pPG = lossesPhotoHadronic(E,p,st.photon.distribution,i,phEmin,phEmax)/E;
			file << "\t" << safeLog10(pPP)
				 << "\t" << safeLog10(pPG)
				 << std::endl;
		}
	},{-1,-1,0});
	file.close();
}


