#include "timescales.h"
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

void timescales(Particle& p, State& st, const std::string& filename)
{
	static const double accEfficiency_PL = GlobalConfig.get<double>("nonThermal.injection.accEfficiency_PL");
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

	p.ps.iterate([&](const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		double tAdv = accretionTime(r);
		double B = st.magf.get(iR);
		double height = height_fun(r);
		double eMaxHillas = electronCharge*B*height;
		p.ps.iterate([&](const SpaceIterator& iRE) {
			double E = iRE.val(DIM_E);
			double tAcc = 1.0/accelerationRate(E,B);
			double tDiff = diffusionTimeIso(E,r,p,B,height);
			double tSyn = 1.0/lossesSyn(E,B,p)/E;

			file << "\t" << safeLog10(r/schwRadius)
				 << "\t" << safeLog10(E/1.6e-12)
				 << "\t" << safeLog10(tAcc)
				 << "\t" << safeLog10(tAdv)
				 << "\t" << safeLog10(tDiff)
				 << "\t" << safeLog10(eMaxHillas/1.6e-12)
				 << "\t" << safeLog10(tSyn);
			
			if (p.id == "ntElectron") {
				double eIC = lossesIC(E,p,st.photon.distribution,iR.coord,st.photon.emin(),
										st.photon.emax())/E;
				double eBrem = lossesBremss(E,st.denf_e.get(iR)+st.denf_i.get(iR),p)/E;  
				file << "\t" << safeLog10(eIC)
					 << "\t" << safeLog10(eBrem)
					 << std::endl;
			} else if(p.id == "ntProton") {
				double pPP = lossesHadronics(E,st.denf_i.get(iR),p)/E;
				double pPG = lossesPhotoHadronic(E,p,st.photon.distribution,iR,st.photon.emin(),
								st.photon.emax())/E;
				file << "\t" << safeLog10(pPP)
					 << "\t" << safeLog10(pPG)
					 << std::endl;
			}
		},{-1,iR.coord[DIM_R],0});
	},{0,-1,0});
	file.close();
}

void radiativeLossesNeutron(Particle& n, State& st, const std::string& filename)
{
	static const double accEfficiency = GlobalConfig.get<double>("nonThermal.injection.accEfficiency");
	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);

	file << "Log(E/GeV)" 
		<< "\t" << "r [M]"
		<< "\t" << "Escape"
		<< "\t" << "Decay"
		<< "\t" << "pp"
		<< "\t" << "pg"
		<< std::endl;
	
	double phEmin = st.photon.emin();
	double phEmax = st.photon.emax();
	n.ps.iterate([&](const SpaceIterator& i) {
		
		double E = i.val(DIM_E);
		double gamma = E/(neutronMass*cLight2);
		double r = i.val(DIM_R);
		double fmtE = log10(E/1.6e-3);
		
		double nEscape = cLight/r;
		double nDecay = 1.0/(gamma*neutronMeanLife);
		double nPP = lossesHadronics(E,st.denf_i.get(i),n)/E;
		double nPG = lossesPhotoHadronic(E,n,st.photon.distribution,i,phEmin,phEmax)/E;

		file << fmtE << "\t" << r/schwRadius
					 << "\t" << safeLog10(nEscape)
					 << "\t" << safeLog10(nDecay)
					 << "\t" << safeLog10(nPP)
				 	 << "\t" << safeLog10(nPG) << endl;
	},{-1,-1,0});
	file.close();
	
}

