#include "radiativeLosses.h"

#include "write.h"
#include "messages.h"
#include "globalVariables.h"

#include <flosses/nonThermalLosses.h>
#include <flosses/lossesSyn.h>
#include <flosses/lossesIC.h>
#include <flosses/lossesHadronics.h>
#include <flosses/lossesPhotoHadronic.h>


#include <fparameters/SpaceIterator.h>
#include <fparameters/Dimension.h>
#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>

void radiativeLosses(Particle& p, State& st, const std::string& filename)
{
	show_message(msgStart, Module_radLosses);


	double rg = gravitationalConstant*blackHoleMass / cLight;

	static const double accEfficiency = GlobalConfig.get<double>("accEfficiency");


	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);

	file << "Log(E/eV)" 
		<< "\t" << "r/rg"
		<< "\t" << "Synchr"
		<< "\t" << "adv"
		<< "\t" << "Acc"
		<< "\t" << "IC / pp"
		<< "\t" << "pg"
		<< std::endl;

	double phEmin = st.photon.emin();
	double phEmax = st.photon.emax();


	p.ps.iterate([&](const SpaceIterator& i) {
		

		double E = i.val(DIM_E);
		double r = i.val(DIM_R);
		const double B = st.magf.get(i);
		
		double rB1   = r/sqrt(paso_r);
		double rB2   = r*sqrt(paso_r);
		double delta_r = rB2-rB1;
		
		const double density = st.denf_e.get(i)+st.denf_i.get(i); //ver

		double beta = sqrt(1.0 - 1.0 / P2(E/(p.mass*cLight2)));

		double v = 0.1*cLight*beta;
		
		double fmtE = log10(E / 1.6e-12);

		

		double eAdv = v/delta_r;  //ver deberia ser v_r XXX

		double eSyn = lossesSyn(i.val(DIM_E), B, st.electron) / i.val(DIM_E);

		double eAcc = accelerationRate(E, B, accEfficiency);
		
		
		file << fmtE << '\t' << r/rg
					<< "\t" << safeLog10(eSyn)
					<< "\t" << safeLog10(eAdv)
					<< "\t" << safeLog10(eAcc);
		
		
		if(p.id == "ntElectron"){
			double eIC = lossesIC(E, p, st.photon.distribution, i.coord, phEmin, phEmax)/ i.val(DIM_E);  //ver unidades del distributino XXX
			
			//lossesIC(i.val(DIM_E), st.electron,
			//[&E, &r, &st.photon](double E) {
			//return photon.distribution.interpolate({ {DIM_E, E },{ DIM_R, r } })},Emin, Emax) / i.val(DIM_E); 
			
			file << "\t" << safeLog10(eIC)
				 << std::endl;
			
		}
		else if(p.id == "ntProton"){
			double pPP = lossesHadronics(E, density, p);
			
			double pPG = lossesPhotoHadronic(E, p, st.photon.distribution, i, phEmin, phEmax);
			
			
			file << "\t" << safeLog10(pPP)
				 << "\t" << safeLog10(pPG)
				 << std::endl;
			}
			

	},{-1,-1,0});


	file.close();
}

