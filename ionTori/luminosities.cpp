#include "luminosities.h"

#include "modelParameters.h"
#include "write.h"

#include "functions.h"
#include "maxwellJuttner.h"
#include "thermalCompton.h"

#include <fparameters/SpaceIterator.h>

#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>


#include <fluminosities/luminosityIC.h>
#include <fluminosities/thermalSync.h>
#include <fluminosities/thermalBremss.h>
#include <fluminosities/thermalIC.h>
#include <fmath/mathFunctions.h>
#include <fmath/physics.h>



void luminosities(State& st, const std::string& filename) {
    
    std::ofstream file;
	file.open(filename.c_str(), std::ios::out);
	
	double rmax = 3.0 * GlobalConfig.get<double>("rCenter");
    double rmin = rmax / 3.0;
    rmax = 2.1 * rmin;
	int nR = st.electron.ps[DIM_R].size() - 1; // GlobalConfig.get<double>("model.particle.default.dim.radius.samples");
	double dr = (rmax-rmin)/nR;
	
	double thetamin =0.0;                           
    double thetamax = pi/5.5;
	int nT = st.electron.ps[DIM_THETA].size() - 1; 
	double dt = (thetamax-thetamin)/nT;

	file << "Log(E/eV)" 
      //    << "\t" << "r [M]"
      //   << "\t" << "Theta"
		  << "\t" << "Synchr"
          << "\t" << "Brems"
		  << "\t" << "IC"
		  << "\t" << "IC_prueba"
          << std::endl;
          
    //st.electron.ps.iterate([&](const SpaceIterator& i){
	//el iterate hay que hacerlo sobre el photon, asi las energias que recorremos son las de los fotones
	st.photon.ps.iterate([&](const SpaceIterator& i){   

        // const double denf_e{ st.denf_e.get(i) };
        // const double denf_i{ st.denf_i.get(i) };
        double fmtE = log10(i.val(DIM_E) / 1.6e-12);

        double energy = i.val(DIM_E);

        double r = i.val(DIM_R);
        double theta = i.val(DIM_THETA);
        
        double norm_temp = boltzmann * st.tempElectrons.get(i) / electronMass / cLight2;
		
        double jBr = st.tpf1.get(i)*(energy*energy);    // jBremss(energy, temp, denf_e, denf_i);
        double jSy = st.tpf2.get(i);                                 //jSync(energy, temp, magf, denf_e);
        
		double ic = comptBremss(energy, norm_temp, r, theta, i, st.denf_e, st.tpf1);
                        
        double jIC1 = jIC_Bremss(energy, norm_temp, r, theta, i, st.denf_e, jBr, rmax);
        double jIC2 = jIC_Sync(energy, norm_temp, r, theta, i, st.denf_e, jSy, rmax);
        
        double jIC = jIC1 + jIC2;
		
		double ic_2 = luminosityIC(energy, st.electron, i, st.tpf1, temp_e(r, theta)/1.0e3);
		
		
		double area_i = 2.0*pi*(r*dt)*dr;
		double vol_i = 2.0*pi*r*(r*dt)*dr;
		
		double emi_to_lumi = cLight*area_i;  //from erg/cm^3 to erg/s
		
    
        file << fmtE << "\t" << r
                            << "\t" << theta
                            << "\t" << safeLog10(jSy*emi_to_lumi)
                            << "\t" << safeLog10(jBr*emi_to_lumi)
							<< "\t" << safeLog10(ic*emi_to_lumi)
							<< "\t" << safeLog10(ic_2*vol_i) 
                            << "\t" << safeLog10(jIC * emi_to_lumi)
                            << std::endl;
                            
    }, { -1, -1, -1 });
    
    file.close();
    
}