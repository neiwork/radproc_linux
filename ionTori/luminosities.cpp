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
#include <fmath/mathFunctions.h>
#include <fmath/physics.h>



void luminosities(State& st, const std::string& filename) {
    
    std::ofstream file;
	file.open(filename.c_str(), std::ios::out);

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

        const double magf{ st.magf.get(i) };
        const double denf_e{ st.denf_e.get(i) };
        const double denf_i{ st.denf_i.get(i) };
        double fmtE = log10(i.val(DIM_E) / 1.6e-12);

        double energy = i.val(DIM_E);

        double r = i.val(DIM_R);
        double theta = i.val(DIM_THETA);
		
		double temp = temp_e(r, theta);
		double theta_e = temp/(electronMass*cLight2);
		
        
        double jSy = jSync(energy, temp, magf, denf_e);
        double jBr = st.tpf.get(i);  // jBremss(energy, temp, denf_e, denf_i);  
		double ic = comptBremss(energy, theta_e, r, theta, i, 
						st.denf_e, st.tpf);
		
		double ic_2 = luminosityIC(energy, st.electron, i, st.tpf, temp/1.0e3);
    
        file << fmtE << "\t" << r
                            << "\t" << theta
                            << "\t" << safeLog10(jSy)
                            << "\t" << safeLog10(jBr)
							<< "\t" << safeLog10(ic)
							<< "\t" << safeLog10(ic_2)
                            << std::endl;
                            
    }, { -1, -1, -1 });
    
    file.close();
    
}