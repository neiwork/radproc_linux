#include "luminosities.h"

#include "modelParameters.h"
#include "write.h"

#include "functions.h"

#include <fparameters/SpaceIterator.h>

#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>


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
          << std::endl;
          
    st.electron.ps.iterate([&](const SpaceIterator& i){

        const double magf{ st.magf.get(i) };
        const double denf_e{ st.denf_e.get(i) };
        const double denf_i{ st.denf_i.get(i) };
        double fmtE = log10(i.val(DIM_E) / 1.6e-12);

        double energy = i.val(DIM_E);

        double r = i.val(DIM_R);
        double theta = i.val(DIM_THETA);
		
		double temp = temp_e(r, theta);
		
        
        double jSy = jSync(energy, temp, magf, denf_e);
        double jBr = jBremss(energy, temp, denf_e, denf_i);  
    
        file << fmtE << "\t" << r
                            << "\t" << theta
                            << "\t" << safeLog10(jSy)
                            << "\t" << safeLog10(jBr)
                            << std::endl;
                            
    }, { -1, -1, -1 });
    
    file.close();
    
}