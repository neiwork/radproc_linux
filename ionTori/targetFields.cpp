#include "targetFields.h"


#include "functions.h"
#include "modelParameters.h"
#include <fluminosities/thermalBremss.h>
#include <fluminosities/thermalSync.h>
#include <fparameters/SpaceIterator.h>

//#include <fmath/physics.h>
//#include <fparameters/parameters.h>

//#include <boost/property_tree/ptree.hpp>


void tpfFill_Bremss(State& st) {	  
	st.tpf1.fill([&](const SpaceIterator& i) {
		 double energy = i.val(DIM_E);
		 
         const double denf_e{ st.denf_e.get(i) };
         const double denf_i{ st.denf_i.get(i) };
         const double temp{ st.tempElectrons.get(i) };
        
		 double jBr = jBremss(energy, temp, denf_e, denf_i);  //[jBremss] = erg cm^-3 ster^-1 s^-1 Hz^-1
		 return jBr*4.0*pi/(energy*energy);
    });
}

void tpfFill_Sync(State& st) {
    st.tpf2.fill([&](const SpaceIterator& i) {
        double energy = i.val(DIM_E);
        
        const double magfield{ st.magf.get(i) };
        const double denf_e{ st.denf_e.get(i) };
        const double temp{ st.tempElectrons.get(i) };
        
        return jSync(energy, temp, magfield, denf_e);
    });
}

	  

