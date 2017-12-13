#include "targetFields.h"


#include "functions.h"
#include "modelParameters.h"
#include <fluminosities/thermalBremss.h>
#include <fparameters/SpaceIterator.h>

//#include <fmath/physics.h>
//#include <fparameters/parameters.h>

//#include <boost/property_tree/ptree.hpp>


void tpfFill(State& st)
{	  
	st.tpf.fill([&](const SpaceIterator& i){
		 double E = i.val(DIM_E);
         double r = i.val(DIM_R);
         double theta = i.val(DIM_THETA);
		 
         const double denf_e{ st.denf_e.get(i) };
         const double denf_i{ st.denf_i.get(i) };
		
		 double temp = temp_e(r, theta);
		 double jBr = jBremss(E, temp, denf_e, denf_i);  //[jBremss] = erg cm^-3 ster^-1
		 return jBr; //4.0*pi*jBr/P2(E);
	  });
}

	  

