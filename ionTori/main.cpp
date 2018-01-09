#include <stdio.h>


#include "modelParameters.h"
#include "targetFields.h"
#include "State.h"
#include "write.h"
#include "radiativeLosses.h"
#include "distribution.h"
//#include "luminosities.h"
#include "luminosities2.h"

#include <inout/ioutil.h>


#include <fparticle/Particle.h>
#include <fparameters/parameters.h>
#include <fparameters/Dimension.h>
#include <fparameters/SpaceIterator.h>
#include <fmath/physics.h>

#include <boost/property_tree/ptree.hpp>
#include <stdexcept>

int main()
{
	std::string folder{ prepareOutputfolder() };

	try {
		GlobalConfig = readConfig();
		prepareGlobalCfg();
		State model(GlobalConfig.get_child("model"));
		
		tpfFill_Bremss(model);  // esto completa la psv con los fotones de Bremsstrahlung
        tpfFill_Sync(model);      // idem Sync
		
		// radiativeLosses(model, folder+"\\electronLosses.txt");
		
		thermalDistribution(model.electron, model);
		writeAllSpaceParam(folder+"\\electronDist.txt", model.electron.distribution);
		
		//writeAllSpaceParam(folder+"\\bremss.txt", model.tpf1);
		
        luminosities2(model, folder+"\\electronLuminosities.txt");
		
		
		//writeRandTParamSpace(getFileName(folder, "\\magf"), model.magf, 0);
       // writeRandTParamSpace("magf.dat", model.magf, 0);
        //writeRandTParamSpace("denf.dat", model.denf_i, 0);
        
		//ParamSpaceValues psv(model.electron.ps);

		//psv.fill([&](const SpaceIterator& i){
		//	double E = i.val(DIM_E);
		//	double z = i.val(DIM_R);
		//	return frad(E, z);
		//});

		
		//injection(model.electron, model);
		
		//distribution(model.electron, model);

		

		//processes(model, getFileName(folder, "luminosity"));

	}
	catch (std::runtime_error& e)
	{
		std::cout << "ERROR: " << e.what() << std::endl;
//		throw;
	}

	return 0;
}




/*int main(int argc, char **argv)
{
	printf("hello world\n");
	return 0;
}*/
