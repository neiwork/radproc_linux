#include <stdio.h>
#include "write.h"
#include "messages.h"
#include "modelParameters.h"
#include "State.h"
#include "torusSampling.h"
#include "targetFields.h"
#include "luminosities.h"

#include "distribution.h"
/*
#include "radiativeLosses.h"


#include <fparticle/Particle.h>

#include <fparameters/Dimension.h>
#include <fparameters/SpaceIterator.h>
#include <fmath/physics.h>
*/
//#include <stdexcept>

#include <fparameters/parameters.h>
#include <inout/ioutil.h>
#include <boost/property_tree/ptree.hpp>

int main()
{
	std::string folder{ prepareOutputfolder() };

	try {        
        //lag();
		
		GlobalConfig = readConfig();
		prepareGlobalCfg();
		show_message(msgStart, Module_state);

		State model(GlobalConfig.get_child("model"));
		show_message(msgEnd, Module_state);
		
		thermalDistribution(model.proton, model);
		
		Matrix a;
		show_message(msgStart, Module_torusSampling);
		torusSampling(model, a);
		show_message(msgEnd, Module_torusSampling);

//		writeMatrix("probMatrix2", model.electron, a);
			
/*		show_message(msgStart, Module_targetField);
		tpfFill_Bremss(model);  // esto completa la psv con los fotones de Bremsstrahlung
        tpfFill_Sync(model);    // idem Sync
		show_message(msgEnd, Module_targetField);
*/ 		
		//writeAllSpaceParam(folder+"\\electronDist.txt", model.electron.distribution);
		
		//writeAllSpaceParam(folder+"\\temp.txt", model.tempElectrons);
		
        //luminosities(model, folder+"\\electronLuminosities.txt", a);
		luminosities(model,"lum.txt",a);
		
		
		//writeRandTParamSpace(getFileName(folder, "\\magf"), model.magf, 0);
        //writeRandTParamSpace(getFileName(folder, "\\denf"), model.denf_e, 0);
        //writeRandTParamSpace(getFileName(folder, "\\tempf"), model.tempElectrons, 0);
        //writeRandTParamSpace("magf.dat", model.magf, 0);
        //writeRandTParamSpace("denf.dat", model.denf_i, 0);
        
		//ParamSpaceValues psv(model.electron.ps);

		//psv.fill([&](const SpaceIterator& i){
		//	double E = i.val(DIM_E);
		//	double z = i.val(DIM_R);
		//	return frad(E, z);
		//});

		
		//injection(model.electron, model);
		
		//distribution(model.electron, model);
    
	}
	catch (std::runtime_error& e)
	{
		std::cout << "ERROR: " << e.what() << std::endl;
//		throw;
	}
	return 0;
}


