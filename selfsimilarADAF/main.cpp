#include <ssADAF/ssADAF.h>

#include <stdio.h>
#include "modelParameters.h"
#include "messages.h"
#include "State.h"
#include "write.h"
//#include "radiativeLosses.h"
#include "torusSampling2.h"

//#include "distribution.h"
//#include "luminosities.h"
//#include "functions.h"

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
		show_message(msgStart, Module_state);
		State model(GlobalConfig.get_child("model"));
		show_message(msgEnd, Module_state);
		show_message(msgStart, Module_targetField);
		
		torusSampling2();
        
	}
	
	catch (std::runtime_error& e)
	{
		std::cout << "ERROR: " << e.what() << std::endl;
	}
	return 0;	
}