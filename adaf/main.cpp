#include <stdio.h>
#include "messages.h"
#include "modelParameters.h"
#include "State.h"
#include "icMatrix.h"
#include "thermalLuminosities.h"
#include "globalVariables.h"
#include <fparameters/parameters.h>
#include <inout/ioutil.h>
#include <iostream>
#include <boost/property_tree/ptree.hpp>

using namespace std;
int main()
{
	Matrix scattADAF,scattCD,absCD;
	Vector esc,escCD;
	string folder{prepareOutputfolder()};
	try {
		GlobalConfig = readConfig();
		
		prepareGlobalCfg();
		
		show_message(msgStart, Module_state);
		State model(GlobalConfig.get_child("model"));
		show_message(msgEnd, Module_state);

		show_message(msgStart, Module_torusSampling);
		icMatrix(model,scattADAF,scattCD,absCD,esc,escCD);
		show_message(msgEnd, Module_torusSampling);
		thermalLuminosities(model,"lum.txt",scattADAF,scattCD,absCD,esc,escCD);

	}
	catch (std::runtime_error& e)
	{
		std::cout << "ERROR: " << e.what() << std::endl;
	}
	return 0;
}