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
	Matrix scatt;
	Vector esc;
	string folder{prepareOutputfolder()};
	try {
		GlobalConfig = readConfig();
		
		prepareGlobalCfg();
		
		show_message(msgStart, Module_state);
		State model(GlobalConfig.get_child("model"));
		show_message(msgEnd, Module_state);

		show_message(msgStart, Module_torusSampling);
		icMatrix(model,scatt,esc);
		show_message(msgEnd, Module_torusSampling);
		
		thermalLuminosities(model,"lum.txt",scatt,esc);

	}
	catch (std::runtime_error& e)
	{
		std::cout << "ERROR: " << e.what() << std::endl;
	}
	return 0;
}