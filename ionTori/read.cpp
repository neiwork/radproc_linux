#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <fparameters/parameters.h>

#include "globalVariables.h"
#include "read.h"

void readThermalProcesses(int flags[])
{
	for (int i=0;i<numProcesses;i++) {
		flags[i] = GlobalConfig.get<int>("ThermalProcessNumber."+std::to_string(i));
	}
}