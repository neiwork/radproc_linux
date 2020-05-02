#include <iostream>
#include <fstream>
#include <boost/property_tree/ptree.hpp>
#include <fparameters/parameters.h>
#include <fparameters/ParamSpaceValues.h>
#include "globalVariables.h"
#include "read.h"
#include "adafFunctions.h"
#include "write.h"
#include <fparameters/SpaceIterator.h>
#include <fparameters/Dimension.h>
#include <fparameters/parameters.h>
#include <boost/property_tree/ptree.hpp>

void readThermalProcesses(int flags[])
{
	for (int i=0;i<numProcesses;i++) {
		flags[i] = GlobalConfig.get<int>("thermal.processNumber."+std::to_string(i));
	}
}

void readEandRParamSpace(const std::string& filename, ParamSpaceValues& data, int t, int vol)
{
	std::ifstream file;
	file.open(dataName(filename).c_str(), std::ios::in);
	
	data.ps.iterate([&](const SpaceIterator& i){
		double voll = (vol == 1) ? volume(i.val(DIM_R)) : 1.0;
		double aux1,aux2,aux3,dist;
		file >> aux1 >> aux2 >> aux3 >> dist;
		data.set(i,pow(10,dist)/voll);
	},{-1,-1,t});  
	file.close();
}