#include "blob.h"
#include "adafFunctions.h"
#include "globalVariables.h"
#include "State.h"
#include <fparameters/parameters.h>
#include <boost/property_tree/ptree.hpp>

void blob(State& st)
{
	factorDensity = GlobalConfig.get<double>("factorDensity");
	tAccBlob = GlobalConfig.get<double>("tAccBlob");    // [lightcurve variability time in sec]
	double sumt = 0.0;
	double pasoRaux = pow(100.0,1.0/1000);
	rBlob = schwRadius*1.1;
	while ( sumt < tAccBlob) {
		double dr = rBlob*(pasoRaux-1.0);
		sumt += dr / (-radialVel(rBlob));
		rBlob *= pasoRaux;
	}
}