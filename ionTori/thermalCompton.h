#pragma once

#include <fparameters/SpaceIterator.h>
#include <fparameters/Dimension.h>
//#include <fparameters/parameters.h>

double comptBremss(double E, double theta_e, double r, double theta, const SpaceCoord& distCoord, 
						ParamSpaceValues& denf_e,const ParamSpaceValues& jBr);