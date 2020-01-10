#pragma once



#include "State.h"

//double emiToLumi(const ParamSpace& pps, ParamSpaceValues& psv, double E, int t_ix);

void nonThermalRadiation(State& st, const std::string& filename);

//double Llab(double Lint);