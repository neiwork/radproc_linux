#pragma once



#include "State.h"


/*genera una tabla con los valores de tau(Egamma,r)*/
void ggIntAbsorption(State& st, const std::string& filename);
double ggOpticalDepth(int,State&,int);

/*devuelve el valor de tau(E_ix,r_current*/
double internalAbs(int E_ix, State& st, double r_current) ;