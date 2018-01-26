#include "blackBody.h"
#include <fmath/physics.h>

double bb_RJ(double frecuency, double temp) {
	return pi*2.0*frecuency*frecuency*boltzmann*temp / cLight2;
}

double bb(double frecuency, double temp) {	
	return pi*2.0*planck*frecuency*frecuency*frecuency / cLight2 / (1.0-exp(planck*frecuency/(boltzmann*temp)));
}