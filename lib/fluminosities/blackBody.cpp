#include "blackBody_RJ.h"
#include <fmath/physics.h>

double bb_RJ(double energy, double temp, double radius) {
	double frecuency = energy / planck;
	
	double flux = pi*2.0*frecuency*frecuency*boltzmann*temp / cLight2;
	return flux * 4.0*pi*radius*radius;
}

double bb_P(double energy, double temp, double radius) {
	double frecuency = energy / planck;
	
	double flux = 2.0*planck*frecuency*frecuency*frecuency / cLight2 / (1.0-exp(energy/(boltzmann*temp)));
	return flux * 4.0*pi*radius*radius;
}