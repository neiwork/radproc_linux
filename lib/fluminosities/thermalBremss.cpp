#include "thermalBremss.h"

#include <fparameters/parameters.h>
#include <fmath/physics.h>

double eiAuxFunction(double norm_temp) {
    if (norm_temp < 1.0)
        return 4.0 * sqrt(2.0 * norm_temp / P3(pi)) * (1.0 + 1.781 * pow(norm_temp, 1.34));
    else if (norm_temp >= 1.0)
        return 9.0 * norm_temp / (2.0 * pi) * (log(1.123 * norm_temp + 0.48) + 1.5);
}

double eiCoolingRate(double norm_temp, double dens_e, double dens_i) {
    
    return dens_i * dens_e * thomson * fineStructConst * electronMass * P3(cLight) * eiAuxFunction(norm_temp);
}    // pasar solo la temperatura en un punto.

double eeAuxFunction(double norm_temp) {
    
    if (norm_temp < 1.0)
        return 20.0/(9.0 * sqrt(pi)) * (44.0-3.0*P2(pi)) * pow(norm_temp, 1.5) * 
        (1.0 + 1.1*norm_temp + P2(norm_temp) - 1.25 * pow(norm_temp, 2.5));
    else if (norm_temp >= 1.0)
        return 24.0 * norm_temp * (log(2.0 * 0.5616 * norm_temp) + 1.28);
}

double eeCoolingRate(double norm_temp, double dens_e) {
    return dens_e * dens_e * electronRadius * electronRadius *
    fineStructConst * electronMass * P3(cLight) * eeAuxFunction(norm_temp);
}

double gauntFactor(double temp_aux) {
    if (temp_aux < 1.0)
        return sqrt(3.0 / pi * temp_aux );
    else if (temp_aux >= 1.0)
        return sqrt(3.0) / pi * log(4.0 / 0.576965 * temp_aux);
}

double jBremss(double energy, double temp, double dens_i, double dens_e) {
    
    double norm_temp = boltzmann * temp / (electronMass * cLight2);
	
    double aux = eeCoolingRate(norm_temp, dens_e) + eiCoolingRate(norm_temp, dens_e, dens_i);
    double temp_aux = boltzmann * temp / energy;
    
    return aux/(4.0*pi) * planck/(boltzmann*temp) * exp(- (1.0/temp_aux)) * gauntFactor(temp_aux);
}