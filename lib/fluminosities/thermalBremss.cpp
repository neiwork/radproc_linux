#include <thermalBremss.h>

#include <fparameters/parameters.h>
#include <fmath/physics.h>

double eiCoolingRate() {
    
    norm_temp = boltzmann * temp_e(r, theta) / (electronMass * cLight2);
    return denf_i * denf_e * thomson * fineStructConst * electronMass * P3(cLight) * eiAuxFunction(norm_temp);
}

double eiAuxFunction(double temp) {
    
    if (temp < 1.0)
        return 4.0 * sqrt(2.0 * temp / P3(PI)) * (1.0 + 1.781 * pow(temp, 1.34));
    else if (temp >= 1.0)
        return 9.0 * temp / (2.0 * PI) * (log(1.123 * temp + 0.48) + 1.5);
}

double eeCoolingRate() {
    
    norm_temp = boltzmann * temp_e(r, theta) / (electronMass * cLight2);
    return P2(denf_e *electronRadius) * fineStructConst * electronMass * P3(cLight) * eeAuxFunction(norm_temp);
}

double eeAuxFunction(double temp) {
    
    if (temp < 1.0)
        return 20.0/(9.0 * sqrt(PI)) * (44.0-3.0*P2(PI)) * pow(temp, 1.5) * (1.0 + 1.1*temp + P2(temp) - 1.25*pow(temp,2.5));
    else if (temp >= 1.0)
        return 24.0 * temp * (log(2.0 * 0.5616 * temp) + 1.28);
}

double jBremss(double frecuency, double r, double theta) {
    
    double aux = eeCoolingRate() + eiCoolingRate();
    double aux2 = boltzmann * temp_e(r, theta) / (planck * frecuency);
    
    return aux/(4.0*PI) * planck/(boltzmann*temp_e(r, theta)) * exp(- (1.0/aux2)) * gauntFactor(aux2, frecuency);
}

double gauntFactor(double temp, double frecuency) {
    if (temp < 1.0)
        return sqrt(3.0 / PI * temp );
    else if (temp >= 1.0)
        return sqrt(3.0) / PI * log(4.0 / 0.576965 * temp);
}