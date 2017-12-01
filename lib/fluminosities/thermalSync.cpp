#include <thermalSync.h>

#include <fparameters/parameters.h>
#include <fmath/physics.h>
#include <fmath/mathFunctions.h>

double mAux(double frecuency, double norm_temp, double magfield) {
    
    double nu0 = (electronCharge * magfield) / (2.0 * pi * electronMass * cLight);
    double xM = (2.0 * frecuency) / (3.0 * nu0 * (P2(norm_temp)));
    
    double alpha = 1.0;
    double beta   = 1.0;
    double gamma = 1.0;
    
    return (4.0505 * alpha / pow(xM, 1.0/6.0) ) * (1.0 + 0.4*beta / pow(xM, 1.0/4.0) + 
                0.5316 * gamma / sqrt(xM) ) * exp(-1.8899 * pow(xM, 1.0/3.0));
}

double jSync(double energy, double temp, double magfield, double denf_e) {
    
    double frecuency = energy / planck;
    double norm_temp = boltzmann * temp / (electronMass * cLight2);
    
    return (1.0/4.0*pi) * P2(electronCharge)/(sqrt(3.0)*cLight) * (4.0 * pi * denf_e * frecuency) / 
    bessk(2, 1.0/norm_temp) * mAux(frecuency, norm_temp, magfield);
}   // esto debería tener unidades de erg cm^-3 ster^-1
