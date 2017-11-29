#include <thermalSync.h>

#include <fparameters/parameters.h>
#include <fmath/physics.h>

double mAux(double frecuency, double r, double theta, const ParamSpaceValues& magf, const Particle& c) {
    
    norm_temp = boltzmann * temp_e(r, theta) / (c.mass * cLight2);
    nu0 = (electronCharge * magf) / (2.0 * PI * c.mass * cLight);
    xM = (2.0 * frecuency) / (3.0 * nu0 * (P2(norm_temp)));
    
    alpha = 1.0;
    beta   = 1.0;
    gamma = 1.0;
    
    return (4.0505 * alpha / pow(xM, 1.0/6.0) ) * (1.0 + 0.4*beta / pow(xM, 1.0/4.0) + 
                0.5316 * gamma / sqrt(xM) ) * exp(-1.8899 * pow(xM, 1.0/3.0));
}

double jSync(double frecuency, const Particle& c, const SpaceCoord& psc, const ParamSpaceValues& magf) {
    r =
    theta = 
    norm_temp = boltzmann * temp_e(r, theta) / (c.mass * cLight2);
    
    return (1.0/4.0*PI) * P2(electronCharge)/(sqrt(3.0)*cLight) * (4.0 * PI * denf_e * frecuency) / 
    bessel2(1.0/norm_temp) * mAux(frecuency, r, theta, magf);    // falta definir esta bessel
}   // esto debería tener unidades de erg cm^-3 ster^-1
