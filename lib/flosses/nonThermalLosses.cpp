#include "nonThermalLosses.h"
#include <fparameters/parameters.h>
#include <fmath/physics.h>
#include <gsl/gsl_math.h>
#include <boost/property_tree/ptree.hpp>
#include <gsl/gsl_sf_bessel.h>

double adiabaticLosses(double E, double z, double vel_lat, double gamma)  //en [erg/s]
{
	static const double openingAngle = GlobalConfig.get<double>("openingAngle");

	double jetRadius = z*openingAngle;

	return gamma*2.0*(vel_lat*E / (3.0*jetRadius));  
	//en el sist lab es sin Gamma
	//termina quedando return 2.0*cLight*E / (3.0*z);
}

double BohmDiffusionCoeff(double E, double B)
{
	double larmorR = E/(electronCharge*B);
	return 1.0/3.0 * larmorR * cLight;
}

double diffusionTimeParallel(double E, double height, double B)
{
	double zeda = GlobalConfig.get<double>("nonThermal.injection.SDA.fractionTurbulent");
	double diffCoeff = zeda * BohmDiffusionCoeff(E,B);
	return height*height/diffCoeff;
}

double diffusionTimePerpendicular(double E, double height, double B)
{
	double larmorR = E/(electronCharge*B);
	double zeda = GlobalConfig.get<double>("nonThermal.injection.SDA.fractionTurbulent");
	double q = GlobalConfig.get<double>("nonThermal.injection.SDA.powerSpectrumIndex");
	double meanFreePath = larmorR/(3.0*zeda) * pow(height/larmorR,q-1.0);
	double diffCoeff = zeda * BohmDiffusionCoeff(E,B) / (1.0 + P2(meanFreePath/larmorR));
	return height*height/diffCoeff;
}

double diffCoeff_p(double E, Particle& p, double height, double B, double rho)
{
	double zeda = GlobalConfig.get<double>("nonThermal.injection.SDA.fractionTurbulent");
	double q = GlobalConfig.get<double>("nonThermal.injection.SDA.powerSpectrumIndex");
	double g = E / (p.mass*cLight2);
	double kMin = 1.0/height;
	double vA = B / sqrt(4.0*pi*rho);
	double rL = E/(electronCharge*B);
	return zeda * gsl_pow_2(p.mass*cLight)* (cLight*kMin) *gsl_pow_2(vA/cLight) * pow(rL*kMin,q-2) * g*g;
}

double diffCoeff_g(double g, Particle& p, double height, double B, double rho)
{
	double zeda = GlobalConfig.get<double>("nonThermal.injection.SDA.fractionTurbulent");
	double q = GlobalConfig.get<double>("nonThermal.injection.SDA.powerSpectrumIndex");
	double kMin = 1.0/height;
	double vA = B / sqrt(4.0*pi*rho);
	double rL = g*p.mass*cLight2/(electronCharge*B);
	return zeda * (cLight*kMin) *gsl_pow_2(vA/cLight) * pow(rL*kMin,q-2) * g*g;
}

double diffCoeff_r(double g, Particle& p, double height, double B)
{
	double zeda = GlobalConfig.get<double>("nonThermal.injection.SDA.fractionTurbulent");
	double q = GlobalConfig.get<double>("nonThermal.injection.SDA.powerSpectrumIndex");
	double kMin = 1.0/height;
	double rL = g*p.mass*cLight2/(electronCharge*B);
	return 1.0/9.0 * cLight / zeda * rL * pow(kMin*rL,1-q);
}

double diffLength(double g, Particle& p, double r, double height, double B, double vR)
{
	return diffCoeff_r(g,p,height,B)/abs(vR);
}


double diffusionTimeTurbulence(double E, double height, Particle& p, double B)   //en [s]
{
	double zeda = GlobalConfig.get<double>("nonThermal.injection.SDA.fractionTurbulent");
	double q = GlobalConfig.get<double>("nonThermal.injection.SDA.powerSpectrumIndex");
	double larmorRadius = E/(electronCharge*B);
	double tEscape = height/cLight;
	return 9.0*zeda*tEscape*pow(larmorRadius/height, q-2.0);
}

double accelerationTimeSDA(double E, Particle& p, double B, double height, double rho)
{
	double zeda = GlobalConfig.get<double>("nonThermal.injection.SDA.fractionTurbulent");
	double q = GlobalConfig.get<double>("nonThermal.injection.SDA.powerSpectrumIndex");
	double larmorRadius = E/(electronCharge*B);
	double tEscape = height/cLight;
	double alfvenVel = B/sqrt(4.0*pi*rho);
	return (1.0/zeda) / P2(alfvenVel/cLight) * tEscape * pow(larmorRadius/height,2.0-q);
}

double accelerationRateSDA(double E, Particle& p, double B, double height, double rho)
{
	double gamma = E / (p.mass*cLight2);
	return diffCoeff_g(gamma,p,height,B,rho) / (gamma*gamma);
}


double relaxTime_e(double E, double temp, double dens) {
	
	double gamma = E / electronRestEnergy;
	double lnLambda = 20.0;
	double theta = boltzmann*temp / electronRestEnergy;
	
	if (gamma > 2.0 && theta > 0.3) {
		double k1 = gsl_sf_bessel_K1(1.0/theta);
		double k2 = gsl_sf_bessel_Kn(2, 1.0/theta);
		double factor = abs(k1/k2 - 1.0/gamma);
		return 2.0/3.0 * gamma / (dens*thomson*cLight*(lnLambda+9.0/16.0-0.5*log(2.0))) /
				factor;
	} else
		return 4.0*sqrt(pi)*pow(theta, 1.5) / (dens*thomson*cLight*lnLambda);
}

double relaxTime_p(double E, double temp, double dens) {
	
	double gamma = E / (protonMass*cLight2);
	double lnLambda = 20.0;
	double theta = boltzmann*temp / (protonMass*cLight2);
	double b1 = 70e6*EV_TO_ERG / (protonMass*cLight2);
	double b2 = 500e6*EV_TO_ERG / (protonMass*cLight2);
	if ( gamma-1.0 > b1 && gamma-1.0 < b2 ) {
		double beta = sqrt(1.0-1.0/(gamma*gamma));
		double sigma_h = 2.3e-26;
		return 4.0 * gamma*gamma * beta / (gamma*gamma - 1.0) / (dens*sigma_h*cLight);
	} else
		return 4.0*sqrt(pi)*pow(theta, 1.5) * P2(protonMass/electronMass) / (dens*thomson*cLight*lnLambda);
}

double accelerationRate(double E, double B) //en [s]^-1
{
	double accEff = GlobalConfig.get<double>("nonThermal.injection.PL.accEfficiency");
	return accEff*cLight*electronCharge*B/E;
}


double escapeRate(double size, double vel) //en [1/s]
{
	//static const double Gamma = GlobalConfig.get<double>("Gamma");

	return vel /size;  //ver si necesito algun gamma para escribirlo en el FF
}