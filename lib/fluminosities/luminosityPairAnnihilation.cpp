#include "luminosityPairAnnihilation.h"

#include <fparameters/parameters.h>
#include <fmath/RungeKutta.h>
#include <fmath/physics.h>
#include <gsl/gsl_math.h>

double rateAnnihilation(double Eg, double Ep, double Em)
{
	double gamma_ph = Eg / electronRestEnergy;
	double gamma_p = Ep / electronRestEnergy;
	double gamma_m = Em / electronRestEnergy;
	double gUcm = 0.0;
	if ( (gamma_ph > gamma_p && gamma_ph < gamma_m) || (gamma_ph < gamma_p && gamma_ph > gamma_m) )
		gUcm = gamma_ph * (gamma_p + gamma_m - gamma_ph);
	else
		gUcm = sqrt(0.5*(1.0+gamma_p*gamma_m+sqrt((gsl_pow_2(gamma_p)-1.0)*(gsl_pow_2(gamma_m)-1.0))));
	
	double constant = pi * gsl_pow_2(electronRadius) * cLight / gsl_pow_2(gamma_p*gamma_m);
	return constant * ( gUcm*gUcm / (abs(gamma_ph-gamma_p)+2.0/pi) + gUcm*gUcm / (abs(gamma_ph-gamma_m)+2.0/pi) );
}

double luminosityPairAnnihilation(double Eg, Particle& p, const SpaceCoord& psc)
{
	double gamma_ph = Eg / electronRestEnergy;
	double gammaA = (4*gsl_pow_2(gamma_ph)+1.0)/(4*gamma_ph);
	double E1min = (gamma_ph > 0.5 ? p.emin() : max(gammaA*electronRestEnergy,p.emin()));
	double E1max = p.emax();
	return integSimpsonLog(E1min,E1max,[&p,Eg,gamma_ph,&psc](double Ep)
			{
				double gamma_p = Ep / electronRestEnergy;
				double Np = (Ep > p.emin() && Ep < p.emax()) ? 
						p.distribution.interpolate({{0,Ep}},&psc) / 2.0 : 0.0;
				double gammaB = (2*gsl_pow_2(gamma_ph)-2*gamma_ph+1) / (2*gamma_ph-1.0);
				double E2min = p.emin();
				if (gamma_p > gammaB && gamma_ph >= 0.5)
					E2min = p.emin();
				else if (gamma_p < gammaB && gamma_ph >= 1.0) {
					double beta_p = sqrt(1.0 - 1.0/(gamma_p*gamma_p));
					double Fp = 2.0*gamma_ph - gamma_p * (1.0 + beta_p);
					double gamma_aux_p = 0.5 * (Fp + 1.0/Fp);
					E2min = gamma_aux_p * electronRestEnergy;
				} else {
					double beta_p = sqrt(1.0 - 1.0/(gamma_p*gamma_p));
					double Fm = 2.0*gamma_ph - gamma_p * (1.0 - beta_p);
					double gamma_aux_m = 0.5 * (Fm + 1.0/Fm);
					E2min = gamma_aux_m * electronRestEnergy;
				}
				double E2max = p.emax();
				double integral = integSimpsonLog(E2min,E2max,[&p,Eg,Ep,&psc](double Em)
						{
							double Nm = (Em > p.emin() && Em < p.emax()) ? 
									p.distribution.interpolate({{0,Em}},&psc) / 2.0 : 0.0;
							double Rate = rateAnnihilation(Eg,Ep,Em);
							return Rate*Nm;
						},50);
				return Np*integral;
			},50)*Eg;
}

