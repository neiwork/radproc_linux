#include "luminosityPhotoHadronic.h"

#include <fparameters/parameters.h>
#include <fmath/RungeKutta.h>
//#include <fmath/interpolation.h>
#include <flosses/crossSectionInel.h>
#include <finjection/pgammaPionInj.h>
#include <fmath/physics.h> 

//double lossesPhotoHadronic(double E, Particle& particle, const ParamSpaceValues& tpf, const SpaceCoord& psc, double phEmin, double phEmax)
//{  //E=Ep

double luminosityPhotoHadronic(double E, Particle& p, const ParamSpaceValues& tpf, const SpaceCoord& psc, double phEmin, double phEmax) //fun1 tpf)
{

	double diezE = 10.0*E;
	double distCreator = p.distribution.interpolate({ { 0, diezE } }, &psc);//(diezE);// interpol(diezE, creator.energyPoints, Ncreator, Ncreator.size() - 1);

	double t_1   = t_pion_PH(diezE, p, 
						[&tpf, &psc](double x) {return tpf.interpolate({ {0, x } }, &psc); },
						phEmin, phEmax);     //esto no es lossesPH porque son perdidas solo del canal de produccion de piones
			
	double omega = omegaPH(diezE, p, //tpf, tpEmin, tpEmax);
						[&tpf, &psc](double x) {return tpf.interpolate({ {0, x } }, &psc); },
						phEmin, phEmax);

	if (omega > 0.0)	{
		double averageInel = t_1/omega;

		double k1 = 0.2;

		double k2 = 0.6;

		double p1 = (k2-averageInel)/(k2-k1);
	
		double nNeutralPion = 1.0-0.5*p1;

		double luminosity = 20.0*nNeutralPion*omega*distCreator;

		return luminosity*P2(E);
	}
	else	{return 0.0;}

	
}
