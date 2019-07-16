#include "neutronPgamma.h"


#include "pgammaPionInj.h"
#include <fmath/RungeKutta.h>
#include <fparameters/parameters.h>
#include <flosses/crossSectionInel.h>
#include <flosses/lossesPhotoHadronic.h>
//#include <fmath\interpolation.h>


//con incluir #include "pgammaPionInj.h" ya incluyo las siguientes funciones

//double t_pion_PH(double E, Particle& particle, fun1 tpf);   

//double omegaPH(double E, Particle& particle, fun1 tpf);


/*

double tauEsc(double E, Particle& particle, const ParamSpaceValues& tpf, const SpaceCoord& distCoord, double tpEmin, double tpEmax)  //E=Ep
{

	double mass = particle.mass;

	double psi_pn = 0.5;  //probabilidad de que el n se convierta en proton psi_np = 1-psi_nn

	double tauDec = neutronMeanLife*E/(neutronMass*cLight2);

	double tcross = 0.0;//radius / cLight; //VER

	double psiEsc = exp(-tcross*(psi_pn*omegaPH(E,particle,tpf,distCoord, tpEmin, tpEmax)+1.0/tauDec));

	return psiEsc*psi_pn*omegaPH(E,particle,tpf,distCoord, tpEmin, tpEmax);
}*/


double neutronPgamma(double E, Particle& particle, Particle& proton, const ParamSpaceValues& tpf, const SpaceCoord& distCoord, double phEmin, double phEmax)
{
	
	//double protonDist = proton.dist(E);// interpol(E, proton.energyPoints, Nproton, Nproton.size() - 1);
	double protonDist;
	if (E < proton.emin() || E> proton.emax()){
		protonDist = 0.0;
	}
	else{
		protonDist = proton.distribution.interpolate({ { 0, E } }, &distCoord); 
	}
	
	double t_1   = t_pion_PH(E, proton, tpf, distCoord, phEmin, phEmax); //esto no es lossesPH porque son perdidas solo del canal de produccion de piones
	double omega = omegaPH(E, proton, tpf, distCoord, phEmin, phEmax);
	
	double emissivity;

	if (t_1 > 0.0 && omega > 0.0)	{
		double averageInel = t_1/omega;

		double seda_pn = 0.5;
//		double k2 = 0.6;
//		double p1 = (k2-averageInel)/(k2-k1);	
//		double nChargedPion = 2.0-1.5*p1;

		emissivity = seda_pn*protonDist*omega/(1-averageInel);
	}
	else	{
		emissivity = 0;
	}
	
	return emissivity;
}
