#include "neutronPgamma.h"


#include "pgammaPionInj.h"

#include <fparameters/parameters.h>

//con incluir #include "pgammaPionInj.h" ya incluyo las siguientes funciones

//double t_pion_PH(double E, Particle& particle, fun1 tpf);   
//double omegaPH(double E, Particle& particle, fun1 tpf);



double psiEsc(double E, double tcross, Particle& neutron, const ParamSpaceValues& tpf, const SpaceCoord& distCoord, double tpEmin, double tpEmax)  //E=Ep
{

	double seda_pn = 0.5;  //probabilidad de que el n se convierta en proton psi_np = 1-psi_nn

	double tauDec = neutronMeanLife*E/(neutronMass*cLight2);

	double L = 0.0;

	int nT = 10;
	double dt = tcross/nT;
	double t = 0.0;

	for (size_t i=0; i < nT; ++i){

		L += tcross*(seda_pn*omegaPH(E,neutron,tpf,distCoord,tpEmin,tpEmax)+1.0/tauDec)*dt;

		t += dt;

	}

	double result = exp(-L);

	return result; //psiEsc*psi_pn*omegaPH(E, particle, tpf, distCoord, tpEmin, tpEmax);
}


double neutronPgamma(double E, double tcross, Particle& neutron, Particle& proton, const ParamSpaceValues& tpf, const SpaceCoord& distCoord, double tpEmin, double tpEmax)  
{
	//double protonDist = proton.dist(E);// interpol(E, proton.energyPoints, Nproton, Nproton.size() - 1);

	double t_1   = t_pion_PH(E, proton, tpf, distCoord, tpEmin, tpEmax);     //esto no es lossesPH porque son perdidas solo del canal de produccion de piones
	double omega = omegaPH(E, proton, tpf, distCoord, tpEmin, tpEmax);
	

	double emissivity;

	if (t_1 > 0.0 && omega > 0.0)	{
		double averageInel = t_1/omega;

		double seda_pn = 0.5;

		double Ep = E/(1.0-averageInel);

		double protonDist;
		if (Ep < proton.emin() || Ep > proton.emax()){protonDist = 0.0;}
		else{ protonDist = proton.distribution.interpolate({ { 0, Ep } }, &distCoord);}

//		double k2 = 0.6;
//		double p1 = (k2-averageInel)/(k2-k1);	
//		double nChargedPion = 2.0-1.5*p1;

		double pEsc = psiEsc(Ep, tcross, neutron, tpf, distCoord, tpEmin, tpEmax);
		double tauEsc = 1.0/(pEsc*seda_pn*omega);

		emissivity = protonDist/(1.0-averageInel)/tauEsc;
	}
	else	{ emissivity = 0; }
	
	return emissivity;
}