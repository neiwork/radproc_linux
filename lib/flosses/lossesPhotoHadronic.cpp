#include "lossesPhotoHadronic.h"

#include "crossSectionInel.h"
#include <fmath/RungeKutta.h>
#include <fparameters/parameters.h>
#include <fmath/physics.h>
//#include <algorithm>


double cPionPH(double u)   //limite inferior
{
	return pionThresholdPH; //145 MeV  threshold energy for pion production
}

double cPairPH(double u)   //limite inferior
{
	return pairThresholdPH; 
}

double dPH(double u, double E, double mass)   //limite superior     
{

	return 2*u*E/(mass*cLight2);   //d1 = 2*u*ep/(masa*cluz**2)
}

//double fPHPion(double u,double t, double E, double mass, fun1 tpf)   //funcion a integrar
double fPHPion(double u,double t, double E, double mass, const ParamSpaceValues& tpf, const SpaceCoord& psc)   //funcion a integrar
{
	double Nph = tpf.interpolate({ { 0, u } }, &psc	); 
	//return tpf(u)*crossSectionPHPion(t)*inelasticityPHPion(t)*t/P2(u);
	return Nph*crossSectionPHPion(t)*inelasticityPHPion(t)*t/P2(u);
}

double fPHPair(double u,double t, double E, double mass, const ParamSpaceValues& tpf, const SpaceCoord& psc)   //funcion a integrar
{
	double Nph = tpf.interpolate({ { 0, u } }, &psc	); 
	return Nph
		   *crossSectionBetheHeitler(t)*inelasticityBetheHeitler(t)*t/P2(u);
}

double fBHPair2(double Ee, double Ep, const ParamSpaceValues& tpf, const SpaceCoord& psc)
{
	double gamma_p = Ep / (protonMass*cLight2);
	double Nph = tpf.interpolate({{0,Ee}},&psc);
	double a = 2.0*gamma_p * Ee / (electronMass*cLight2);
	return (pow(a,1.5)*(log(a)-2.0/3.0)+2.0/3.0) * Nph / P2(Ee);
}

double lossesPhotoHadronic(double E, Particle& particle, const ParamSpaceValues& tpf, const SpaceCoord& psc, double phEmin, double phEmax)
{  //E=Ep

	double mass = particle.mass;

	double cte	=	0.5*P2(mass*cLight2)*cLight;

	double b   = phEmax;   //energia maxima de los fotones en erg

	double a1  = mass*cLight2*pionThresholdPH/(2*E);
	double a2  = mass*cLight2*pairThresholdPH/(2*E);

	double a_pi = std::max(a1,phEmin); //(a1,targetPhotonEmin); 
	double a_pa = std::max(a2,phEmin); //targetPhotonEmin);

	double integral = RungeKutta(a_pi, b, &cPionPH, 
	[E,mass](double u){
		return dPH(u,E,mass); 
	}, [&](double u, double t){
		return fPHPion(u,t,E,mass,tpf,psc);
	});

	/*
	integral += RungeKutta(a_pa, b, &cPairPH,
		[E, mass](double u){
		return dPH(u, E, mass);
	}, [E, mass, tpf, &psc](double u, double t){
		return fPHPair(u, t, E, mass, tpf,psc);
	});
	*/
	
	double constant = 7.0*P3(electronRestEnergy)*fineStructConst*thomson*cLight / 
						(9.0*sqrt(2.0)*pi*E*E/protonMass/cLight2);
	double integral2 = constant*RungeKuttaSimple(P2(mass*cLight2)/E,phEmax,[&](double Ee) 
	{return fBHPair2(Ee,E,tpf,psc);});
	
	return cte*integral/E + integral2;
	//return integral2*E;
}

