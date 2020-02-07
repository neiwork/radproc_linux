#include "lossesPhotoHadronic.h"

#include "crossSectionInel.h"
#include <fmath/RungeKutta.h>
#include <fparameters/parameters.h>
#include <finjection/pgammaPionInj.h>
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

double dPH(double u, double g)   //limite superior     
{
	return 2*u*g;   //d1 = 2*u*ep/(masa*cluz**2)
}

//double fPHPion(double u,double t, double E, double mass, fun1 tpf)   //funcion a integrar
double fPHPion(double u,double t, const ParamSpaceValues& tpf, const SpaceCoord& psc)   //funcion a integrar
{
	double Nph = tpf.interpolate({ { 0, u } }, &psc	); 
	return Nph*crossSectionPHPion(t)*inelasticityPHPion(t)*t/P2(u);
}

double fPHPair(double u,double t, const ParamSpaceValues& tpf, const SpaceCoord& psc)   //funcion a integrar
{
	double Nph = tpf.interpolate({ { 0, u } }, &psc	); 
	return Nph*crossSectionBetheHeitler(t)*inelasticityBetheHeitler(t)*t/P2(u);
}

double fBHPair2(double Ee, double Ep, const ParamSpaceValues& tpf, const SpaceCoord& psc)
{
	double gamma_p = Ep / (protonMass*cLight2);
	double Nph = tpf.interpolate({{0,Ee}},&psc);
	double a = 2.0*gamma_p * Ee / (electronMass*cLight2);
	return (pow(a,1.5)*(log(a)-2.0/3.0)+2.0/3.0) * Nph / P2(Ee);
}

double lossesPhotoHadronic(double E, Particle& p, const ParamSpaceValues& tpf, 
							const SpaceCoord& psc, double phEmin, double phEmax)
{  //E=Ep
	double g = E/(p.mass*cLight2);
	double cte	=	0.5*cLight/(g*g);
	double a1  = pionThresholdPH/(2*g);
	double a2  = pairThresholdPH/(2*g);

	double a_pi = std::max(a1,phEmin); //(a1,targetPhotonEmin); 
	double a_pa = std::max(a2,phEmin); //targetPhotonEmin);

	double integral = RungeKutta(a_pi, phEmax, &cPionPH, 
	[&](double u){return 2.0*u*g;},
	[E,tpf,&psc](double u, double t){return fPHPion(u,t,tpf,psc);});

	/*
	integral += RungeKutta(a_pa, b, &cPairPH,
		[E, mass](double u){
		return dPH(u, E, mass);
	}, [E, mass, tpf, &psc](double u, double t){
		return fPHPair(u, t, E, mass, tpf,psc);
	});
	*/
	
	/* double constant = 7.0*P3(electronRestEnergy)*fineStructConst*thomson*cLight / 
						(9.0*sqrt(2.0)*pi*E*E/protonMass/cLight2);
	double integral2 = constant*RungeKuttaSimple(P2(mass*cLight2)/E,phEmax,[&](double Ee) 
	{return fBHPair2(Ee,E,tpf,psc);}); */
	
	return cte*integral*E;
	//return integral2*E;
}

double lossesPhotoMeson(double E, Particle& p, const ParamSpaceValues& tpf, 
							const SpaceCoord& psc, double phEmin, double phEmax)
{
	double g = E/(p.mass*cLight2);
	double cte	=	0.5*cLight/(g*g);
	double a1  = pionThresholdPH/(2*g);
	double a_pi = std::max(a1,phEmin);

	double integral = RungeKutta(a_pi, phEmax, &cPionPH, 
	[&](double u){return 2.0*u*g;},
	[tpf,&psc](double u, double t){return fPHPion(u,t,tpf,psc);});

	return cte*integral*E;
}

double lossesPhotoPair(double E, Particle& p, const ParamSpaceValues& tpf, 
							const SpaceCoord& psc, double phEmin, double phEmax)
{
	double g = E/(p.mass*cLight2);
	double cte	=	0.5*cLight/(g*g);
	double a2  = pairThresholdPH/(2*g);
 
	double a_pa = std::max(a2,phEmin);

	double integral = RungeKutta(a_pa, phEmax, &cPairPH,
		[&](double u){return 2.0*u*g;},
		[tpf,&psc](double u, double t){
		return fPHPair(u, t, tpf,psc);
	});
	
	return cte*integral*E;
}

double lossesPhotoHadronic_simple(double E, Particle& p, const ParamSpaceValues& tpf, const SpaceCoord& psc,
							double phEmin, double phEmax)
{
	return E * t_pion_PHsimple(E,p,tpf,psc,phEmin,phEmax);
}