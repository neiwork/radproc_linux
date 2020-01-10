#include "lossesIC.h"

#include "crossSectionInel.h"
//#include "dataLosses.h"
#include <fmath/RungeKutta.h>
#include <fparameters/parameters.h>
#include <fmath/physics.h>

double cIC(double u)   //limite inferior
{
	return u;
}

double dIC(double u, double E, double mass)   //limite superior     
{
	//DataLosses* data = (DataLosses*)voiddata;
	//const double E = data->E;
	//const double mass = data->mass;

	double s = 4*u*E/P2(mass*cLight2);
	return s*E/(1+s);
}

//double fIC(double u,double t, double E, double mass, fun1 tpf)   //funcion a integrar
double fIC(double u,double t, double E, double mass, const ParamSpaceValues& tpf, const SpaceCoord& psc)   //funcion a integrar
{  

	double Erep = mass*cLight2;
	double r    = t*P2(Erep)/(4*u*E*(E-t));
	
	double Nterm = tpf.interpolate({ { 0, u } }, &psc	); 
	//double Nterm = tpf(u);

	double result = (Nterm/u)*(t-u)*(2*r*log(r)+
       	(1+2*r)*(1-r)+(P2(t/(E-t))*(1-r))/(2*(1+(t/(E-t)))));

	return  result;
}

//double lossesIC(double E, Particle& particle, fun1 tpf, double phEmin, double phEmax)
double lossesIC(double E, Particle& p, const ParamSpaceValues& tpf, const SpaceCoord& psc, double phEmin, double phEmax)
{
	double constant  = 3*crossSectionThomson(p.mass)*P2(p.mass)*pow(cLight,5)/4;
	constant = 3.0*thomson*p.mass*p.mass*pow(cLight,5)/4.0;

	double a  = phEmin;     //energia minima de los fotones en erg
	double b  = phEmax;     //energia maxima de los fotones en erg
	double mass = p.mass;
	double integral = RungeKutta(a,b,&cIC,[mass, E](double u){return dIC(u, E, mass);},
		[&](double u, double t){ return fIC(u, t, E, mass, tpf, psc); });    //le asigno a la variable integral el resultado de la integracion

	double de = constant*integral/P2(E);
	return de;
}

double factorQ(double e1, double e, double Ee)
{
	double gamma_e = 4.0*e*Ee/P2(electronRestEnergy);
	double q = e1/(gamma_e*(Ee-e1));
	return 2*q*log(q)+(1+2*q)*(1-q)+0.5*P2(gamma_e*q)*(1-q)/(1+gamma_e*q);
}

double dN_dtde1(double e1, double Ee, Particle& p, const ParamSpaceValues& tpf,const SpaceCoord& psc,
				double logphEmin, double phEmax)
{
	double g = Ee/electronRestEnergy;
	double Emin = 0.25*e1*electronRestEnergy/(g*(Ee-e1));
	return integSimpson(log(Emin),log(e1),[&e1,&Ee,&tpf,&psc](double loge)
			{
				double e = exp(loge);
				return tpf.interpolate({{0,e}},&psc) * factorQ(e1,e,Ee);
			},40);
}
double lossesICnew(double E, Particle& p, const ParamSpaceValues& tpf, const SpaceCoord& psc,
					double phEmin, double phEmax)
{
	double constant = 2.0*pi*P2(electronRadius*electronRestEnergy)*cLight/(E*E);
	return constant*integSimpson(log(phEmin),log(min(0.99*E,phEmax)),[E,phEmin,phEmax,&p,&tpf,&psc](double loge1)
				{
					double e1 = exp(loge1);
					return e1*e1*dN_dtde1(e1,E,p,tpf,psc,phEmin,phEmax);
				},50);
}

double lossesIC_Th(double E, Particle& p, const ParamSpaceValues& tpf, const SpaceCoord& psc,
					double phEmin, double phEmax)
{
	double g = E / (p.mass*cLight2);
	double Uph = integSimpson(phEmin,phEmax,[tpf,&psc](double Eph)
					{return tpf.interpolate({{0,Eph}},&psc)*Eph;},100);
	return 4.0/3.0 * thomson * cLight * Uph * P2(electronMass/p.mass) * g*g;
}
