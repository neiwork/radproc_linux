#include "luminositySynchrotron.h"

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_synchrotron.h>
#include "opticalDepthSSA.h"
#include <fparameters/parameters.h>
#include <fmath/RungeKutta.h>
#include <fmath/physics.h>


double fSyn(double x, double E, const Particle& creator, const ParamSpaceValues& magf, const SpaceCoord& psc)         //funcion a integrar   x=Ee; L=L(Ega)
{

	const double magneticField = magf.get(psc);

	double distCreator;
	if (x < creator.emin() || x > creator.emax()){
		distCreator = 0.0;
	}
	else{
		distCreator = creator.distribution.interpolate({ { 0, x } }, &psc); 
	}
	double Erest = creator.mass*cLight2;
	double cte = sqrt(3.0)*P3(electronCharge)*magneticField / (planck*Erest);
	double Echar = 3.0*electronCharge*planck*magneticField*P2(x/Erest) / 
					(4.0*pi*creator.mass*cLight);
	
	double aux = E/Echar;  //aca el aux es el x real

	double result = cte*1.85*distCreator*pow(aux,1.0/3.0)*exp(-aux);  

	return result;
}


double luminositySynchrotron(double E, const Particle& c, const SpaceCoord& psc, const ParamSpaceValues& magf)
{
	
	double integralS = RungeKuttaSimple(c.emin(), c.emax(), [&](double e){
		return fSyn(e, E, c, magf, psc);
	});

	double luminosityS = integralS*E; //multiplico por E asi obtengo luminosidad

	if (luminosityS > 0.0){ return luminosityS; }
	else { return 0.0; }

}



double fSyn2(double x, double E, const Particle& creator, double magf, const SpaceCoord& psc)         //funcion a integrar   x=Ee; L=L(Ega)
{
	double distCreator;
	if (x < creator.emin() || x > creator.emax()){
		distCreator = 0.0;
	}
	else{
		distCreator = creator.distribution.interpolate({ { 0, x } }, &psc); 
	}
	double Erest = creator.mass*cLight2;
	double cte = sqrt(3.0)*P3(electronCharge)*magf / (planck*Erest);
	double Echar = 3.0*electronCharge*planck*magf*P2(x/Erest) / 
					(4.0*pi*creator.mass*cLight);
	
	double aux = E/Echar;  //aca el aux es el x real

	double result = cte*1.85*distCreator*pow(aux,1.0/3.0)*exp(-aux);

	return result;
}
double luminositySynchrotron2(double E, const Particle& c, const SpaceCoord& psc, double magf)
{
	double integralS = RungeKuttaSimple(c.emin(), c.emax(), [&](double e){
		return fSyn2(e, E, c, magf, psc);
	});

	double luminosityS = integralS*E; //multiplico por E asi obtengo luminosidad

	if (luminosityS > 0.0){ return luminosityS; }
	else { return 0.0; }

}

double fSyn3(double Ee, double E, const Particle& creator, double magf, const SpaceCoord& psc)
{
	double distCreator = (Ee > creator.emin() && Ee < creator.emax()) ? 
				creator.distribution.interpolate({{0,Ee}},&psc) : 0.0;
	double gamma_e = Ee/(electronMass*cLight2);
	double muMin = -0.999;
	double muMax = 0.999;
	size_t nMu = 30;
	double dMu = (muMax-muMin)/nMu;
	double sum = 0.0;
	double mu = muMin;
	double Constant = 3.0*planck*electronCharge*magf*gamma_e*gamma_e / (4.0*pi*electronMass*cLight);
	for (size_t jMu=0;jMu<nMu;jMu++) {
		double Ec = Constant * sqrt(1.0-mu*mu);
		double x = E/Ec;
		double integ = (x < 10.0) ? integSimpson(log(x),log(10.0),[&](double logxp)
						{
							double xp = exp(logxp);
							return gsl_sf_bessel_Knu(5.0/3.0,xp)*xp;
						},50) : 0.0;
		sum += dMu*integ;
		mu += dMu;
	}
	return distCreator/(gamma_e*gamma_e) * sum;
}

double phiSy(double x, double si)
{
	return (x < 800) ? si*gsl_sf_synchrotron_1(x) : 0.0;
}


double luminositySynchrotronExact(double Eph, const Particle& c, const SpaceCoord& psc, double magf)
{
	double integralS = integSimpsonLog(c.emin(),c.emax(),[Eph,magf,&c,&psc](double Ee)
			{
				double g = Ee/(c.mass*cLight2);
				double integAng = integSimpson(0.0,0.9999,[&](double mu)
				//double integAng = qromo(0.0,0.9999,[&](double mu)
						{
							double si = sqrt(1.0-mu*mu);
							double Ec = 3.0*electronCharge*planck*magf*si*g*g/(4*pi*c.mass*cLight);
							double x = Eph/Ec;
							return phiSy(x,si);
						//});
						},10);
				return c.distribution.interpolate({{0,Ee}},&psc)*integAng;
			},50);
	double constant = sqrt(3.0)*gsl_pow_3(electronCharge)*magf/(planck*c.mass*cLight2);
	return (integralS > 0.0) ? constant*integralS*Eph : 0.0;
}


/*
double luminositySynchrotron_conSSA(double E, const Particle& creator)
{

	double Emax = creator.emax();
	double Emin = creator.emin();

	double tau = opticalDepthSSA(E, creator.mass, Emin, Emax, creator);  //E=Eph

	double factorSSA = (1.0-exp(-tau))/tau;

	if (factorSSA > 1.0 || factorSSA == 0.0) //1.0e-3)  //lo cambio solo en el caso que interesa
	{ factorSSA = 1.0; }	
	
	//double integralS = RungeKuttaSimple(Emin, Emax, bind(fSyn,_1,E,creator));

	double integralS = RungeKuttaSimple(creator.emin(), creator.emax(),   //RungeKuttaSimple(double a, double b, fun1 f)
		[E, &creator](double x){
		return fSyn(x, E, creator);  //double fSyn(double x, double E, const Particle& creator)
	});


	double luminosityS = factorSSA*integralS*E; //multiplico por E asi obtengo luminosidad

	if (luminosityS > 0.0){return luminosityS ; }
	else {return 0.0;} 
	 
}*/


////////////////////////////////////////////////////////////////////////////////////////////////////
/*
double fSynS(double x, double E, const Particle& creator)         //funcion a integrar   x=Ee; L=L(Ega)
{
	double distCreator = creator.dist(x);//interpol(x,Ecreator,Ncreator,Ncreator.size()-1);

	double cte	= pow(3.0,0.5)*P3(electronCharge)*magneticField/(planck*creator.mass*cLight2);

	double Echar = 3*electronCharge*planck*magneticField*P2(x)/(4*pi*P3(creator.mass)*cLight*P2(cLight2));
	
	double tau = E/x;
	double aux = E/(Echar*(1-tau));  //aca el aux es el x real

	double result = cte*1.85*(1-tau)*distCreator*pow(aux,(1.0/3.0))*exp(-aux);  

	return result;    // esta condicion la puse en el limite inferior tau<1 ? result : 0.0;
}



double luminositySynchrotronSec(double E, const Particle& c)
{
	double Emax = c.emax();
	double Emin = c.emin();
	double inf  = std::max(Emin,E);   //esto lo agrego asi le saco la condicion sobre tau < 1

	double tau = opticalDepthSSA(E, c.mass, inf, Emax, c);  //E=Eph

	double factorSSA = 1.0;

	if (tau > 1.0e-3)  //lo cambio solo en el caso que interesa
	{ factorSSA = (1.0-exp(-tau))/tau;}

	double integralS = RungeKuttaSimple(Emin, Emax, [&E, &c](double e){
		return fSynS(e, E, c);
	});

	double luminosityS = factorSSA*integralS*E*volume; //multiplico por E asi obtengo luminosidad
	                                 //divido por E asi obtengo emisividad y no luminosidad

	return luminosityS ; 

}

*/