#include "lossesSyn.h"


#include <fmath/RungeKutta.h>
#include <fmath/physics.h>


double lossesSyn(double E, double magf, Particle& p)
{
	double gamma = E / (p.mass*cLight2);
	double wmag = magf*magf / (8.0 * pi);
	return (4.0/3.0)*thomson*cLight*wmag*P2(electronMass/p.mass) * gamma*gamma;
}	

/////////////////////////////////////////////////////////////////
double fSynLosSec(double x, double magneticField, double E)         //funcion a integrar   x=Ega
{
	double mass = electronMass;
	double cte = pow(3.0, 0.5)*P3(electronCharge)*magneticField / (planck*mass*cLight2);

	double Echar = 3 * electronCharge*planck*magneticField*P2(E) / (4 * pi*P3(mass)*cLight*P2(cLight2));

	double tau = x / E;
	double aux = x / (Echar*(1 - tau));  //aca el aux es el x real

	double result = cte*1.85*(1 - tau)*pow(aux, (1.0 / 3.0))*exp(-aux);

	return tau < 1 ? result : 0.0;
} 
///////////////////synchr losses for secondary pairs

double lossesSynSec(double E, double magneticField, Particle& particle)
{
	double EmaxL = 1.6e-12*pow(10.0,17.0);  
	double EminL = 1.6e-12*pow(10.0,-2.0); 

	double result = RungeKuttaSimple(EminL, EmaxL, [&E, &magneticField](double x){
		return fSynLosSec(x, magneticField, E);
	});  //integra entre Emin y Emax
	
	return result;
}
