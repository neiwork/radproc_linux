#include "luminosityHadronic.h"


#include <fparameters/parameters.h>
#include <fmath/RungeKutta.h>
#include <fmath/interpolation.h>
#include <flosses/crossSectionInel.h>
#include <fmath/physics.h>
#include <algorithm>

double fHadron(double x, const Particle& p,
	const double density, const SpaceCoord& psc) //funcion a integrar   x=Ecreator; L=L(Ega)
{	
	double Kpi = 0.17;
	double eval = p.mass*cLight2+x/Kpi;
	double Ekin = x/Kpi;
	double distCreator=0.0;
	if (x < p.emax() && x < p.emin()) {
		distCreator = p.distribution.interpolate({ { 0, x } }, &psc);
	}
	
	double thr = 0.0016; //1GeV
	//double sigma = 30e-27*(0.95+0.06*log(Ekin/thr));
	
	double l = log10((protonMass*cLight2+x/Kpi)/1.6);
	double sigma = 1.e-27 * (34.3+1.88*l+0.25*l*l);

	//const double density = denfdensity.get(psc);

	double pionEmiss = cLight*density*sigma*distCreator/Kpi;  //sigma = crossSectionHadronicDelta(Ekin)
							//lo saco asi pongo la condicion Ekin > Ethr en el limite de la int
	double result = pionEmiss/sqrt(P2(x)-P2(neutralPionMass*cLight2));
	return result;
}


double luminosityHadronic(double E, const Particle& creator,
	const double density, const SpaceCoord& psc)
{
	double Kpi = 0.17;
	double thr = 0.0016; //1GeV

	double Max  = 1.6e-12*pow(10.0,17.0);   //esto es un infinito 
	double Min  = std::max(E+P2(neutralPionMass*cLight2)/(4*E),thr*Kpi); //== Ekin > Ethr

	double integral = RungeKuttaSimple(Min, Max, 
		[&](double x) {return fHadron(x, creator, density, psc); }
	);    //integra entre Emin y Emax

	double luminosity = integral*E*planck*0.25/pi; // [erg s^-1 Hz^-1 cm^-3]
	return luminosity; //P2(Dlorentz);  

	//Dlorentz es el factor que transforma las distribuciones en el caso de jets
	// con /P2(Dlorentz) paso del sist de lab al sist comoving con el jet ;   
}