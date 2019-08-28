#include "pairInjection.h"


#include <fparameters/parameters.h>
#include <fmath/RungeKutta.h>
#include <fmath/physics.h>
//#include <algorithm>

double cAnnihilation(double x,double E)
{
	const double Erep = electronMass*cLight2;
	return x*Erep*Erep / (4.0*E*(x-E));
}

double cAnnihilationPrueba(double x,double g)
{
	return x / (4.0*g*(x-g));
}

/*double dAnnihilation(double x )
{
	return targetPhotonEmax;        //este es el infinito del limite superior para la integral en Eph
}*/

double fAnnihilation(double x, double y, double E, const ParamSpaceValues& ntPh, const ParamSpaceValues& tpf, const SpaceCoord& distCoord, double tpEmin, double tpEmax)   
{ 
	//x=Eg; y=Eph
	const double Erest = electronMass*cLight2;
	double photonDist_x(0.0), photonDist_y(0.0);
	double secondInfLimit = Erest*Erest * x / (4.0*E*(x-E));

	photonDist_x = (x > tpEmin && x < tpEmax) ?
		ntPh.interpolate({ { 0, x } }, &distCoord) : 0.0;
	photonDist_y = (y > secondInfLimit && y < tpEmax) ?
		tpf.interpolate({ { 0, y } }, &distCoord) : 0.0;

	double result  =   photonDist_x*photonDist_y / (x*x*x*y*y) *     
					   ( (4.0*x*x / (E*(x-E))) * log(4.0*E*y*(x-E)/(Erest*Erest*x))
					   - 8.0*x*y/(Erest*Erest) + 2.0*(2.0*x*y-Erest*Erest)*P2(x/Erest)/(E*(x-E))
					   - (1.0-Erest*Erest/(x*y)) * pow(x,4)/(P2(E*(x-E))));
	return (y < Erest && y > 2*P2(Erest)/x ) ? result : 0.0;   //pido que de algo solo si epsilon < Erep   //&& Erest<x
}

double fAnnihilationPrueba(double eg,double om,double g,const ParamSpaceValues& tpf,
							const SpaceCoord& distCoord,double tpEmin, double tpEmax)   
{ 
	double electronRestEnergy = electronMass*cLight2;
	double photonDist_x(0.0);
	if (om >= tpEmin/electronRestEnergy && om <= tpEmax/electronRestEnergy)
		photonDist_x = tpf.interpolate({{0,om*electronRestEnergy}},&distCoord);

	double result = photonDist_x/(om*om) *     
					   ( 4.0*eg*eg/(g*(eg-g)) * log(4.0*g*om*(eg-g)/eg)
					   - 8.0*eg*om + 2.0*eg*eg*(2.0*eg*om-1.0)
					   - (1.0-1.0/(eg*om))*eg*eg*eg*eg/(g*g*(eg-g)*(eg-g)) );
	return (om < 0.1) ? result : 0.0;
}


double pairInjection(double E, const ParamSpaceValues& ntPh, const ParamSpaceValues& tpf,
					const SpaceCoord& distCoord, double tpEmin, double tpEmax)
{	
	using std::bind; using namespace std::placeholders; // para _1, _2, etc.
	double cte = 3.0*cLight*thomson*pow((electronMass*cLight2),4)/32.0;
	double sup = tpEmax;  //este es el infinito del limite superior para la integral en Egamma
	double inf = E;  //Ega_min < Ee_min  --> la condicion esta asegurada

	double integral  = RungeKutta(inf,sup, 
		[E,tpEmin](double Eg) {return max(tpEmin,cAnnihilation(Eg, E)); },  //limite inferior
		[tpEmax](double Eg) {return tpEmax;},							  //limite superior
		[E,tpf,ntPh,&distCoord,tpEmin,tpEmax](double Eg, double Eph) 
		{return fAnnihilation(Eg,Eph,E,ntPh,tpf,distCoord,tpEmin,tpEmax); });
		
	double emissivityA = cte*integral;
	return emissivityA;
}

double pairInjectionPrueba(double E, const ParamSpaceValues& ntPh, const ParamSpaceValues& tpf,
					const SpaceCoord& distCoord, double tpEmin, double tpEmax)
{	
	double electronRestEnergy = electronMass*cLight2;
	double gamma = E / electronRestEnergy;
	double cte = 3.0*cLight*thomson*electronRestEnergy/32.0;
	double sup = tpEmax/electronRestEnergy;  //este es el infinito del limite superior para la integral en Egamma

	double integral  = RungeKuttaSimple(gamma,sup,[&](double eg){ return 
						RungeKuttaSimple(cAnnihilationPrueba(eg,gamma),tpEmax/electronRestEnergy,[&](double om) { return
						((eg >= tpEmin/electronRestEnergy && eg <= tpEmax/electronRestEnergy) ?
						ntPh.interpolate({{0,eg*electronRestEnergy}},&distCoord) : 0.0) / (eg*eg*eg) *
						fAnnihilationPrueba(eg,om,gamma,tpf,distCoord,tpEmin,tpEmax);});});
		
	double emissivityA = cte*integral;
	return emissivityA;
}