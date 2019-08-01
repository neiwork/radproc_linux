#include "pairInjection.h"


#include <fparameters/parameters.h>
#include <fmath/RungeKutta.h>
#include <fmath/physics.h>
//#include <algorithm>

double cAnnihilation(double x,double E)
{
	const double Erep = electronMass*cLight2;
	
	double inf = x*P2(Erep) / (4.0*E*(x-E));  

	return inf;   
}

/*double dAnnihilation(double x )
{
	return targetPhotonEmax;        //este es el infinito del limite superior para la integral en Eph
}*/

double fAnnihilation(double x, double y, double E, const ParamSpaceValues& ntPh, const ParamSpaceValues& tpf, const SpaceCoord& distCoord, double tpEmin, double tpEmax)   
{ 
	//x=Ega; y=Eph
	double photonDist_x(0.0), photonDist_y(0.0);

	if (x >= tpEmin && x <= tpEmax){
		photonDist_x = ntPh.interpolate({ { 0, x } }, &distCoord); //tpf(x);
	}
	double nph;
	if (y >= tpEmin && y <= tpEmax){
		photonDist_y = tpf.interpolate({ { 0, y } }, &distCoord); //tpf(y);
	}

	const double Erest = electronMass*cLight2;
                      
	double result  =   ((photonDist_x*photonDist_y)/(P2(y)*P3(x)))*     
					   ((4.0*P2(x)/(E*(x-E)))*log(4.0*E*y*(x-E)/(Erest*Erest*x))
					   -8.0*x*y/(Erest*Erest) + 2.0*(2.0*x*y-P2(Erest))*P2(x/Erest)/(E*(x-E))
					   -(1.0-P2(Erest)/(x*y))*pow(x,4)/(P2(E*(x-E))));
//&& (y > 2*P2(Erest)/x)
	return ((y<Erest) && (y > 2*P2(Erest)/x)  ) ? result : 0.0;   //pido que de algo solo si epsilon < Erep   //&& Erest<x
}


double pairInjection(double E, const ParamSpaceValues& ntPh, const ParamSpaceValues& tpf, const SpaceCoord& distCoord, double tpEmin, double tpEmax)  //E=Ee
{	
	using std::bind; using namespace std::placeholders; // para _1, _2, etc.

	double cte = 3.0*cLight*thomson*pow((electronMass*cLight2),4)/32.0;

	double sup = 10.0*tpEmax;  //este es el infinito del limite superior para la integral en Egamma
   
	double inf = E;  //Ega_min < Ee_min  --> la condicion esta asegurada


	double integral  = RungeKutta(inf,sup, 
		[E,tpEmin](double u) {return max(tpEmin,cAnnihilation(u, E)); },  //limite inferior
		[tpEmax](double u) {return tpEmax; },							  //limite superior
		[E, tpf, ntPh, &distCoord, tpEmin, tpEmax](double u, double t) 
		{return fAnnihilation(u, t, E, ntPh, tpf, distCoord, tpEmin, tpEmax); });
		
	double emissivityA = cte*integral;

	return emissivityA;

}