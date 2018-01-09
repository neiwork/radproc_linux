#include "thermalCompton.h"

#include "modelParameters.h"
#include <fparameters/parameters.h>
#include <fmath/RungeKutta.h>
#include <fmath/physics.h>


double cICt(double u, double E, double gMin)  //limite inferior
{

	double inf = 0.5*sqrt(E/u);

	return std::max(gMin,inf);  
}

double dICt(double u, double gMax)//double E)         //limite superior   
{          

	return gMax; //E;  //esta es la condicion epsilon < Ega	                                 
}



double fICt(double u, double t, double E, const Particle& creator, const SpaceCoord& distCoord, const ParamSpaceValues& tpf)   //funcion a integrar  u=Ee
{    
	
	
	double Erest = creator.mass*cLight2;
	double Eval = t*Erest;
	double distCreator;
	if (Eval < creator.emin() || Eval> creator.emax()){
		distCreator = 0.0;
	}
	else{
		distCreator = creator.distribution.interpolate({ { 0, Eval } }, &distCoord); 
	}
	  

	double Nph = tpf.interpolate({ { 0, u*Erest } }, &distCoord);  

	double f = 2.0*E*log(E/(4.0*P2(t)*u)) + E + 4.0*P2(t)*u - P2(E)/(2.0*P2(t)*u);

	double function = (distCreator*Erest/(2.0*pow(t,4.0)))*(Nph*Erest/P2(u))*f*3.0*thomson*cLight;


	return function;
}


double thCompton(double E, const Particle& creator, const SpaceCoord& distCoord, const ParamSpaceValues& tpf, double phEmin)
{
	
	double Erest = creator.mass*cLight2;
	double E_norm = E/Erest;
	
	double gMin = creator.emin()/Erest;
	double gMax = creator.emax()/Erest;
	double eMin = phEmin/Erest;
	double eMax = E/Erest;
	//double mass = creator.mass;
	
	//double cte  = 3.0*crossSectionThomson(creator.mass)*P2(creator.mass)*pow(cLight,5)/4;

	double integral = RungeKutta(eMin,eMax,
		[E_norm,gMin](double u){
			return cICt(u,E_norm,gMin);
		}, 
		[gMax](double u){
			return dICt(u,gMax);
		}, 
		[E_norm,&creator,&distCoord, tpf](double u, double t){
			return fICt(u, t,E_norm,creator, distCoord, tpf); 
		});    //le asigno a la variable integral el resultado de la integracion   

	double luminosity = integral*E; //cte*P2(E);

	return luminosity;

	}

