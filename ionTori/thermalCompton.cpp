#include "thermalCompton.h"

#include "modelParameters.h"
#include <fparameters/parameters.h>
#include <fmath/RungeKutta.h>
#include <fmath/physics.h>



double eta_1(double P, double A)
{
	return P*(A-1.0)/(1.0-A*P);
	
}

double eta_3(double P, double A)
{
	return -1.0-log(P)/log(A);
}
		
double eta_2(double P, double A)
{
	return pow(3.0,-eta_3(P,A))*eta_1(P,A);
}


double comptBremss(double E, double theta_e, double r, double theta, const SpaceCoord& distCoord, 
						ParamSpaceValues& denf,const ParamSpaceValues& jBr)
{
	
	double lim_inf = r;
	double lim_sup = r*1000.0; // ver como obtengo el r_out desde aca
	
	double tau = RungeKuttaSimple(lim_inf, lim_sup, 
	[&E,&r,&theta,&denf,&distCoord](double r) 
		{return denf.interpolate({{ DIM_R, r},{ DIM_THETA, theta} },&distCoord)* 4.0*pi/P2(E); } );  

			
	double P = 1.0-exp(-tau);
	
	double A = 1.0+4.0*theta_e+16.0*P2(theta_e);
	
	
	double e1 = eta_1(P,A);
	double e3 = eta_3(P,A);
	
	double xc = E/(electronMass*cLight2);
	//distCoord.dims;
	
	double jBr_i = jBr.get(distCoord);
	
	double res = jBr_i*3.0*e1*theta_e*((1.0/3.0-xc/(3.0*theta_e)) - 1.0/(e3+1.0)*( pow(3.0,-e3-1.0) -pow(xc/(3.0*theta_e) , e3+1.0) ) );
	
	
	
	}