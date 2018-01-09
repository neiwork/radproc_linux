#include "maxwellJuttner.h"

#include <fmath/mathFunctions.h>
#include <fmath/physics.h>
#include <boost/math/special_functions/bessel.hpp>


double f_norm(double g, double norm_temp)
{
	double beta = sqrt(1.0-1.0/P2(g));
	
	return P2(g)*beta*exp(-g/norm_temp); 
	
}

double maxwellRel(double gamma, double norm_temp, double norm)
{
		
	double beta = sqrt(1.0-1.0/P2(gamma));
	
	double K2 =  boost::math::cyl_bessel_k(2, 1.0/norm_temp); //bessk(2, 1.0/norm_temp);

	//double dist_g = P2(gamma)*beta/(norm_temp*K2*exp(gamma/norm_temp));
	
	double dist_g = P2(gamma)*beta*exp(-gamma/norm_temp);

	
//	if( K2 > 0.0 )
	//{
		return  norm*dist_g;
	//}
	//else
	//{ 
	//	return 0.0;
	//}
}
