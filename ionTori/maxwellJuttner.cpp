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
		
	double beta = sqrt(1.0-1.0/(gamma*gamma));
	
	double K2 =  boost::math::cyl_bessel_k(2, 1.0/norm_temp); //bessk(2, 1.0/norm_temp);
	
	double dist_g = gamma*gamma*beta*exp(-gamma/norm_temp) / (K2*norm_temp);

	return  dist_g*norm;

}
