#include "maxwellJuttner.h"

#include <fmath/mathFunctions.h>
#include <fmath/physics.h>


double maxwellRel(double E, double T, double mass)
{
	double Erest = mass*cLight2;
	
	double theta = boltzmann*T/Erest;
	
	double gamma = E/Erest;
	
	double beta = sqrt(1.0-1.0/P2(gamma));
	
	double K2 =  bessk(2, 1.0/theta);

	double dist_g = P2(gamma)*beta/(theta*K2*exp(gamma/theta));

	return dist_g/Erest;  //paso de N(g) -> N(E)  [N(E)] = erg^-1
}
