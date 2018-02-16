#include "photonEscape.h"

#include "modelParameters.h"
#include <flosses/crossSectionInel.h>
#include <fparameters/SpaceIterator.h>
#include <fparameters/parameters.h>
#include <fmath/physics.h>

#include <boost/property_tree/ptree.hpp>


double photonEscape_f(double x)
{
	if (0.1 < x && x < 1.0){ return (1.0 - x) / 0.9; }
	else if (x <= 0.1){ return 1.0; }
	else { return 0.0; }
}

double escapePhoton(double Eph, Particle& electron, ParamSpaceValues& denf, const SpaceIterator& i) //en [s]
{

	double rg = GlobalConfig.get<double>("rg");

	double rMax = GlobalConfig.get<double>("rmax");

	//st.photon.ps.iterate([&](const SpaceIterator& i) { 
	
		double average = 0.0;   
		
		electron.ps.iterate([&](const SpaceIterator& j) {  //integro sobre la energia de los electrones
			
			double E = j.val(DIM_E);
			double gamma = E/(electronMass*cLight2);
			double beta = sqrt(1.0-1.0/P2(gamma));
			
			average = average + angleAveragedKN(Eph/(electronMass*cLight2))*beta;
			
		}, { -1, i.coord[1], i.coord[2] } );
	
		double radius = (rMax - i.val(DIM_R))*rg;
		double opticalDepth = 2.0*radius*average*denf.get(i);
	
	
	double result = radius*(1.0+opticalDepth*photonEscape_f(Eph/(electron.mass*cLight2)))/cLight;
	
	return result;
	


}
