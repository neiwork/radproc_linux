#include "distribution.h"


#include "maxwellJuttner.h"
#include "functions.h"
#include "modelParameters.h"

#include <fparameters/Dimension.h>
#include <fparameters/SpaceIterator.h>
#include <fparameters/parameters.h>
//#include <boost/property_tree/ptree.hpp>
#include <fmath/physics.h>

void thermalDistribution(Particle& p, State& st)
{

	p.distribution.fill ([&](const SpaceIterator& i){
		
		const double E = i.val(DIM_E);
		double r = i.val(DIM_R);
		double theta = i.val(DIM_THETA);

		double temp = temp_e(r, theta);
		return maxwellRel(E, temp, p.mass);
	});
	
}