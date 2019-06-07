#include "thermalDistribution.h"


#include "maxwellJuttner.h"
#include "modelParameters.h"

//#include <fmath/RungeKutta.h>
#include <fparameters/Dimension.h>
#include <fparameters/SpaceIterator.h>
#include <fparameters/parameters.h>

#include <fmath/physics.h>

#include <boost/math/special_functions/bessel.hpp>

void thermalDistribution(Particle& p, State& st)
{
	double Erest = p.mass*cLight2;
	
	double gMin = p.emin()/Erest;
	double gMax = p.emax()/Erest;
	
	p.ps.iterate([&](const SpaceIterator& i) {
	//p.distribution.fill ([&](const SpaceIterator& i){
		
		//const double E = i.val(DIM_E);
		double r = i.val(DIM_R);
		
		double norm_temp, A;

		if(p.id == "electron"){
			norm_temp = boltzmann*st.tempElectrons.get(i)/(Erest);
			A = st.denf_e.get(i); 
		}
		else if(p.id == "proton"){
			norm_temp = boltzmann*st.tempIons.get(i)/(Erest);
			A = st.denf_i.get(i); 
		}
		
		//double K2 =  boost::math::cyl_bessel_k(2, 1.0/norm_temp); //bessk(2, 1.0/norm_temp);
		
		
		p.ps.iterate([&](const SpaceIterator& j) {
				
			const double E = j.val(DIM_E);
			double g = E/Erest;
			
			double result =  maxwellRel(g, norm_temp, A/Erest);
			
			p.distribution.set(j,result); //en unidades de cm^-3 erg^-1
						
		}, { -1, i.coord[DIM_R]} );
		
	}, { 0, -1} );
	
	
}