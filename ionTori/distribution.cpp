#include "distribution.h"


#include "maxwellJuttner.h"
#include "functions.h"
#include "modelParameters.h"

#include <fmath/RungeKutta.h>
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
		double theta = i.val(DIM_THETA);

		double norm_temp = boltzmann*st.tempElectrons.get(i)/(Erest);
		
		double K2 =  boost::math::cyl_bessel_k(2, 1.0/norm_temp); //bessk(2, 1.0/norm_temp);
		
		double integral = RungeKuttaSimple(gMin,gMax,
			[&norm_temp](double g){
			return f_norm(g, norm_temp); 
		});
		
		double A=0.0;
		if (integral > 0.0)
		{
			A = st.denf_e.get(i)/integral;  //*theta*K2
		}
		
		p.ps.iterate([&](const SpaceIterator& j) {
				
			const double E = j.val(DIM_E);
			double g = E/Erest;
			
			double result =  maxwellRel(g, norm_temp, A/Erest);
			
			p.distribution.set(j,result); //en unidades de cm^-3 erg^-1
			
		//	return ;  //paso de N(g) -> N(E)  [N(E)] = erg^-1
			
		}, { -1, i.coord[DIM_R], i.coord[DIM_THETA] } );
		
	}, { 0, -1, -1 } );
	
	
}