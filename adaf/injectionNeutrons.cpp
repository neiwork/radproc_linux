#include <finjection/neutronPp.h>
#include <finjection/neutronPgamma.h>

#include "injectionNeutrons.h"
#include "messages.h"
#include "globalVariables.h"
#include "adafFunctions.h"
#include <fparameters/parameters.h>
#include <fparameters/SpaceIterator.h>
#include <fparameters/Dimension.h>
#include <fmath/fbisection.h>
#include <fmath/RungeKutta.h>
#include <fmath/physics.h>
#include <boost/property_tree/ptree.hpp>
#include <iostream>

void injectionNeutrons(Particle& n, Particle& p, State& st)
{
	n.ps.iterate([&](const SpaceIterator& i) {
		double dens = st.denf_i.get(i);
		const double rB1 = i.val(DIM_R)/sqrt(paso_r);
		const double rB2 = rB1*paso_r;
		const double thetaH = st.thetaH.get(i);
		const double vol  = (rB2*rB2*rB2-rB1*rB1*rB1)*(4.0/3.0)*pi*cos(thetaH);;
		n.ps.iterate([&](const SpaceIterator& jE) {
			const double E = jE.val(DIM_E);
			double total = neutronPp(E,p,dens,i);
			n.injection.set(jE,total); //en unidades de erg^-1 s^-1 cm^{-3}
		},{-1,i.coord[DIM_R],0});
	},{0,-1,0});
}