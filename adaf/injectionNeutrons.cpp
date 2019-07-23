#include "injectionNeutrons.h"

#include "messages.h"
#include "globalVariables.h"
#include "adafFunctions.h"
//#include <fmath/fbisection.h>
//#include <fmath/RungeKutta.h>
//#include <fmath/physics.h>
//#include <boost/property_tree/ptree.hpp>
#include <iostream>
#include <finjection/neutronPp.h>
#include <finjection/neutronPgamma.h>
#include <fparameters/SpaceIterator.h>
#include <fparameters/Dimension.h>
#include <fparameters/parameters.h>
#include "write.h"

void injectionNeutrons(Particle& n, Particle& p, State& st)
{
	show_message(msgStart,Module_neutronInjection);
	
	std::ofstream file;
	file.open("neutronInj.txt",std::ios::out);
	file 	<< "log(E/eV)"
			<< '\t' << "r/rS"
			<< '\t' << "pp"
			<< "\t" << "pY"
			<< std::endl;

	n.ps.iterate([&](const SpaceIterator& i) {
		double dens = st.denf_i.get(i);
		const double rB1 = i.val(DIM_R)/sqrt(paso_r);
		const double rB2 = rB1*paso_r;
		const double thetaH = st.thetaH.get(i);
		
		double tcross = i.val(DIM_R) / cLight; //corregir porque deberia ser (Rmax(theta)-r)/cLight;
		n.ps.iterate([&](const SpaceIterator& jE) {
			const double E = jE.val(DIM_E);
			double n_pg = neutronPgamma(E,tcross,n,p,st.photon.distribution,jE,
							st.photon.emin(),st.photon.emax());
			double n_pp = neutronPp(E,p,dens,i);
			double total = n_pg + n_pp;
			file << safeLog10(E/1.6e-12) << "\t"
				 << i.val(DIM_R)/schwRadius << "\t"
				 << safeLog10(n_pp) << "\t"
				 << safeLog10(n_pg) << std::endl;
			n.injection.set(jE,total); //en unidades de erg^-1 s^-1 cm^{-3}
		},{-1,i.coord[DIM_R],0});
	},{0,-1,0});
	file.close();
	
	show_message(msgEnd,Module_neutronInjection);
}