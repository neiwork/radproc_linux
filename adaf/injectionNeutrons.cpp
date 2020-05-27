#include "injectionNeutrons.h"

#include "messages.h"
#include "globalVariables.h"
#include "adafFunctions.h"
//#include <fmath/fbisection.h>
#include <fmath/RungeKutta.h>
//#include <fmath/physics.h>
//#include <boost/property_tree/ptree.hpp>
#include <iostream>
#include <finjection/pgammaPionInj.h>
#include <flosses/lossesPhotoHadronic.h>
#include <finjection/neutronPp.h>
#include <finjection/neutronPgamma.h>
#include <fparameters/SpaceIterator.h>
#include <fparameters/Dimension.h>
#include <fparameters/parameters.h>
#include <fparticle/Particle.h>
#include "write.h"

void injectionNeutrons(State& st)
{
	show_message(msgStart,Module_neutronInjection);
	
    Particle &p = st.ntProton;
    Particle &n = st.ntNeutron;
    
	std::ofstream file;
	file.open("neutronInjection.dat",std::ios::out);
	file 	<< "log(E/eV)"
			<< '\t' << "r/rS"
			<< '\t' << "pp"
			<< "\t" << "pY"
			<< std::endl;

	double sumQ = 0.0;
	n.ps.iterate([&](const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		double dens = st.denf_i.get(iR);
		const double vol = volume(r);
		double height = height_fun(r);
		double tcross = height / cLight; //corregir porque deberia ser (Rmax(theta)-r)/cLight;
		n.ps.iterate([&](const SpaceIterator& iRE) {
			const double En = iRE.val(DIM_E);
			//double n_pg = neutronPgamma(En,tcross,n,p,st.photon.distribution,iR,
			//				st.photon.emin(),st.photon.emax());
			/*double Ep = 4.0/3.0 * En;
			double t_pg = t_pion_PHsimple(Ep,st.ntProton,st.photon.distribution,iR, 
					st.photon.emin(),st.photon.emax());
			t_pg = Ep/lossesPhotoMeson(Ep,st.ntProton,st.photon.distribution,iR,st.photon.emin(),
									st.photon.emax());
			double n_pg = 4.0/3.0 * st.ntProton.distribution.interpolate({{DIM_E,Ep}},&iR.coord) / t_pg;*/
			double n_pg = neutronPgamma(En,n,p,st.photon.distribution,iR,st.photon.emin(),st.photon.emax());
			double n_pp = neutronPp(En,p,dens,iR);
			double total = n_pp + n_pg;
			file << safeLog10(En/EV_TO_ERG) << "\t"
				 << iR.val(DIM_R)/schwRadius << "\t"
				 << safeLog10(n_pp*vol)  << "\t"
				 << safeLog10(n_pg*vol) << endl;
			n.injection.set(iRE,total); //en unidades de erg^-1 s^-1 cm^{-3}
		},{-1,iR.coord[DIM_R],0});
		sumQ += integSimpsonLog(n.emin(), n.emax(), [&] (double En)
				{
					return En * n.injection.interpolate({{DIM_E,En}},&iR.coord);
				},30) * vol;
	},{0,-1,0});
	cout << "Total power injected in " << n.id << " = " << sumQ << endl;

	file.close();
	show_message(msgEnd,Module_neutronInjection);
}