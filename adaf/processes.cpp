#include "processes.h"



#include "modelParameters.h"
#include "write.h"
#include "messages.h"

#include <fluminosities/luminositySynchrotron.h>
#include <fluminosities/luminosityIC.h>
#include <fluminosities/luminosityNTHadronic.h>
#include <fluminosities/luminosityPhotoHadronic.h>

#include <fparameters/Dimension.h>
#include <fparameters/SpaceIterator.h>

//#include <fparameters/parameters.h>

//#include <fmath/physics.h>
//#include <boost/property_tree/ptree.hpp>



/*double Llab(double Lint, double gamma)
{
	double Dlorentz = computeDlorentz(gamma);
	double boost = pow(Dlorentz, 4.0);
	return Lint*boost;
}*/


/* Takes [emi] =  E^2*[Q(E)] and calculates int(2.0*pi*P2(jetR)*emi dz); 
for [N(E)] = 1/erg, then it just sums over all z and returns erg/s  */




void processes(State& st, const std::string& filename)
{

	show_message(msgStart, Module_luminosities);

	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);
	
	double Emin = st.photon.emin();
	double Emax = st.photon.emax();
	

	file << "log(E/eV)"
			<< '\t' << "eSyn"
			<< '\t' << "pSyn"
			<< "\t" << "eIC"
			<< "\t" << "pPP"
			<< "\t" << "pPG"
			<< std::endl;

		
	int nE = st.photon.ps[DIM_E].size();
	
	for (int E_ix = 0; E_ix < nE; E_ix++) {


		double E = st.photon.distribution.ps[DIM_E][E_ix];
		double fmtE = log10(E / 1.6e-12);
		
		double eSyn = 0.0;
		double eIC = 0.0;
		double pSyn = 0.0;
		double pPP = 0.0;
		double pPG = 0.0;

		st.photon.ps.iterate([&](const SpaceIterator &i) {

			//const double E = i.val(DIM_E);
			
			double density = st.denf_e.get(i)+st.denf_i.get(i);
			
			
			eSyn += luminositySynchrotron(E, st.electron, i, st.magf);
			eIC  += luminosityIC(E,st.electron,i.coord,st.photon.distribution,Emin); //ver unidades del distribution XXX
			
			//luminosityIC(E, st.electron, i.coord, [&st, &E, &r](double E) {
						//	return st.photon.distribution.interpolate({ {DIM_E, E },{ DIM_R, r } });},Emin);
							
			
			pSyn += luminositySynchrotron(E, st.proton, i, st.magf);
				
			pPP  += luminosityNTHadronic(E, st.proton, density, i);
			
			pPG  += luminosityPhotoHadronic(E, st.proton, st.photon.distribution, i, Emin, Emax);				
							
		
		}, { E_ix, -1 , 0});
		
			file << fmtE
				<< '\t' << safeLog10(eSyn)
				<< '\t' << safeLog10(pSyn)
				<< "\t" << safeLog10(eIC)
				<< "\t" << safeLog10(pPP)
				<< "\t" << safeLog10(pPG)
				<< std::endl;

	}  
}
	
	





///////////////////
//#   pragma omp parallel for \
	//		private(i, eSyn, eIC) \
	//		shared(st, Qsyn, Qic) \
	//		default(none) \
	//		schedule(static, 1) \
	//		num_threads(2)

//#pragma omp parallel sections
//{
//#pragma omp section
//	{
//		Qsyn.fill([&st](const SpaceIterator &i){
//			return luminositySynchrotron(i.val(DIM_E), st.electron); //estos devuelven erg/s/cm^3, integrar!
//		});
//	}

//#pragma omp section
//	{
//		Qic.fill([&st](const SpaceIterator &i){
//			return luminosityAnisotropicIC(i.val(DIM_E), st.electron, i.val(DIM_R));
//		});
//	}
//}



//double dz = z[i]*(z_int - 1);
//volumen de la celda i
//double vol_i = pi*P2(jetRadius(z[i], openingAngle))*dz;;
//double E = pps[0][E_ix];