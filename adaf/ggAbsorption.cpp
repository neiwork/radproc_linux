#include "ggAbsorption.h"

#include "write.h"
#include "modelParameters.h"
#include "globalVariables.h"

#include <fmath/RungeKutta.h>
#include <fparameters/Dimension.h>
#include <fparameters/SpaceIterator.h>
#include <fmath/physics.h>



//E = Ega; E_ph = energia del photon target

double ggCrossSection(double x, double E_ph, double E)
{

	double Erest = electronMass*cLight2;

	double beta = sqrt( 1.0 - (2.0 * Erest * Erest) / ( (1.0 - x) * E_ph * E) );

	double sigma = (3.0*thomson/16.0)*(1.0-beta*beta)*( 2.0*beta*(beta*beta-2.0) + 
					(3.0-pow(beta,4.0))*log((1+beta)/(1-beta)) );
					
	return sigma;
}


double c_gg(double x, double E)  
{    
		
	double Erest = electronMass*cLight2;
	
	double result = 2.0*Erest*Erest/(E*(1.0-x));
	
	return result;
	
}


double f_gg(double x, double E_ph, double E, State& st, const SpaceCoord& distCoord)  
{    
	double Nph;
	if (E_ph < st.photon.emin() || E_ph > st.photon.emax()){
		Nph = 0.0;
	}
	else{
		Nph = st.photon.distribution.interpolate({ { 0, E_ph } }, &distCoord); 
	}
	
	double sigma = ggCrossSection(x, E_ph, E);
	
	double result = 2.0*pi*(1-x)*sigma*Nph;
	
	return result;
	
}


double internalAbs(int E_ix, State& st, double r_current)  
{
	
	double E = st.photon.ps[DIM_E][E_ix];  //E=Egamma

	double opacity = 0.0;


	st.photon.ps.iterate([&](const SpaceIterator &i) {
		
		double r = i.val(DIM_R);
		
		if (r >= r_current){ 
		
			//double thetaH = st.thetaH.get(i);
			double rB1=r/sqrt(paso_r);
			double rB2=r*sqrt(paso_r);
			double delta_r = rB2-rB1; //(r/sqrt(paso_r))*(paso_r-1.0);
			//double vol=(rB2*rB2*rB2-rB1*rB1*rB1)*(4.0/3.0)*pi*cos(thetaH);
		
			//double magneticField = st.magf.get(i);
				
			double phEmin = st.photon.emin();
			double phEmax = st.photon.emax();
				
//***************integral
			double xmin = -1.0;
			double xmax = 1.0;
			int n_x = 10;
			double xstep = (xmax-xmin)/n_x;
			
			int n_Eph = 30;
			
			double x = xmin;
			
			double L1 = 0;

			for (int i_x = 0; i_x < n_x; ++i_x)     //le saco el n para que se multiplique n veces y no n+1
			{
				double sup = phEmax;
				double inf = max(phEmin,c_gg(x,E));
				
				if(sup > inf){
					double Eph_int = pow((sup / inf), (1.0 / n_Eph));

					double E_ph = inf;

					for (int j_E = 0; j_E < n_Eph; ++j_E)
					{
						double dEph = E_ph*(Eph_int - 1);

						L1 += f_gg(x, E_ph, E, st, i)*dEph*xstep;

						E_ph = E_ph*Eph_int;
					}
					x = x + xstep;
				}
			}
//***********************************

			double integral = L1;
			/*RungeKutta(-1.0, 1.0,
				[E,phEmin](double x){
					return max(phEmin,c_gg(x,E));
				}, 
				[phEmax](double x){
					return phEmax;
				}, 
				[E,&st,i](double x, double E_ph){
					return f_gg(x, E_ph, E, st, i);  
					//return fICemi(u, t,E,creator, distCoord, tpf); 
				});  */
		
			opacity += integral*delta_r; 
		}
		
	},{E_ix,-1,0});
		
	return opacity;
}


void ggIntAbsorption(State& st, const std::string& filename)
{
	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);

	file << "Log(E/eV)" 
		<< "\t" << "r [M]"
		<< "\t" << "tau"
		<< "\t" << "exp(-tau)"
		<< std::endl;

	//double phEmin = st.photon.emin();
	//double phEmax = st.photon.emax();
	
	st.photon.ps.iterate([&](const SpaceIterator& i) {
		
		double E = i.val(DIM_E);
		double r = i.val(DIM_R);

		double tau = internalAbs(i.coord[DIM_R] , st, r); 

		file << log(E/1.6e-12) << "\t" << r/schwRadius
					 << "\t" << tau
					 << "\t" << exp(-tau)
					 << std::endl;
		
	},{-1,-1,0});
	
	file.close();
}