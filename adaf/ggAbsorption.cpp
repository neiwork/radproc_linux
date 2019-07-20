#include "ggAbsorption.h"

#include "write.h"
#include "modelParameters.h"
#include "globalVariables.h"
#include "adafFunctions.h"

#include <fmath/RungeKutta.h>
#include <fparameters/Dimension.h>
#include <fparameters/SpaceIterator.h>
#include <fmath/physics.h>



//E = Ega; E_ph = energia del photon target

double ggCrossSection(double x, double E_ph, double E)
{
	double Erest = electronMass*cLight2;
	double aux = 2.0*Erest*Erest / ((1.0-x)*E_ph*E);
	double beta = sqrt( abs(1.0 -  aux)); 

	return (beta < 1.0) ? (3.0*thomson/16.0)*(1.0-beta*beta)*( 2.0*beta*(beta*beta-2.0) + 
					(3.0-pow(beta,4.0))*log((1.0+beta)/(1.0-beta))) : 0.0;
}

double c_gg(double x, double E)  
{    
	double Erest = electronMass*cLight2;
	return 2.0*Erest*Erest/(E*(1.0-x));
}

double f_gg(double x, double E_ph, double E, State& st, const SpaceCoord& distCoord)  
{    
	double Nph = (E_ph > st.photon.emin() && E_ph < st.photon.emax()) ?
				st.photon.distribution.interpolate({{0,E_ph}},&distCoord) : 0.0;
	double sigma = ggCrossSection(x, E_ph, E);
	return (1-x)*sigma*Nph;
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
			opacity += integral*delta_r; 
		}
		
	},{E_ix,-1,0});
		
	return 2.0*pi*opacity;
}

double ggAbsortionCoeff(double E, double rlocal, State& st, const SpaceCoord& distCoord)  
{
	double xmin = -0.999;
	double xmax = 0.999;
	int n_x = 10;
	double dx = (xmax-xmin)/(n_x-1.0);
	int nEph = 20;
	double x = xmin;
	double integral2 = 0;
	for (int i_x=0;i_x<n_x;i_x++) {
		double inf = max(st.photon.emin(),c_gg(x,E));
		double sup = st.photon.emax();
		double integral1 = 0.0;
		if(sup > inf) {
			double Eph_int = pow(sup/inf,1.0/nEph);
			double Eph = inf;
			for (int jE=0;jE<nEph;++jE) {
				double dEph = (Eph/sqrt(Eph_int))*(Eph_int-1.0);
				double nPh = (rlocal < st.photon.ps[DIM_R][nR-1] && rlocal > st.photon.ps[DIM_R][0]) ?
								st.photon.distribution.interpolate({{0,Eph},{1,rlocal}},&distCoord) : 0.0;
				integral1 += nPh*ggCrossSection(x,Eph,E)*dEph;
				Eph *= Eph_int;
			}
		}
		integral2 += integral1*(1.0-x)*dx;
		x += dx;
	}
	return integral2;
}

double ggOpticalDepth(int E_ix, State& st, int iR)
{
	SpaceCoord psc = {E_ix,iR,0};
	
	double E = st.photon.ps[DIM_E][E_ix];
	double rMax = st.photon.ps[DIM_R][nR-1]*sqrt(paso_r);
	double rMin = st.photon.ps[DIM_R][0]/sqrt(paso_r);
	double pasoprim = pow(rMax/rMin,1.0/(nR*5.0));

	double r0 = st.photon.ps[DIM_R][iR];
	double drprim = r0*(pasoprim-1.0);
	double thetaprim = inclination*(pi/180.0);

	size_t nPhi = 5;
	double opticalDepth = 0.0;
	double theta0 = 0.5*pi;
	double phi0 = 0.0;
	double dPhi = 2.0*pi/nPhi;
	for (size_t jPhi=0;jPhi<nPhi;jPhi++) {

		double x0=r0*sin(theta0)*cos(phi0);
		double y0=r0*sin(theta0)*sin(phi0);
		double z0=r0*cos(theta0);
		
		double rprim = drprim;
		double r1 = r0;
		double theta1 = theta0;
		double thetaMinLocal = acos(costhetaH(r1));
		while (r1 < rMax && r1 > rMin && theta1 > thetaMinLocal && theta1 < pi-thetaMinLocal) {
			double xprim = rprim*sin(thetaprim);
			double zprim = rprim*cos(thetaprim);
			double x1 = x0+xprim;
			double z1 = z0+zprim;
			r1=sqrt(x1*x1+y0*y0+z1*z1);
			thetaMinLocal = acos(costhetaH(r1));
			
			theta1 = atan(sqrt(x1*x1+y0*y0)/abs(z1));

			opticalDepth += ggAbsortionCoeff(E,r1,st,psc)*drprim;
			
			drprim = r1*(pasoprim-1.0);
			rprim += drprim;
		}
		phi0 += dPhi;
	}
	return 2.0*pi*opticalDepth/nPhi;
}

void ggIntAbsorption(State& st, const std::string& filename)
{
	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);

	file << "Log(E/eV)" 
		<< "\t" << "r [M]"
		<< "\t" << "tau1"
		<< "\t" << "tau2"
		<< std::endl;

	for (int r_ix=0;r_ix< st.photon.ps[DIM_R].size();r_ix++) {
		double r = st.photon.ps[DIM_R][r_ix];
		st.photon.ps.iterate([&](const SpaceIterator& i) {
			double tau1 = internalAbs(i.coord[DIM_E] , st, r);
			double tau2 = ggOpticalDepth(i.coord[DIM_E],st,r_ix);

			file << log10(i.val(DIM_E) / 1.6e-12) << "\t" << r/schwRadius
						 << "\t" << tau1
						 << "\t" << tau2
						 << std::endl;
		},{-1,r_ix,0});
		file << std::endl;
	}
	file.close();
}