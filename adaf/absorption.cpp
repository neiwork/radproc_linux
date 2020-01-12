#include "absorption.h"

#include "write.h"
#include "modelParameters.h"
#include "globalVariables.h"
#include "adafFunctions.h"

#include <fmath/RungeKutta.h>
#include <fparameters/Dimension.h>
#include <fparameters/SpaceIterator.h>
#include <fluminosities/thermalSync.h>
#include <fluminosities/luminositySynchrotron.h>
#include <fluminosities/blackBody.h>
#include <fmath/physics.h>
#include <fluminosities/opticalDepthSSA.h>

double ggCrossSection(double x, double E_ph, double E)
{
	double Erest = electronMass*cLight2;
	double aux = 2.0*Erest*Erest / ((1.0-x)*E_ph*E);
	double beta = sqrt( abs(1.0 -  aux)); 
	return (beta < 1.0) ? (3.0*thomson/16.0)*(1.0-beta*beta)*( 2.0*beta*(beta*beta-2.0) + 
					(3.0-pow(beta,4.0))*log((1.0+beta)/(1.0-beta))) : 0.0;
}

double ggCrossSection2(double E, double Eph)
{
	double s0 = E*Eph/(electronRestEnergy*electronRestEnergy);
	double term1 = (s0+0.5*log(s0)-1.0/6.0+0.5/s0)*log(sqrt(s0)+sqrt(s0-1));
	double term2 = (s0+4.0/9.0-1.0/(9.0*s0))*sqrt(1.0-1.0/s0);
	return (s0 > 1.0) ? 1.5*thomson/(s0*s0) * (term1-term2) : 0.0;
}

double c_gg(double x, double E)  
{    
	double Erest = electronMass*cLight2;
	return 2.0*Erest*Erest/(E*(1.0-x));
}

double f_gg(double x, double E_ph, double E, State& st, const SpaceCoord& distCoord)  
{    
	double Nph = (E_ph > st.photon.emin() && E_ph < st.photon.emax()) ?
				st.photon.distribution.interpolate({{DIM_E,E_ph}},&distCoord) : 0.0;
	double sigma = ggCrossSection(x, E_ph, E);
	return (1-x)*sigma*Nph;
}

/*
double internalAbs(int E_ix, State& st, double r_current)  
{
	double E = st.photon.ps[DIM_E][E_ix];
	double opacity = 0.0;

	st.photon.ps.iterate([&](const SpaceIterator &i) {		
		double r = i.val(DIM_R);
		if (r >= r_current){ 
			//double thetaH = st.thetaH.get(i);
			double rB1=r/sqrt(paso_r);
			double rB2=r*sqrt(paso_r);
			double delta_r = rB2-rB1;

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
}*/

double ggAbsorptionCoeff2(double E, double rlocal, State& st, const SpaceCoord& distCoord)
{
	return (rlocal < st.photon.ps[DIM_R][nR-1] && rlocal > st.photon.ps[DIM_R][0]) ?
			RungeKuttaSimple(st.photon.emin(),st.photon.emax(),[&E,&rlocal,&st,&distCoord](double Eph)
				{return st.photon.distribution.interpolate({{DIM_E,Eph},{DIM_R,rlocal}},&distCoord)*
						ggCrossSection2(Eph,E);}) : 0.0;
}

double ggAbsorptionCoeff(double E, double rlocal, State& st, const SpaceCoord& distCoord)  
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
								st.photon.distribution.interpolate({{DIM_E,Eph},{DIM_R,rlocal}},&distCoord) : 0.0;
				integral1 += nPh*ggCrossSection(x,Eph,E)*dEph;
				Eph *= Eph_int;
			}
		}
		integral2 += integral1*(1.0-x)*dx;
		x += dx;
	}
	return integral2;
}

double ssaAbsorptionCoeff(double energy, double magf, Particle& p, SpaceCoord& psc)
{
	double integral = integSimpson(log(p.emin()),log(p.emax()),[energy,&p,magf,&psc](double logx)
	{
		double x = exp(logx);
		return fSSA2(x,energy,p,magf,psc)*x;
	},50);
	return -planck*planck*planck*cLight2*integral/(8*pi*energy*energy);
}

void opticalDepth(Vector& tau, int E_ix, State& st, int iR)
{
	SpaceCoord psc = {E_ix,iR,0};
	
	double E = st.photon.ps[DIM_E][E_ix];
	double rMax = st.photon.ps[DIM_R][nR-1];
	double rMin = st.photon.ps[DIM_R][0];
	double pasoprim = pow(rMax/rMin,1.0/(nR*10.0));

	double r0 = st.photon.ps[DIM_R][iR];
	double drprim = r0*(pasoprim-1.0);
	double thetaprim = inclination*(pi/180.0);

	size_t nPhi = 5;
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
		do {
			double xprim = rprim*sin(thetaprim);
			double zprim = rprim*cos(thetaprim);
			double x1 = x0+xprim;
			double z1 = z0+zprim;
			r1=sqrt(x1*x1+y0*y0+z1*z1);
			thetaMinLocal = acos(costhetaH(r1));
			theta1 = atan(sqrt(x1*x1+y0*y0)/abs(z1));
			double magneticField = ((r1 < rMax && r1 > rMin) && 
							(theta1 > thetaMinLocal && theta1 < pi-thetaMinLocal)) ? 
							st.magf.interpolate({{DIM_R,r1}},&psc) : 0.0;
			double temp = (r1 < rMax && r1 > rMin) && (theta1 > thetaMinLocal && 
							theta1 < pi-thetaMinLocal) ? electronTemp(r1) : 0.0;
			double dens_e = (r1 < rMax && r1 > rMin) && (theta1 > thetaMinLocal && 
							theta1 < pi-thetaMinLocal) ? electronDensity(r1) : 0.0;
			double ssaAbsCoeffe = 0.0;
			double ssaAbsCoeffp = 0.0;
			double ssaAbsCoeffeTh = 0.0;
			double ggAbsCoeff = 0.0;
			double fmtE = log10(E/1.6e-12);
			
			double frequency = E/planck;
			if (fmtE < 0.0) {
				ssaAbsCoeffe = ssaAbsorptionCoeff(E,magneticField,st.ntElectron,psc);
				double bbody = bb(frequency,temp);
				double jSy = jSync(E,temp,magneticField,dens_e);
				//double jSy = luminositySynchrotron2(E,st.ntElectron,psc,magneticField)/frequency/(4.0*pi);
				ssaAbsCoeffeTh = (bbody > 0.0) ? jSy/bbody : 0.0;
				ssaAbsCoeffp = ssaAbsorptionCoeff(E,magneticField,st.ntProton,psc);
			}
			//cout << "bb = " << bb_RJ(frequency,temp) << "\t ssaAbsCoeffeTh = " << ssaAbsCoeffeTh << endl;
			if (fmtE > 7.0) ggAbsCoeff = ggAbsorptionCoeff2(E,r1,st,psc);
			
			tau[0] += ggAbsCoeff*drprim;
			tau[1] += (ssaAbsCoeffe+ssaAbsCoeffeTh)*drprim;
			//tau[1] += ssaAbsCoeffe*drprim;
			tau[2] += ssaAbsCoeffp*drprim;
			
			drprim = r1*(pasoprim-1.0);
			rprim += drprim;
		} while (r1 < rMax && r1 > rMin && theta1 > thetaMinLocal && theta1 < pi-thetaMinLocal);
		phi0 += dPhi;
	}
	tau[0] /= nPhi;
	tau[1] /= nPhi;
	tau[2] /= nPhi;
}