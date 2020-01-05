#include "jetEmission.h"

// Standard C++ libraries
#include <iomanip>

// Boost libraries
#include <boost/property_tree/ptree.hpp>
#include <boost/math/tools/roots.hpp>

// Standard user-made libraries
#include <fparameters/SpaceIterator.h>
#include <fparameters/parameters.h>
#include <fmath/physics.h>
#include <fparameters/Dimension.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_gamma.h>
#include "absorption.h"
// Project headers
#include "globalVariables.h"
#include "messages.h"
#include "read.h"
#include "adafFunctions.h"
#include "write.h"
#include <fmath/RungeKutta.h>

#include <fluminosities/thermalSync.h>
#include <fluminosities/blackBody.h>
#include "absorption.h"


double n_pl(double g, double gammaC, double gammaMin, double Npl, double pJet)
{
	if (gammaMin > gammaC) {
		double factor = 0.0;
		if (g > gammaC && g < gammaMin)
			factor = pow(gammaMin,1.0-pJet)/(g*g);
		else if (g > gammaMin)
			factor = pow(g,-(pJet+1.0));
		return Npl*(pJet-1.0)*gammaC*factor;
	} else {
		double factor = 0.0;
		if (g >= gammaMin && g <= gammaC)
			factor = 1.0;
		else if (g > gammaC)
			factor = gammaC/g;
		return Npl*(pJet-1.0)*factor*pow(g,-pJet);
	}
}

double n_pl_2(double g, double Npl, double pJet)
{
	return Npl*(pJet-1.0)*pow(g,-pJet);
}

double jSyJet(double g, double nu, double magf)
{
	double U_B = magf*magf/(8.0*pi);
	double nuB = electronCharge*magf/(2*pi*electronMass*cLight);
	double t = nu / (3.0*g*g*nuB);
	double K_43,K_13;
	if (t > 30) {
		K_43 = K_13 = exp(-t) * sqrt(0.5*pi/t);
	} else if (t < 0.01*sqrt(2)) {
		K_13 = 0.5*gsl_sf_gamma(1.0/3.0) * pow(2.0/t,1.0/3.0);
		K_43 = 0.5*gsl_sf_gamma(4.0/3.0) * pow(2.0/t,4.0/3.0);
	} else {
		K_43 = gsl_sf_bessel_Knu(4.0/3.0,t);
		K_13 = gsl_sf_bessel_Knu(1.0/3.0,t);
	}
	return 3.0*sqrt(3.0)/pi * thomson * cLight * U_B * t*t * (K_43*K_13-3*t/5 * (K_43*K_43-K_13*K_13));
}

double jSyJet2(double g, double nu, double magf)
{
	double nuC = 3.0/(4.0*pi)*electronCharge*magf/(electronMass*cLight) * g*g;
	double x = nu/nuC;
	return sqrt(3.0)*P3(electronCharge)/electronRestEnergy * magf *1.85*pow(x,1.0/3.0)*exp(-x);
}

double sigma_Sy(double g, double nu, double magf)
{
	double nuB = electronCharge*magf/(2*pi*electronMass*cLight);
	double t = nu/(3*g*g*nuB);
	double K_43,K_13;
	if (t > 30) {
		K_43 = K_13 = exp(-t) * sqrt(0.5*pi/t);
	} else if (t < 0.01*sqrt(2)) {
		K_13 = 0.5*gsl_sf_gamma(1.0/3.0) * pow(2.0/t,1.0/3.0);
		K_43 = 0.5*gsl_sf_gamma(4.0/3.0) * pow(2.0/t,4.0/3.0);
	} else {
		K_43 = gsl_sf_bessel_Knu(4.0/3.0,t);
		K_13 = gsl_sf_bessel_Knu(1.0/3.0,t);
	}
	return 8.0*sqrt(3.0)/15 *pi*pi*electronCharge / magf * t / pow(g,5.0) * (K_43*K_43-K_13*K_13);
}

double tau_Sy(double deltaR, double nu, double gammaC, double gammaMin, double Npl, double pJet, double magf)
{
	return deltaR * RungeKuttaSimple(gammaMin,1.0e10,[&nu,&gammaC,&gammaMin,&Npl,&pJet,&magf]
						(double g) {return sigma_Sy(g,nu,magf)*n_pl(g,gammaC,gammaMin,Npl,pJet);});
}

void jetProcesses(State& st, const string& filename)
{
	double accRateFraction = GlobalConfig.get<double>("nonThermal.jet.accRateFraction");
	double accRate = accRateFraction*accRateADAF(schwRadius);
	double openingAngle = GlobalConfig.get<double>("nonThermal.jet.openingAngle");
	double bulkLorentzFactor = GlobalConfig.get<double>("nonThermal.jet.lorentzFactor");
	double inclinationAngle = GlobalConfig.get<double>("nonThermal.jet.inclination")*(pi/180.0);
	double nt_electronFraction = GlobalConfig.get<double>("nonThermal.jet.etaInj");
	double eB = GlobalConfig.get<double>("nonThermal.jet.magneticEnergyFraction");
	double eElectrons = GlobalConfig.get<double>("nonThermal.jet.electronEnergyFraction");
	double pJet = GlobalConfig.get<double>("nonThermal.jet.pIndex");
	double zMin = GlobalConfig.get<double>("nonThermal.jet.zMin");
	
	double z0 = P2(bulkLorentzFactor)*schwRadius*zMin;
	double zMax = z0*1.0e4;
	size_t nZ = 100;
	
	double gamma2 = sqrt(0.5*(P2(bulkLorentzFactor)+1.0));
	double betaJet = sqrt(1.0-1.0/P2(bulkLorentzFactor));
	double vJet = betaJet*cLight;
	double n1 = accRate/(pi*z0*z0*openingAngle*openingAngle*vJet)/protonMass;
	double n2 = (4.0*gamma2+3.0)*n1;
	double Ti = 6.3e11;
	double Te = 1.0e9;
	
	double pasoZ = pow(zMax/z0,1.0/nZ);
	double z = z0;
	size_t jZ = 0;
	Matrix lumBeamed;
	matrixInit(lumBeamed,nE,nZ,0.0);
	while (z < zMax) {
		z *= pasoZ;
		
		double dz = z*(pasoZ-1.0);
		double area = pi*P2(openingAngle*z);
		double vol = area*dz;
		double dens = n2 * P2(z0/z);
		double magf = sqrt(8.0*pi*eB*(gamma2-1.0)*dens*protonMass*cLight2);
		
		// NONTHERMAL POPULATION
		double gammaMin = max((gamma2-1.0)*(pJet-2.0)/(pJet-1.0) * (protonMass/electronMass) * 
								(eElectrons/nt_electronFraction),1.0);
		gammaMin = 2.0;
		double tDyn = z/cLight;
		double Npl = nt_electronFraction*dens*pow(gammaMin,pJet-1.0);
		double gammaC = 6*pi*electronRestEnergy/(magf*magf*thomson*z);
		
		size_t jE = 0;
		st.photon.ps.iterate([&](const SpaceIterator& i) {
			double nu = i.val(DIM_E)/planck;
			double beta = sqrt(1.0-1.0/P2(bulkLorentzFactor));
			double doppler = 1.0/(bulkLorentzFactor*(1.0-beta*cos(inclinationAngle)));
			double nuDoppler = nu/doppler;
			
			double lumLocal = 0.0;
			double gamma1 = 2.0;
			double gammaMax = 1.0e7;
			size_t nnE = 1000;
			double pasoG = pow(gammaMax/gamma1,1.0/nnE);
			double gamma = gamma1;
			for (size_t jjE=0;jjE<nnE;jjE++) {
				double dgamma = gamma*(pasoG-1.0);
				//if (jE==0)
				//	cout << z/schwRadius << "\t" << gamma << "\t" << n_pl(gamma,gammaC,gammaMin,Npl,pJet)*dgamma << endl;
				lumLocal += dgamma*n_pl(gamma,gammaC,gammaMin,Npl,pJet)*jSyJet2(gamma,nuDoppler,magf);
				//lumLocal += dgamma*n_pl_2(gamma,Npl,pJet)*jSyJet2(gamma,nuDoppler,magf);
				gamma *= pasoG;
			}
			lumLocal *= vol;
			
			//double lumLocal = vol*RungeKuttaSimple(gammaC,gammaC*1.0e5,[&gammaC,&gammaMin,&Npl,
			//					&pJet,&magf,&nuDoppler](double g){return n_pl(g,gammaC,gammaMin,Npl,pJet)*
			//					jSyJet2(g,nuDoppler,magf);});
			
			double a_th = jSync(planck*nu,Te,magf,dens)/bb(nu,Te);
			double a_pl = RungeKuttaSimple(gammaMin,gammaMin*1.0e10,[&](double g)
						{return n_pl(g,gammaC,gammaMin,Npl,pJet)*sigma_Sy(g,nu,magf);});
			double tau = z*openingAngle*(a_pl+a_th);
			//tau = tau_Sy(z*openingAngle,nu,gammaC,gammaMin,Npl,pJet,magf);
			lumLocal *= (tau > 1.0e-10 ? (1.0-exp(-tau))/tau : 1.0);
			lumBeamed[jE][jZ] = lumLocal * doppler*doppler*doppler;
			jE++;
		},{-1,0,0});
		jZ++;
	}
	
	ofstream file;
	file.open(filename.c_str(),ios::out);
	size_t jE = 0.0;
	st.photon.ps.iterate([&](const SpaceIterator& i) {
		double nu = i.val(DIM_E)/planck;
		double lum = 0.0;
		for (size_t jZ=0;jZ<nZ;jZ++) {
			lum += lumBeamed[jE][jZ];
		}
		file << nu << "\t" << nu*lum << endl;
		jE++;
	},{-1,0,0});
	file.close();
}


void jetProcesses2(State& st, const string& filename)
{
	double accRateFraction = GlobalConfig.get<double>("nonThermal.jet.accRateFraction");
	double accRate = accRateFraction*accRateADAF(schwRadius);
	double openingAngle = GlobalConfig.get<double>("nonThermal.jet.openingAngle");
	double bulkLorentzFactor = GlobalConfig.get<double>("nonThermal.jet.lorentzFactor");
	double inclinationAngle = GlobalConfig.get<double>("nonThermal.jet.inclination")*(pi/180.0);
	double nt_electronFraction = GlobalConfig.get<double>("nonThermal.jet.etaInj");
	double pJet = GlobalConfig.get<double>("nonThermal.jet.pIndex");
	double mB = 1.5;
	size_t nZ = 100;
	
	double z0 = 50*schwRadius;
	double r0 = openingAngle*z0;
	double betaJet = sqrt(1.0-1.0/P2(bulkLorentzFactor));
	double vJet = betaJet*cLight;
	double B0 = sqrt(accRate*cLight2*8/(r0*r0*vJet));
	double n0 = accRate/(bulkLorentzFactor*pi*r0*r0*vJet*protonMass);
	double zMax = z0*1.0e5;
	
	double pasoZ = pow(zMax/z0,1.0/nZ);
	double z = z0;
	size_t jZ = 0;
	Matrix lumBeamed;
	matrixInit(lumBeamed,nE,nZ,0.0);
	while (z < zMax) {
		z *= pasoZ;
		double dz = z*(pasoZ-1.0);
		double area = pi*P2(openingAngle*z);
		double vol = area*dz;
		double dens = n0 * P2(z0/z);
		double magf = B0 * pow(z0/z,mB);
		
		// NONTHERMAL POPULATION
		double gammaMin = 2.0;
		double Npl = nt_electronFraction*dens*pow(gammaMin,pJet-1.0);
		double gammaC = 6*pi*electronRestEnergy/(magf*magf*thomson*z);
		
		size_t jE = 0;
		st.photon.ps.iterate([&](const SpaceIterator& i) {
			double nu = i.val(DIM_E)/planck;
			double beta = sqrt(1.0-1.0/P2(bulkLorentzFactor));
			double doppler = 1.0/(bulkLorentzFactor*(1.0-beta*cos(inclinationAngle)));
			double nuDoppler = nu/doppler;
			double lumLocal = vol*RungeKuttaSimple(gammaMin,gammaMin*1.0e10,[&gammaC,&gammaMin,&Npl,
									&pJet,&magf,&nuDoppler](double g){return n_pl(g,gammaC,gammaMin,Npl,pJet)*
										jSyJet2(g,nuDoppler,magf);});
			double a_pl = RungeKuttaSimple(gammaMin,gammaMin*1.0e10,[&](double g)
							{return n_pl(g,gammaC,gammaMin,Npl,pJet)*sigma_Sy(g,nu,magf);});
			double tau = z*openingAngle*a_pl;
			//tau = tau_Sy(z*openingAngle,nu,gammaC,gammaMin,Npl,pJet,magf);
			lumLocal *= (tau > 1.0e-10 ? (1.0-exp(-tau))/tau : 1.0);
			lumBeamed[jE][jZ] = lumLocal * doppler*doppler*doppler;
			jE++;
		},{-1,0,0});
		jZ++;
	}
	
	ofstream file;
	file.open(filename.c_str(),ios::out);
	size_t jE = 0.0;
	st.photon.ps.iterate([&](const SpaceIterator& i) {
		double nu = i.val(DIM_E)/planck;
		double lum = 0.0;
		for (size_t jZ=0;jZ<nZ;jZ++) {
			lum += lumBeamed[jE][jZ];
		}
		file << nu << "\t" << nu*lum << endl;
		jE++;
	},{-1,0,0});
	file.close();
}