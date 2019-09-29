#include "distributionNeutrons.h"
#include <fparameters/parameters.h>
#include "globalVariables.h"
#include <fparameters/Dimension.h>
#include <fparameters/parameters.h>
#include <fparameters/SpaceIterator.h>
#include <finjection/neutronDecay.h>
#include "write.h"
#include "adafFunctions.h"

//extern "C" {
//	#include <nrMath/random.h>
//}

using namespace std;

void distributionNeutrons(State& st)
{
	ofstream fileNeutrons,fileProtons,fileElectrons;
	fileNeutrons.open("neutronProp.txt",ios::out);
	fileProtons.open("protonInjByNeutrons.txt",ios::out);
	fileElectrons.open("electronInjByNeutrons.txt",ios::out);
	Vector Q0n(nE,0.0);
	size_t jE = 0;
	st.ntNeutron.ps.iterate([&](const SpaceIterator &iE) {
		st.ntNeutron.ps.iterate([&](const SpaceIterator &iER) {
			double rB2 = iER.val(DIM_R)*sqrt(paso_r);
			double rB1 = rB2/paso_r;
			double vol = (4.0/3.0)*pi*cos(st.thetaH.get(iER))*(rB2*rB2*rB2-rB1*rB1*rB1);
			Q0n[jE] += st.ntNeutron.injection.get(iER) * vol;
		},{iE.coord[DIM_E],-1,0});
		jE++;
	},{-1,0,0});

	size_t nRn = nE;
	double r = cLight*neutronMeanLife*1.0e-4;
	double pasoRn = pow(st.ntNeutron.emax()*1.0e4/(neutronMass*cLight2),1.0/nRn);
	
	double Qp_tot = 0.0;
	double Qe_tot = 0.0;
	for (size_t jR=0;jR<nRn;jR++) {
		r *= pasoRn;
		double NnFactor = 1.0 / (4.0 * pi * cLight) / (r*r);
		double vol = 4.0*pi*r*r*r*(pasoRn-1.0);
		size_t jE = 0;
		st.ntNeutron.ps.iterate([&](const SpaceIterator &iE) {
			double E = iE.val(DIM_E);
			double gamma = E / (neutronMass*cLight2);
			double tDecay = neutronMeanLife*gamma;
			double rDec = cLight*tDecay;
			double Nn = NnFactor * Q0n[jE]*exp(-r/rDec);
			st.ntNeutron.distribution.set(iE,Nn);
			jE++;
			fileNeutrons << safeLog10(r) << "\t" << safeLog10(E/1.6e-12) << "\t"
						 << safeLog10(Nn) << "\t"
						 << safeLog10(Nn*E*E*4.0*pi*r*r) << endl;
		},{-1,0,0});
		
		double Qp_local = 0.0;
		double pasoEp = pow(st.ntProton.emax()/st.ntProton.emin(),1.0/nE);
		st.ntProton.ps.iterate([&](const SpaceIterator &iE) {
			double E = iE.val(DIM_E);
			double tDecay = neutronMeanLife*(E/(neutronMass*cLight2));
			double Qp = st.ntNeutron.distribution.get(iE)/tDecay;
			Qp_local += Qp*E*E*(pasoEp-1.0);
			fileProtons << safeLog10(r) << "\t" << safeLog10(E/1.6e-12) << "\t"
						<< safeLog10(Qp) << "\t"
						<< safeLog10(Qp*E*E*4.0*pi*r*r*r) << endl;
		},{-1,0,0});
		
		double Qe_local = 0.0;
		double pasoEe = pow(st.ntElectron.ps[DIM_E][nE-1]/st.ntElectron.ps[DIM_E][0],1.0/nE);
		st.ntElectron.ps.iterate([&](const SpaceIterator &iE) {
			double E = iE.val(DIM_E);
			double Qe = injElectronNeutronDecay(E,st.ntNeutron,iE);
			Qe_local += Qe*E*E*(pasoEe-1.0);
			fileElectrons << safeLog10(r) << "\t" << safeLog10(E/1.6e-12) << "\t"
						  << safeLog10(Qe) << "\t"
						  << safeLog10(Qe*E*E*4.0*pi*r*r*r) << endl;
		},{-1,0,0});
		
		Qp_tot += Qp_local*vol;
		Qe_tot += Qe_local*vol;
	}
	
	cout << "Total power injected in " << st.ntProton.id << " by neutron decay = " << Qp_tot << endl;
	cout << "Total power injected in " << st.ntElectron.id << " by neutron decay = " << Qe_tot << endl;
	
	fileNeutrons.close();
	fileProtons.close();
	fileElectrons.close();
}

#include <nrMath/ran1.h>

void jetNeutronDecay(State& st)
{
	double openAngle = GlobalConfig.get<double>("nonThermal.neutrons.jetDecay.openingAngle")*(pi/180.0);
	double m2 = P2(tan(openAngle));
	double zMin = GlobalConfig.get<double>("nonThermal.neutrons.jetDecay.zMin");
	double zMax = GlobalConfig.get<double>("nonThermal.neutrons.jetDecay.zMax");
	
	size_t nNeutrons = GlobalConfig.get<size_t>("nonThermal.neutrons.jetDecay.nNeutrons");
	size_t nTheta = GlobalConfig.get<size_t>("nonThermal.neutrons.jetDecay.nTheta");
	Vector Qjet(nE,0.0);

	long *idum;
	long j = -5;
	idum = &j; 
	
	ofstream fileRandoms;
	fileRandoms.open("randoms.txt",ios::out);
	
	st.ntNeutron.ps.iterate([&](const SpaceIterator& iR) {
		double rB2 = iR.val(DIM_R)*sqrt(paso_r);
		double rB1 = rB2/paso_r;
		double r0 = iR.val(DIM_R)/schwRadius;
		double thetaMin = st.thetaH.get(iR);
		double vol = (4.0/3.0)*pi*cos(thetaMin)*(rB2*rB2*rB2-rB1*rB1*rB1);
		double dyaux = (1.0-sin(thetaMin))/nTheta;
		Vector Pdecay(nE,0.0);
		size_t count = 0;
        for(size_t kTh=1;kTh<=nTheta;kTh++) {
            double yaux = sin(thetaMin)+kTh*dyaux;  // Theta distributed uniformly in sin(theta).
            double theta0 = asin(yaux);
			double y0 = r0*sin(theta0);
            double z0 = r0*cos(theta0);
			//y0 = r0;
			//z0 = 0.0;
			for(size_t jN=1;jN<=nNeutrons;jN++) {
                //double random_number = gsl_rng_uniform(RandomNumberGenerator);
				double random_number = ran1(idum);
                double phiprim = 2.0*pi*random_number;
                //random_number = gsl_rng_uniform(RandomNumberGenerator);
				random_number = ran1(idum);
                double muprim = 1.0-2.0*random_number;
				double a = 1.0-(1.0+m2)*P2(muprim);
				double b = 2.0*(y0*sqrt(1-P2(muprim))*sin(phiprim)-z0*muprim*m2);
				double c = y0*y0 - z0*z0*m2;
				double b_2a = 0.5*b/a;
				double c_a = c/a;
				double disc = b_2a*b_2a - c_a;
				if (disc > 0.0) {
					int indicator = 0;
					double r1,r2;
					r1 = r2 = 0.0;
					double root1 = -b_2a - sqrt(disc);
					double root2 = -b_2a + sqrt(disc);
					double z_root1 = root1*muprim+z0;
					double z_root2 = root2*muprim+z0;
					double abszr1 = abs(z_root1);
					double abszr2 = abs(z_root2);
					if (root1 > 0.0 && root2 > 0.0) {
						double r_zMin = (z_root1 > 0.0) ? (zMin-z0)/muprim : (-zMin-z0)/muprim;
						double r_zMax = (z_root1 > 0.0) ? (zMax-z0)/muprim : (-zMax-z0)/muprim;
						if (min(abszr1,abszr2) > zMin && max(abszr1,abszr2) < zMax) {
							r1 = root1;
							r2 = root2;
							indicator = 1;
						} else if (abszr1 < zMin && abszr2 > zMin && abszr2 < zMax) {
							r1 = r_zMin;
							r2 = root2;
							indicator = 2;
						} else if (z_root1 > zMax && z_root2 > zMin && z_root2 < zMax) {
							r1 = r_zMax;
							r2 = root2;
							indicator = 3;
						} else if (abszr1 > zMin && abszr1 < zMax && abszr2 > zMax) {
							r1 = root1;
							r2 = r_zMax;
							indicator = 4;
						} else if (z_root1 > zMin && z_root1 < zMax && z_root2 < zMin) {
							r1 = root1;
							r2 = r_zMin;
							indicator = 5;
						} else if (abszr1 < zMin && abszr2 > zMax) {
							r1 = r_zMin;
							r2 = r_zMax;
						} else if (z_root1 > zMax && z_root2 < zMin) {
							r1 = r_zMax;
							r2 = r_zMin;
							indicator = 6;
						} else if (min(abszr1,abszr2) > zMax || max(abszr1,abszr2) < zMin) {
							r1 = 0.0;
							r2 = 0.0;
							indicator = 0;
						}
					} else if (root1*root2 < 0.0 && abszr2 < zMax) {
						double r_zMin = (muprim > 0.0) ? (zMin-z0)/muprim : (-zMin-z0)/muprim;
						double r_zMax = (muprim > 0.0) ? (zMax-z0)/muprim : (-zMax-z0)/muprim;
						r1 = max(r_zMin,root2);
						r2 = r_zMax;
						indicator = 7;
					}
					
					if (indicator) {
						count++;
						size_t jE = 0;
						st.ntNeutron.ps.iterate([&](const SpaceIterator& iRE) {
							double gamma = iRE.val(DIM_E)/(neutronMass*cLight2);
							double tDecay = neutronMeanLife*gamma;
							double rDecay = cLight*tDecay/schwRadius;
							double Plocal = exp(-r1/rDecay) * (1.0-exp(-(r2-r1)/rDecay));
							Pdecay[jE] += Plocal;
							jE++;
						},{-1,iR.coord[DIM_R],0});
					}
				}
            }
        }
		size_t jE = 0;
		st.ntNeutron.ps.iterate([&](const SpaceIterator& iRE) {
			double fE = Pdecay[jE] / (nNeutrons*nTheta);
			Qjet[jE] += ( fE * st.ntNeutron.injection.get(iRE) * vol );
			jE++;
		},{-1,iR.coord[DIM_R],0});
	},{0,-1,0});
	
	
	double pasoEn = pow(st.ntNeutron.emax()/st.ntNeutron.emin(),1.0/nE);
	double Qtot = 0.0;
	size_t jE = 0;
	st.ntNeutron.ps.iterate([&](const SpaceIterator& iE) {
		double E = iE.val(DIM_E);
		Qtot += Qjet[jE]*E*E*(pasoEn-1.0);
		jE++;
	},{-1,0,0});
	cout << "Total power injected in the jet by neutron decays = " << Qtot << endl;
	fileRandoms.close();
}