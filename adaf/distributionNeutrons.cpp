#include "distributionNeutrons.h"
#include <fparameters/parameters.h>
#include "globalVariables.h"
#include <fmath/RungeKutta.h>
#include <fparameters/Dimension.h>
#include <fparameters/parameters.h>
#include <fparameters/SpaceIterator.h>
#include <finjection/neutronDecay.h>
#include <finjection/pgammaPionInj.h>
#include <finjection/pairBH.h>
#include <flosses/lossesPhotoHadronic.h>
#include <flosses/lossesSyn.h>
#include "absorption.h"
#include "write.h"
#include "adafFunctions.h"

//extern "C" {
//	#include <nrMath/random.h>
//}

using namespace std;

void distributionNeutrons(State& st)
{
	ofstream fileNeutrons,fileProtons,fileElectrons;
	fileNeutrons.open("neutronProp.dat",ios::out);
	fileProtons.open("protonInjByNeutrons.dat",ios::out);
	fileElectrons.open("electronInjByNeutrons.dat",ios::out);
	
	fileNeutrons << "Rs = " << schwRadius << " cm" << endl;
	fileProtons << "Rs = " << schwRadius << " cm" << endl;
	fileElectrons << "Rs = " << schwRadius << " cm" << endl;
	
	Vector Q0n(nE,0.0);
	size_t jE = 0;
	st.ntNeutron.ps.iterate([&](const SpaceIterator &iE) {
		st.ntNeutron.ps.iterate([&](const SpaceIterator &iER) {
			double r = iER.val(DIM_R);
			Q0n[jE] += st.ntNeutron.injection.get(iER) * volume(r);
		},{iE.coord[DIM_E],-1,0});
		jE++;
	},{-1,0,0});

	size_t nRn = nE;
	double RnMin = cLight*neutronMeanLife*1.0e-4;
	double gammaMax = st.ntNeutron.emax() / (st.ntNeutron.mass * cLight2);
	double pasoRn = pow(gammaMax,1.0/nRn);
	
	double Qp_tot = 0.0;
	double Qe_tot = 0.0;
	double Rn = RnMin;
	for (size_t jR=0;jR<nRn;jR++) {
		Rn *= pasoRn;
		double dRn = Rn * (pasoRn - 1.0);
		double NnFactor = 1.0 / (4*pi*cLight) / (Rn*Rn);
		double vol = 4*pi*Rn*Rn*dRn;
		size_t jE = 0;
		st.ntNeutron.ps.iterate([&](const SpaceIterator &iE) {
			double En = iE.val(DIM_E);
			double g = En / (neutronMass*cLight2);
			double tDecay = neutronMeanLife*g;
			double rDec = cLight*tDecay;
			double Nn = NnFactor * Q0n[jE]*exp(-Rn/rDec);
			st.ntNeutron.distribution.set(iE,Nn);
			jE++;
			fileNeutrons << Rn / schwRadius << "\t" << En/EV_TO_ERG/1e9 << "\t"
						 << Nn << "\t"
						 << Nn*En*En * 4*pi*Rn*Rn << endl;
		},{-1,0,0});
		
		double Qp_local = 0.0;
		double pasoEp = pow(st.ntProton.emax()/st.ntProton.emin(),1.0/nE);
		st.ntProton.ps.iterate([&](const SpaceIterator &iE) {
			double Ep = iE.val(DIM_E);
			double En = Ep * 1.001;
			double tDecay = neutronMeanLife*(En/(neutronMass*cLight2));
			double Qp = P2(1.001) * st.ntNeutron.distribution.interpolate({{DIM_E,En}},&iE.coord) / tDecay;
			Qp_local += Qp*Ep*Ep*(pasoEp-1.0);
			fileProtons << safeLog10(Rn / schwRadius) << "\t" << safeLog10(Ep/1.6e-12) << "\t"
						<< safeLog10(Qp) << "\t"
						<< safeLog10(Qp*Ep*Ep*4.0*pi*Rn*Rn) << endl;
		},{-1,0,0});
		
		double Qe_local = 0.0;
		double pasoEe = pow(st.ntElectron.ps[DIM_E][nE-1]/st.ntElectron.ps[DIM_E][0],1.0/nE);
		st.ntElectron.ps.iterate([&](const SpaceIterator &iE) {
			double Ee = iE.val(DIM_E);
			double Qe = injElectronNeutronDecay(Ee,st.ntNeutron,iE);
			Qe_local += Qe*Ee*Ee*(pasoEe-1.0);
			fileElectrons << safeLog10(Rn) << "\t" << safeLog10(Ee/1.6e-12) << "\t"
						  << safeLog10(Qe) << "\t"
						  << safeLog10(Qe*Ee*Ee*4.0*pi*Rn*Rn) << endl;
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


void nonThermalPhotonDensity(State& st)
{
	ofstream file;
	file.open("nonThermalPhotonDensity_gap.dat",ios::out);
	
	size_t nZ = nR;
	double zMin = schwRadius;
	double pasoZ = pow(1e6,1.0/nZ);
	double z = zMin;
	st.photon.ps.iterate([&](const SpaceIterator& iZ) {
		z *= pasoZ;
		double Uph = 0.0;
		size_t jE=0;
		double pasoE = pow(st.ntPhoton.emax()/st.ntPhoton.emin(),1.0/(nE-1));
		st.ntPhoton.ps.iterate([&](const SpaceIterator& iZE) {
			double nPh = 0.0;
			double energy = iZE.val(DIM_E);
			double dE = energy * (pasoE - 1.0);
			st.ntPhoton.ps.iterate([&](const SpaceIterator& iZER) {
				double r = iZER.val(DIM_R);
				double height = height_fun(r);
				double dist2 = z*z+r*r;
				double tau_es = st.denf_e.get(iZER)*thomson*height;
				double tescape = height/cLight * (1.0 + tau_es);
				double lum = st.ntPhoton.distribution.get(iZER) * pow(tescape,-1) * volume(r);
				nPh += ( lum / (4*pi*dist2*cLight) );
			},{iZE.coord[DIM_E],-1,0});
			Uph += nPh * energy * dE;
			file << z/schwRadius << "\t" << energy << "\t" << nPh << endl;
			st.ntPhoton.injection.set(iZE,nPh);
			jE++;
		},{-1,iZ.coord[DIM_R],0});
		cout << "z = " << z/schwRadius << " Rs,\t Uph = " << Uph << " erg cm^-3" << endl;
	},{0,-1,0});
	file.close();
}

void distributionNeutronsAGN(State& st)
{
	ofstream fileNeutrons, fileProtons, fileProtonsPhotoMeson, fileElectrons, filePairsPhotoMeson,
				filePairsBH_dec, filePairsBH_photoMeson, filePairs_gg, filePairs_phph;
	fileNeutrons.open("neutronProp.dat",ios::out);
	fileProtons.open("protonInjByNeutrons.dat",ios::out);
	fileProtonsPhotoMeson.open("protonInjByNgamma.dat",ios::out);
	fileElectrons.open("electronInjByNeutrons.dat",ios::out);
	filePairsPhotoMeson.open("pairsInjByNgamma.dat", ios::out);
	filePairsBH_dec.open("pairsInjByBH_dec.dat", ios::out);
	filePairsBH_photoMeson.open("pairsInjByBH_photoMeson.dat", ios::out);
	filePairs_gg.open("pairsInj_gg.dat", ios::out);
	filePairs_phph.open("pairsInj_phph.dat", ios::out);
	
	fileNeutrons << "Rs = " << schwRadius << " cm" << endl;
	fileProtons << "Rs = " << schwRadius << " cm" << endl;
	fileProtonsPhotoMeson << "Rs = " << schwRadius << " cm" << endl;
	fileElectrons << "Rs = " << schwRadius << " cm" << endl;
	filePairsPhotoMeson << "Rs = " << schwRadius << " cm" << endl;
	filePairsBH_dec << "Rs = " << schwRadius << " cm" << endl;
	filePairsBH_photoMeson << "Rs = " << schwRadius << " cm" << endl;
	filePairs_gg << "Rs = " << schwRadius << " cm" << endl;
	filePairs_phph << "Rs = " << schwRadius << " cm" << endl;

	nonThermalPhotonDensity(st);
	
	size_t nZ = nR;
	double zMin = schwRadius;
	double gammaMax = st.ntNeutron.emax() / (st.ntNeutron.mass * cLight2);
	double pasoZ = pow(1.0e6,1.0/nZ);
	double z = zMin;
	for (size_t jZ=0;jZ<nZ;jZ++) {
		z *= pasoZ;
		st.ntNeutron.ps.iterate([&](const SpaceIterator &iE) {
			double En = iE.val(DIM_E);
			double g = En / (neutronMass*cLight2);
			double tDecay = neutronMeanLife*g;
			double rDec = cLight*tDecay;
			double Nn = 0.0;
			st.ntNeutron.ps.iterate([&](const SpaceIterator &iER) {
				double r = iER.val(DIM_R);
				double rEscape = sqrt(r*r+z*z);
				double NnFactor = 1.0 / (4*pi*cLight*rEscape*rEscape);
				Nn += NnFactor * st.ntNeutron.injection.get(iER) * volume(r) * exp(-rEscape/rDec);
			},{iE.coord[DIM_E],-1,0});
			st.ntNeutron.distribution.set(iE,Nn);
			fileNeutrons << z / schwRadius << "\t" << En/EV_TO_ERG/1e9 << "\t"
						 << Nn << "\t"
						 << Nn*En*En << endl;
		},{-1,0,0});
		
		st.ntNeutron.ps.iterate([&](const SpaceIterator& iR) {
			st.ntNeutron.ps.iterate([&](const SpaceIterator& iRE) {
				SpaceCoord i0 = {iRE.coord[DIM_E],0,0};
				double Nn = st.ntNeutron.distribution.get(i0);
				st.ntNeutron.distribution.set(iRE,Nn);
			},{-1,iR.coord[DIM_R],0});
		},{0,-1,0});
		
		double Nenergy = integSimpsonLog(st.ntNeutron.emin(), st.ntNeutron.emax(),
						[&](double En) {
							SpaceCoord distCoord = {0,0,0};
							return st.ntNeutron.distribution.interpolate({{DIM_E,En}},&distCoord)*En;
						},50);
		cout << "z = " << z/schwRadius << "\t Energy in neutrons = " << Nenergy << endl;
		
		st.ntProton.ps.iterate([&](const SpaceIterator &iE) {
			double Ep = iE.val(DIM_E);
			double En = Ep * 1.001;
			double tDecay = neutronMeanLife*(En/(neutronMass*cLight2));
			double Nn = (En > st.ntNeutron.emin() && En < st.ntNeutron.emax()) ?
						st.ntNeutron.distribution.interpolate({{DIM_E,En}},&iE.coord) : 0.0;
			double Qp = (En/Ep) * Nn / tDecay;
			double tEscape = z * 0.5 / cLight;
			st.ntProton.injection.set(iE,Qp*tEscape);
			fileProtons << z / schwRadius << "\t" << Ep / (EV_TO_ERG*1e9) << "\t"
						<< Qp << "\t"
						<< Qp*Ep*Ep << endl;
		},{-1,0,0});
		
		st.ntProton.ps.iterate([&](const SpaceIterator &iE) {
			SpaceCoord iEZ = {0,jZ,0};
			double Ep = iE.val(DIM_E);
			double En = Ep * 1.001;
			double tDecay = neutronMeanLife*(En/(neutronMass*cLight2));
			double Nn = (En > st.ntNeutron.emin() && En < st.ntNeutron.emax()) ?
						st.ntNeutron.distribution.interpolate({{DIM_E,En}},&iE.coord) : 0.0;
			double Qp = 0.5 * (En/Ep) * Nn * omegaPHsimple(En,st.ntNeutron,st.photon.injection,iEZ,
					st.photon.emin(),st.photon.emax());
			double tEscape = z * 0.5 / (0.9*cLight);
			st.ntProton.distribution.set(iE,Qp*tEscape);
			fileProtonsPhotoMeson << z / schwRadius << "\t" << Ep / (EV_TO_ERG*1e9) << "\t"
						<< Qp << "\t"
						<< Qp*Ep*Ep << endl;
		},{-1,0,0});
		
		st.ntElectron.ps.iterate([&](const SpaceIterator &iE) {
			double Ee = iE.val(DIM_E);
			double Qe = injElectronNeutronDecay(Ee,st.ntNeutron,iE);
			fileElectrons << z / schwRadius << "\t" << Ee / (EV_TO_ERG*1e6) << "\t"
						  << Qe << "\t"
						  << Qe*Ee*Ee <<  endl;
		},{-1,0,0});
		
		st.ntPair.ps.iterate([&](const SpaceIterator &iE) {
			double Ee = iE.val(DIM_E);
			SpaceCoord iEZ = {0,jZ,0};
			// VIA PHOTOMESON PRODUCTION
			/*double En = 20 * Ee;
			SpaceCoord iEZ = {0,jZ,0};
			double Nn = (En > st.ntNeutron.emin() && En < st.ntNeutron.emax()) ?
					st.ntNeutron.distribution.interpolate({{DIM_E,En}},&iE.coord) : 0.0;
			double rate_ng = lossesPhotoMeson(En,st.ntNeutron,st.photon.injection,iEZ,st.photon.emin(),
									st.photon.emax()) / En;
			double Qe_photoMeson = 1.0/8.0 * (En/Ee) * Nn * rate_ng;*/
			
			double magf = magneticField(1.5*schwRadius) * pow(z/schwRadius,-1.5);
			double Emu = 3.0 * Ee;
			double Epi = (4.0/3.0) * Emu;
			double En = 5.0 * Epi;
			double Nn = (En > st.ntNeutron.emin() && En < st.ntNeutron.emax()) ?
						st.ntNeutron.distribution.interpolate({{DIM_E,En}},&iE.coord) : 0.0;
			double Qpi = (En/Epi) * Nn * omegaPHsimple(En,st.ntNeutron,st.photon.injection,iEZ,
											st.photon.emin(), st.photon.emax()) * 0.875;
			double rateDecayPi = chargedPionMass*cLight2 / (Epi*chargedPionMeanLife);
			double rateLossesPi = lossesSyn(Epi,magf,st.ntChargedPion) / Epi;
			double Npi = (Epi > chargedPionMass*cLight2) ? Qpi * pow(rateDecayPi + rateLossesPi, -1) : 0.0;
			
			double Qmu = (Epi/Emu) * Npi * rateDecayPi;
			double rateDecayMu = muonMass*cLight2 / (Emu*muonMeanLife);
			double rateLossesMu = lossesSyn(Emu,magf,st.ntMuon) / Emu;
			double Nmu = (Emu > muonMass*cLight2) ? Qmu * pow(rateDecayMu + rateLossesPi, -1) : 0.0;
			
			double Qe_photoMeson = (Emu/Ee) * Nmu * rateDecayMu;
			
			filePairsPhotoMeson << z / schwRadius << "\t" << Ee / (EV_TO_ERG*1e6) << "\t"
						  << Qe_photoMeson << "\t"
						  << Qe_photoMeson*Ee*Ee << endl;
		},{-1,0,0});
		
		st.ntPair.ps.iterate([&](const SpaceIterator &iE) {
			SpaceCoord iEZ = {0,jZ,0};
			double Ee = iE.val(DIM_E);
			double ratio = protonMass/electronMass;
			double Ep = ratio * Ee;
			double Np = (Ep > st.ntProton.emin() && Ep < st.ntProton.emax()) ?
						st.ntProton.injection.interpolate({{DIM_E,Ep}},&iE.coord) : 0.0;
			double Qe_BH = 2.0 * ratio * Np * omegaBH2(Ep,st.ntProton,st.photon.injection,iEZ,
												st.photon.emin(), st.photon.emax());
			
			filePairsBH_dec << z / schwRadius << "\t" << Ee / (EV_TO_ERG*1e6) << "\t"
						  << Qe_BH << "\t"
						  << Qe_BH*Ee*Ee << endl;
		},{-1,0,0});
		
		st.ntPair.ps.iterate([&](const SpaceIterator &iE) {
			SpaceCoord iEZ = {0,jZ,0};
			double Ee = iE.val(DIM_E);
			double ratio = protonMass/electronMass;
			double Ep = ratio * Ee;
			double Np = (Ep > st.ntProton.emin() && Ep < st.ntProton.emax()) ?
						st.ntProton.distribution.interpolate({{DIM_E,Ep}},&iE.coord) : 0.0;
			double Qe_BH = 2.0 * ratio * Np * omegaBH2(Ep,st.ntProton,st.photon.injection,iEZ,
												st.photon.emin(), st.photon.emax());
			
			filePairsBH_photoMeson << z / schwRadius << "\t" << Ee / (EV_TO_ERG*1e6) << "\t"
						  << Qe_BH << "\t"
						  << Qe_BH*Ee*Ee << endl;
		},{-1,0,0});

		SpaceCoord iEZ = {0,jZ,0};
		double mu = 10.0 / sqrt( P2(10.0) + P2(z/schwRadius) );
		double cosAngle = 1.0-2.0*mu*mu;
		double nDot_gg = cLight * integSimpsonLog(st.ntPhoton.emin(),st.ntPhoton.emax(),
				[&st,cosAngle,&iEZ] (double Eg)
				{
					double Ng = (Eg > st.ntPhoton.emin() && Eg < st.ntPhoton.emax()) ?
									0.5 * st.ntPhoton.injection.interpolate({{DIM_E,Eg}},&iEZ) : 0.0;
					double aux = 2.0*P2(electronRestEnergy) / ( (1.0-cosAngle) * Eg );
					double EphMin = max(st.photon.emin(), aux);
					double result = (EphMin < st.photon.emax()) ?
							integSimpsonLog(st.photon.emin(), st.photon.emax(),
							[&st,cosAngle,Eg,&iEZ] (double Eph)
							{
								// El factor 0.5 para solo tener en cuenta la mitad que viene del otro lado.
								double Nph = (Eph > st.photon.emin() && Eph < st.photon.emax()) ?
								0.5*st.photon.injection.interpolate({{DIM_E,Eph}},&iEZ) : 0.0;
								double sigma = ggCrossSection(cosAngle,Eph,Eg);
								return Nph*sigma;
							},60) : 0.0;
					return Ng*result;
				},60);
		double nDot_phph = cLight * integSimpsonLog(st.photon.emin(),st.photon.emax(),
				[&st,cosAngle,&iEZ] (double Eg)
				{
					double Ng = (Eg > st.photon.emin() && Eg < st.photon.emax()) ?
									0.5 * st.photon.injection.interpolate({{DIM_E,Eg}},&iEZ) : 0.0;
					double aux = 2.0*P2(electronRestEnergy) / ( (1.0-cosAngle) * Eg );
					double EphMin = max(st.photon.emin(), aux);
					double result = (EphMin < st.photon.emax()) ?
							integSimpsonLog(st.photon.emin(), st.photon.emax(),
							[&st,cosAngle,Eg,&iEZ] (double Eph)
							{
								// El factor 0.5 para solo tener en cuenta la mitad que viene del otro lado.
								double Nph = (Eph > st.photon.emin() && Eph < st.photon.emax()) ?
								0.5*st.photon.injection.interpolate({{DIM_E,Eph}},&iEZ) : 0.0;
								double sigma = ggCrossSection(cosAngle,Eph,Eg);
								return Nph*sigma;
							},60) : 0.0;
					return 0.5*Ng*result;
				},60);
		filePairs_gg << z / schwRadius << "\t" << nDot_gg << endl;
		filePairs_phph << z / schwRadius << "\t" << nDot_phph << endl;
	}
	
	fileNeutrons.close();
	fileProtons.close();
	fileProtonsPhotoMeson.close();
	fileElectrons.close();
	filePairsPhotoMeson.close();
	filePairsBH_dec.close();
	filePairsBH_photoMeson.close();
	filePairs_gg.close();
	filePairs_phph.close();
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