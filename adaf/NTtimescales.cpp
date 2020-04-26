#include "NTtimescales.h"
#include "write.h"
#include "messages.h"
#include "globalVariables.h"
#include "adafFunctions.h"

#include <flosses/nonThermalLosses.h>
#include <flosses/lossesSyn.h>
#include <flosses/lossesIC.h>
#include <flosses/lossesBrem.h>
#include <flosses/lossesHadronics.h>
#include <flosses/lossesPhotoHadronic.h>
#include <fparameters/SpaceIterator.h>
#include <fparameters/Dimension.h>
#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>

void nonThermalTimescales(Particle& p, State& st, const std::string& filename)
{
	std::ofstream file,file2;
	file.open(filename.c_str(), std::ios::out);
	if (p.id == "ntProton") file2.open("diffusion_adv.dat", std::ios::out);

	file << "r [2M]" 
		<< "\t" << "Log(gamma)"
		<< "\t" << "Acc"
		<< "\t" << "Adv"
		<< "\t" << "Diff"
        << "\t" << "EmaxHillas"
		<< "\t" << "Sy"
		<< "\t" << "IC/pp"
		<< "\t" << "pg/Bremss"
		<< std::endl;
	
	if (p.id == "ntProton") {
		file2   << "r [2M]"
				<< "\t" << "Log(gamma)"
				<< "\t" << "Diff_length"
				<< "\t" << "dR"
				<< "\t" << "Dlenght/dR"
				<< std::endl;
	}
	
	int flag1,flag2,flag3,flag4,flag5;
	flag1 = flag2 = flag3 = flag4 = flag5 = 0;
	double logr1,logr2,logr3,logr4,logr5;
	logr1 = log10(1.5);
	double aux = log10(st.denf_e.ps[DIM_R].last()/schwRadius)/4.0;
	logr2 = logr1+aux;
	logr3 = logr2+aux;
	logr4 = logr3+aux;
	logr5 = logr4+aux;
	p.ps.iterate([&](const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		double logr = log10(r/schwRadius);
		if (logr > logr1) flag1++;
		if (logr > logr2) flag2++;
		if (logr > logr3) flag3++;
		if (logr > logr4) flag4++;
		if (logr > 0.9*logr5) flag5++;
		
		if (flag1 == 1 || flag2 == 1 || flag3 == 1 || flag4 == 1 || flag5 == 1) {
			double tAdv = accretionTime(r);
			double vR = radialVel(r/sqrt(paso_r));
			double dR = r * (sqrt(paso_r)-1.0/sqrt(paso_r));
			double B = st.magf.get(iR);
			double height = height_fun(r);
			double rho = massDensityADAF(r);
			double eMaxHillas = electronCharge*B*height;
			p.ps.iterate([&](const SpaceIterator& iRE) {
				double E = iRE.val(DIM_E);
				double tAcc = 1.0/accelerationRate(E,B);
				if (accMethod == 1) tAcc = accelerationTimeSDA(E,p,B,height,rho);
				double tDiff = diffusionTimeTurbulence(E,height,p,B);

				file << (int)(r/schwRadius)
					 << "\t" << safeLog10(E/(p.mass*cLight2))
					 << "\t" << safeLog10(tAcc)
					 << "\t" << safeLog10(tAdv)
					 << "\t" << safeLog10(tDiff)
					 << "\t" << safeLog10(eMaxHillas/1.602e-12/(p.mass*cLight2));
				if (p.id == "ntProton") {
					double diff_length = diffLength(E/(p.mass*cLight2),p,r,height,B,vR);
					file2   << (int)(r/schwRadius)
							<< "\t" << safeLog10(E/(p.mass*cLight2))
							<< "\t" << safeLog10(diff_length/schwRadius)
							<< "\t" << safeLog10(dR/schwRadius)
							<< "\t" << diff_length/dR << std::endl;
				}
				
				if (p.id == "ntElectron") {
					double tSyn = E/lossesSyn(E,B,p);
					double tIC = E/lossesIC(E,p,st.photon.distribution,iR.coord,st.photon.emin(),
											st.photon.emax());
					double tBrem = E/lossesBremss(E,st.denf_e.get(iR)+st.denf_i.get(iR),p);  
					file << "\t" << safeLog10(tSyn)
						 << "\t" <<safeLog10(tIC)
						 << "\t" << safeLog10(tBrem)
						 << std::endl;
				} else if(p.id == "ntProton") {
					double tSyn = E/lossesSyn(E,B,p);
					double tIC_Th = E/lossesIC_Th(E,p,st.photon.distribution,iR.coord,st.photon.emin(),
											st.photon.emax());
					double tPP = E/lossesHadronics(E,st.denf_i.get(iR),p);
					double tPG = E/lossesPhotoMeson(E,p,st.photon.distribution,iR,st.photon.emin(),
									st.photon.emax());
					double tBH = E/lossesPhotoPair(E,p,st.photon.distribution,iR,st.photon.emin(),
									st.photon.emax());
					file << "\t" << safeLog10(tSyn)
						 << "\t" << safeLog10(tIC_Th)
						 << "\t" << safeLog10(tPP)
						 << "\t" << safeLog10(tPG)
						 << "\t" << safeLog10(tBH)
						 << std::endl;
				}
			},{-1,iR.coord[DIM_R],0});
		}
	},{0,-1,0});
	file.close();
	if (p.id == "ntProton") file2.close();
}

using namespace std;

void secondariesTimescales(Particle& p, State& st, const std::string& filename)
{
	ofstream file;
	file.open(filename.c_str(), ios::out);

	file << "r [M]" 
		<< "\t" << "Log(gamma)"
		<< "\t" << "Adv"
		<< "\t" << "Diff"
		<< "\t" << "Decay"
		<< "\t" << "Sy"
		<< "\t" << "IC/pp"
		<< "\t" << "pg/Bremss"
		<< endl;

	int flag1,flag2,flag3,flag4,flag5;
	flag1 = flag2 = flag3 = flag4 = flag5 = 0;
	double logr1,logr2,logr3,logr4,logr5;
	logr1 = log10(1.5);
	double aux = log10(st.denf_e.ps[DIM_R].last()/schwRadius)/4.0;
	logr2 = logr1+aux;
	logr3 = logr2+aux;
	logr4 = logr3+aux;
	logr5 = logr4+aux;
	p.ps.iterate([&](const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		double logr = log10(r/schwRadius);
		if (logr > logr1) flag1++;
		if (logr > logr2) flag2++;
		if (logr > logr3) flag3++;
		if (logr > logr4) flag4++;
		if (logr > 0.9*logr5) flag5++;
		
		if (flag1 == 1 || flag2 == 1 || flag3 == 1 || flag4 == 1 || flag5 == 1) {
			double tAdv = accretionTime(r);
			double B = st.magf.get(iR);
			double height = height_fun(r);
			p.ps.iterate([&](const SpaceIterator& iRE) {
				double E = iRE.val(DIM_E);
				double tAcc = 1.0/accelerationRate(E,B);
				double tDiff = diffusionTimeTurbulence(E,height,p,B);
				double tSyn = E/lossesSyn(E,B,p);
				
				file << "\t" << (int)(r/schwRadius)
					 << "\t" << safeLog10(E/(p.mass*cLight2))
					 << "\t" << safeLog10(tAdv)
					 << "\t" << safeLog10(tDiff)
					 << "\t" << safeLog10(tSyn);
				
				double tDecay(0.0),tPP(1.0e30),tPG(0.0),tBH(0.0);
				if (p.id == "ntMuon") {
					tDecay = muonMeanLife*(E/(p.mass*cLight2));
					file << "\t" << safeLog10(tDecay) << endl;
				} else if (p.id == "ntChargedPion") {
					tDecay = chargedPionMeanLife*(E/(p.mass*cLight2));
					tPP = E/lossesHadronics(E,st.denf_i.get(iR),p);
					tPG = E/lossesPhotoMeson(E,p,st.photon.distribution,iR,st.photon.emin(),
									st.photon.emax());
					tBH = E/lossesPhotoPair(E,p,st.photon.distribution,iR,st.photon.emin(),
									st.photon.emax());
					file << "\t" << safeLog10(tDecay)
						 << "\t" << safeLog10(tPP)
						 << "\t" << safeLog10(tPG)
						 << "\t" << safeLog10(tBH) << endl;
				}
			},{-1,iR.coord[DIM_R],0});
		}
	},{0,-1,0});
	file.close();
}




void radiativeLossesNeutron(Particle& n, State& st, const std::string& filename)
{
	static const double accEfficiency = GlobalConfig.get<double>("nonThermal.injection.accEfficiency");
	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);

	file << "Log(E/GeV)" 
		<< "\t" << "r [M]"
		<< "\t" << "Escape"
		<< "\t" << "Decay"
		<< "\t" << "pp"
		<< "\t" << "pg"
		<< std::endl;
	
	double phEmin = st.photon.emin();
	double phEmax = st.photon.emax();
	n.ps.iterate([&](const SpaceIterator& i) {
		
		double E = i.val(DIM_E);
		double gamma = E/(neutronMass*cLight2);
		double r = i.val(DIM_R);
		double fmtE = log10(E/1.6e-3);
		
		double nEscape = cLight/r;
		double nDecay = 1.0/(gamma*neutronMeanLife);
		double nPP = lossesHadronics(E,st.denf_i.get(i),n)/E;
		double nPG = lossesPhotoHadronic(E,n,st.photon.distribution,i,phEmin,phEmax)/E;

		file << fmtE << "\t" << r/schwRadius
					 << "\t" << safeLog10(nEscape)
					 << "\t" << safeLog10(nDecay)
					 << "\t" << safeLog10(nPP)
				 	 << "\t" << safeLog10(nPG) << endl;
	},{-1,-1,0});
	file.close();
	
}

