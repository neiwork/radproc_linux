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
	std::ofstream file;
	file.open(filename.c_str(), std::ios::out);

	file << "r [Rs]" 
		<< "\t" << "E [GeV]"
		<< "\t" << "Escape"
		<< "\t" << "Decay"
		<< "\t" << "pp"
		<< "\t" << "pg"
		<< std::endl;
	
	double phEmin = st.photon.emin();
	double phEmax = st.photon.emax();
	
	int flag1,flag2,flag3,flag4,flag5;
	flag1 = flag2 = flag3 = flag4 = flag5 = 0;
	double logr1,logr2,logr3,logr4,logr5;
	logr1 = log10(1.5);
	double aux = log10(st.denf_e.ps[DIM_R].last()/schwRadius)/4.0;
	logr2 = logr1+aux;
	logr3 = logr2+aux;
	logr4 = logr3+aux;
	logr5 = logr4+aux;
	n.ps.iterate([&](const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		double logr = log10(r/schwRadius);
		if (logr > logr1) flag1++;
		if (logr > logr2) flag2++;
		if (logr > logr3) flag3++;
		if (logr > logr4) flag4++;
		if (logr > 0.9*logr5) flag5++;
		
		if (flag1 == 1 || flag2 == 1 || flag3 == 1 || flag4 == 1 || flag5 == 1) {
			n.ps.iterate([&](const SpaceIterator& iRE) {
				double En = iRE.val(DIM_E);
				double gamma_n = En / (n.mass*cLight2);
				double tEscape = height_fun(r) / cLight;
				double tDecay = gamma_n*neutronMeanLife;
				double loss_np = lossesHadronics(En,st.denf_i.get(iR),n);
				double loss_ng = lossesPhotoMeson(En,n,st.photon.distribution,iR,phEmin,phEmax);
				double tNP = (loss_np > 0.0) ? En / loss_np : 1e30;
				double tNG = (loss_ng > 0.0) ? En / loss_ng : 1e30;

				file << r/schwRadius << "\t" << En / (EV_TO_ERG*1e9)
						 << "\t" << tEscape
						 << "\t" << tDecay
						 << "\t" << tNP
						 << "\t" << tNG << endl;
			},{-1,iR.coord[DIM_R],0});
		}
	},{0,-1,0});
	file.close();
}

