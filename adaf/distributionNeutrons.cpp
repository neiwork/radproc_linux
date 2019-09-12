#include "distributionNeutrons.h"
#include <fparameters/parameters.h>
#include "globalVariables.h"
#include <fparameters/Dimension.h>
#include <fparameters/parameters.h>
#include <fparameters/SpaceIterator.h>
#include <finjection/neutronDecay.h>
#include "write.h"
#include "adafFunctions.h"

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
	double r = 100.0*schwRadius;
	double pasoRn = pow(1.0e14,1.0/nRn);
	
	double Qp_tot = 0.0;
	double Qe_tot = 0.0;
	for (size_t jR=0;jR<nRn;jR++) {
		r *= pasoRn;
		double NnFactor = 1.0 / (4.0 * pi * cLight) / (r*r);
		double vol = 4.0*r*r*r*(pasoRn-1.0);
		size_t jE = 0;
		st.ntNeutron.ps.iterate([&](const SpaceIterator &iE) {
			double E = iE.val(DIM_E);
			double gamma = E / (neutronMass*cLight2);
			double tDecay = neutronMeanLife*gamma;
			double rDec = cLight*tDecay;
			double Nn = NnFactor * Q0n[jE]*exp(-r/rDec);
			st.ntNeutron.distribution.set(iE,Nn);
			jE++;
			fileNeutrons << safeLog10(r) << "\t" << safeLog10(E/1.6e-12) << "\t" << safeLog10(Nn) << endl;
		},{-1,0,0});
		
		double Qp_local = 0.0;
		double pasoEp = pow(st.ntProton.ps[DIM_E][nE-1]/st.ntProton.ps[DIM_E][0],1.0/nE);
		st.ntProton.ps.iterate([&](const SpaceIterator &iE) {
			double E = iE.val(DIM_E);
			double Qp = injNeutronDecay(E,st.ntProton,st.ntNeutron,iE);
			double tDecay = neutronMeanLife*(E/(neutronMass*cLight2));
			Qp = st.ntNeutron.distribution.get(iE)/tDecay;
			Qp_local += Qp*E*E*(pasoEp-1.0);
			fileProtons << safeLog10(r) << "\t" << safeLog10(E/1.6e-12) << "\t" << safeLog10(Qp) << endl;
		},{-1,0,0});
		
		double Qe_local = 0.0;
		double pasoEe = pow(st.ntElectron.ps[DIM_E][nE-1]/st.ntElectron.ps[DIM_E][0],1.0/nE);
		st.ntElectron.ps.iterate([&](const SpaceIterator &iE) {
			double E = iE.val(DIM_E);
			double Qe = injNeutronDecay(E,st.ntElectron,st.ntNeutron,iE);
			Qe_local += Qe*E*E*(pasoEe-1.0);
			fileElectrons << safeLog10(r) << "\t" << safeLog10(E/1.6e-12) << "\t" << safeLog10(Qe) << endl;
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

void jetNeutronDecay(State& st)
{
	double openAngle = GlobalConfig.get<double>("jetOpenAngle");
}