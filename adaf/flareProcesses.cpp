#include "flareProcesses.h"
#include "globalVariables.h"
#include "adafFunctions.h"
#include "flareTimeDependent.h"
#include <flosses/lossesSyn.h>

#include "write.h"
#include <fparameters/SpaceIterator.h>
#include <inout/ioutil.h>
#include <boost/property_tree/ptree.hpp>
#include <iostream>
#include <fmath/RungeKutta.h>

using namespace std;

void flareProcesses(State& st) {

	maxRadius = GlobalConfig.get<double>("nonThermal.flare.maxRadius")*schwRadius;
	minRadius = GlobalConfig.get<double>("nonThermal.flare.minRadius")*schwRadius;
	etaInj = GlobalConfig.get<double>("nonThermal.flare.injection.energyFraction");
	pIndex = GlobalConfig.get<double>("nonThermal.flare.injection.primaryIndex");
	
	double t1 = 0.0;
	ofstream vel_file;
	vel_file.open("vel_r.dat",ios::out);
	double r1 = exp(logr.front())*schwRadius;
	double pasoR = pow(maxRadius/exp(logr.front())/schwRadius,1.0/1000);
	for (size_t jT=0;jT<1000;jT++) {
		r1 *= pasoR;
		t1 += r1*(pasoR-1.0)/(-radialVel(r1));
	}
	cout << "tiempo en llegar al horizonte = " << int(t1/3600) << "h" << int(t1/60)-60*int(t1/3600)
			<< "m" << endl;
	
	timeAfterFlare = 0.0;
	double tMax = 3.0*3600;
	size_t nTime = 50;
	double dt = tMax/nTime;
	ofstream file_cm,file_mm,fileNIR,fileX;
	multiZoneInjection(st.ntPair,st);
	file_cm.open("lightCurve_cm.txt");
	file_mm.open("lightCurve_mm.txt");
	fileNIR.open("lightCurveNIR.txt");
	fileX.open("lightCurveX.txt");
	
	
	size_t nEth = 100;
	size_t nRth = 40;
	Vector energies(nEth,0.0);
	Vector radii(nRth,0.0);
	Matrix num_dens;
	if (calculateThermal==0) {
		ifstream file_dens;
		file_dens.open("lumER.txt");
		file_dens >> nEth >> nRth;
		energies.resize(nEth);
		radii.resize(nRth);
		matrixInit(num_dens,nEth,nRth,0.0);
		for (size_t jE=0;jE<nEth;jE++) {
			for (size_t jR=0;jR<nRth;jR++) {
					file_dens >> energies[jE] >> radii[jR] >> num_dens[jE][jR];
			}
		}
		file_dens.close();
	}
	
	st.photon.ps.iterate([&](const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		size_t jR = 0;
		while (r > radii[jR])
			jR++;
		st.photon.ps.iterate([&](const SpaceIterator& iRE) {
			size_t jE = iRE.coord[DIM_E];
			double dist = 0.0;
			if (jR>0) {
				double n1 = num_dens[jE][jR-1];
				double n2 = num_dens[jE][jR];
				double r1 = radii[jR-1];
				double r2 = radii[jR];
				double m = safeLog10(n2/n1)/safeLog10(r2/r1);
				dist = pow(10.0,safeLog10(n1)+m*safeLog10(r/r1));
			} else
				dist = num_dens[jE][0];
			st.photon.distribution.set(iRE,dist);
		},{-1,iR.coord[DIM_R],0});
	},{0,-1,0});
	
	double Eradiated = 0.0;
	for (size_t i=0;i<nTime;i++) {
		timeAfterFlare += dt;
		//////////////////////////////////////////////////////
		
		multiZoneDist(st.ntPair,st,-timeAfterFlare);
		
		if (i % 10 == 0) {
			ofstream fileR;
			fileR.open("energyTime_"+to_string(i)+".txt");
			st.ntPair.ps.iterate([&](const SpaceIterator& iR) {
				double rB2 = iR.val(DIM_R)*sqrt(paso_r);
				double rB1 = rB2/paso_r;
				double vol = (4.0/3.0)*pi*costhetaH(iR.val(DIM_R))*(rB2*rB2*rB2-rB1*rB1*rB1);
				double energy = vol * RungeKuttaSimple(st.ntPair.emin(),st.ntPair.emax(),
					[&](double e) {return st.ntPair.distribution.interpolate({{0,e}},
					&iR.coord)*e;});
				fileR << safeLog10(iR.val(DIM_R)/schwRadius) << "\t" << safeLog10(energy) << endl;
			},{0,-1,0});
			fileR.close();
		}
		
		// TOTAL RADIATED ENERGY
		double dEdt = 0.0;
		st.ntPair.ps.iterate([&](const SpaceIterator& iR) {
			double rB2 = iR.val(DIM_R)*sqrt(paso_r);
			double rB1 = rB2/paso_r;
			double vol = (4.0/3.0)*pi*costhetaH(iR.val(DIM_R))*(rB2*rB2*rB2-rB1*rB1*rB1);
			double B = st.magf.get(iR);
			double dEdtLocal = vol * RungeKuttaSimple(st.ntPair.emin(),st.ntPair.emax(),
				[&](double e) {return st.ntPair.distribution.interpolate({{0,e}},
				&iR.coord)*lossesSyn(e,B,st.ntPair);});
			dEdt += dEdtLocal;
		},{0,-1,0});
		
		Eradiated += dEdt * dt;
		//////////////////////////////////////////////////////
		flareEmission2(st,st.ntPair,file_cm,file_mm,fileNIR,fileX);
			
		if (i==10) writeEandRParamSpace("electronDis",st.ntPair.distribution,0,0);
	}
	cout << "Total energy radiated = " << Eradiated << endl;

	file_cm.close();
	file_mm.close();
	fileNIR.close();
	fileX.close();
}