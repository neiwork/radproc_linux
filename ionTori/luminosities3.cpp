#include "modelParameters.h"
#include "write.h"
#include "luminosities3.h"
#include "messages.h"
#include <fparameters/SpaceIterator.h>
#include <fparameters/parameters.h>
#include <boost/property_tree/ptree.hpp>
#include <fluminosities/compton.h>
#include <fluminosities/thermalSync.h>
#include <fluminosities/thermalBremss.h>
#include <fluminosities/thermalIC.h>
#include <fluminosities/blackBody.h>
#include "metric.h"
#include <fmath/mathFunctions.h>
#include <fmath/physics.h>
#include <boost/math/tools/roots.hpp>

extern "C" {
	#include <nrMath/nrutil.h>
}

#include <iomanip>
using namespace std;

void luminosities3(State& st, const string& filename, Matrix &prob) {
   
	show_message(msgStart, Module_luminosities);
    ofstream file;
	file.open(filename.c_str(), ios::out);
	/*file  << setw(10) << "frecuency [Hz]"
          << setw(10) << "Flux [erg s^-1 cm^-2]"
          << std::endl;
	*/	  
	double rg=GlobalConfig.get<double>("rg");
	int nE=GlobalConfig.get<int>("model.particle.default.dim.energy.samples");
    int nR=GlobalConfig.get<int>("model.particle.default.dim.radius.samples");
	int nTheta = GlobalConfig.get<int>("model.particle.default.dim.theta.samples");
    double rmax=GlobalConfig.get<double>("rEdge");
    double rmin=GlobalConfig.get<double>("rCusp");
    double thetamin = GlobalConfig.get<double>("thetamin");
    double thetamax = GlobalConfig.get<double>("thetamax");
	double dr=(rmax-rmin)/nR;
	double dtheta=(thetamax-thetamin)/nTheta;
    
    double zgap=GlobalConfig.get<double>("zgap");

	double lumBr;
	double lumSync1=0.0;
	double lumSync2=0.0;
	double lumSy;
	Vector lumSyTot(nE,0.0);
	Vector lumBrTot(nE,0.0);
	Vector energies(nE,0.0);
				
	Matrix lumOut;
	Matrix lumOutIC;
	Matrix lumOutSy;
    Matrix lumOutBr;
	matrixInit(lumOutSy,nE,nR,0.0);
    matrixInit(lumOutBr,nE,nR,0.0);
	matrixInit(lumOutIC,nE,nR,0.0);

	int jE=0;
	double areaFactor=2.0*pi*rg*rg;
	double volFactor=areaFactor*rg/3.0;
	
	int jR=0;
	st.photon.ps.iterate([&](const SpaceIterator& itR) { // LOOP EN RADIOS
		double r=itER.val(DIM_R);
		double rB1=r-dr/2.0;
		double rB2=r+dr/2.0;
		double area=rB2*rB2*areaFactor*dtheta;
		double fluxToLum=area;
		double vol=(rB2*rB2*rB2-rB1*rB1*rB1)*volFactor*dtheta;
		double emissToLum=vol*4.0*pi;
		double lumBrCell=0.0;
		double lumSyCell=0.0;
		double lumRJCell=0.0;
		
		st.photon.ps.iterate([&](const SpaceIterator& itRT) {
			double temp=st.tempElectrons.get(itERT);
			double magf=st.magf.get(itERT);
			double dens_i=st.denf_i.get(itERT);
			double dens_e=st.denf_e.get(itERT);
			
			st.photon.ps.iterate([&](const SpaceIterator& itRTE) {
			lumBrCell += jBremss(energies[jE],temp,dens_i,dens_e)*emissToLum;
			double jSy=0.0;
			if (temp >= 1.0e8) {
				jSy = jSync(energies[jE],temp,magf,dens_e);
			}
			lumSyCell += jSy*emissToLum;
			lumRJCell += bb(frecuency,temp)*fluxToLum;
		},{itE.coord[DIM_E],itER.coord[DIM_R],-1});
	
	st.photon.ps.iterate([&](const SpaceIterator& itE) { // LOOP EN ENERGÍAS
		energies[jE]=itE.val(DIM_E);
		double frecuency=energies[jE]/planck;
		int jR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itER) { // LOOP EN RADIOS
			double r=itER.val(DIM_R);
			double rB1=r-dr/2.0;
			double rB2=r+dr/2.0;
			double area=rB2*rB2*areaFactor*dtheta;
			double fluxToLum=area;
			double vol=(rB2*rB2*rB2-rB1*rB1*rB1)*volFactor*dtheta;
			double emissToLum=vol*4.0*pi;
			double lumBrCell=0.0;
			double lumSyCell=0.0;
			double lumRJCell=0.0;
			st.photon.ps.iterate([&](const SpaceIterator& itERT) {
				double temp=st.tempElectrons.get(itERT);
				double magf=st.magf.get(itERT);
				double dens_i=st.denf_i.get(itERT);
				double dens_e=st.denf_e.get(itERT);
				lumBrCell += jBremss(energies[jE],temp,dens_i,dens_e)*emissToLum;
				double jSy=0.0;
				if (temp >= 1.0e8) {
					jSy = jSync(energies[jE],temp,magf,dens_e);
				}
				lumSyCell += jSy*emissToLum;
				lumRJCell += bb(frecuency,temp)*fluxToLum;
			},{itE.coord[DIM_E],itER.coord[DIM_R],-1});
            if (jR>0) {
				double a=1.0-prob[jR-1][jR];
				double b=a*lumSync1;
				double c=b+2.0*lumSyCell;
                lumSync2 = min(c,2.0*lumRJCell);
            } else {
                lumSync2 = 2.0*min(lumSyCell,lumRJCell);
            }
            lumOutBr[jE][jR] *= 2.0;
			lumOutSy[jE][jR] = max(lumSync2-lumSync1,0.0);
			lumSync1=lumSync2;
            if (lumOutBr[jE][jR] > lumOutSy[jE][jR] && frecuency < 1.0e8) lumOutBr[jE][jR] = lumOutSy[jE][jR];
			lumOut[jE][jR]=lumOutSy[jE][jR]+lumOutBr[jE][jR];
			lumSyTot[jE] += lumOutSy[jE][jR];
			lumBrTot[jE] += lumOutBr[jE][jR];
		jR++;
		},{itE.coord[DIM_E],-1,0});
	jE++;	
	},{-1,0,0});
	
	Vector lumInIC(nE,0.0);
	Vector lumOutICvec(nE,0.0);
	
	for (int it=1;it<=3;it++) {      // ITERACIONES PARA CONVERGENCIA DEL COMPTON
	
	if (it > 1) {
		for (int jE=0;jE<nE;jE++){
			lumInIC[jE]=0.0;
			lumOutICvec[jE]=0.0;
		}
	}
	
	int jR=0;
	st.photon.ps.iterate([&](const SpaceIterator& itR) {
		for (int jjE=0;jjE<nE;jjE++) {
			double sumSy=0.0;
			double sumBr=0.0;
			double sumTot=0.0;
			for (int jjR=0;jjR<nR;jjR++) {
				sumTot += prob[jjR][jR]*lumOut[jjE][jjR];
			}
			lumInIC[jjE] = sumTot;      // La suma de ambos es L_C,in,j de Narayan
		}
		double temp=0.0;
		double temp1=0.0;
		st.photon.ps.iterate([&](const SpaceIterator& itRT) {
			temp1 = st.tempElectrons.get(itRT);
			temp = (temp > temp1 ? temp : temp1);
		},{0,itR.coord[DIM_R],-1});
		for (int jE=0;jE<nE;jE++) {
			lumOutIC[jE][jR] = compton(lumInIC,energies,temp,nE,jE);
			lumOutICvec[jE] += lumOutIC[jE][jR];
			lumOut[jE][jR] = lumOutIC[jE][jR]+lumOutBr[jE][jR]+lumOutSy[jE][jR];
		}
		
		jR++;
	},{0,-1,0});

	}
	
	double distance = 100.0*rg;
	Vector flux(nE,0.0);
    Vector nph(nE,0.0);
    double lum=0.0;
	for (int jjE=0;jjE<nE;jjE++) {
		double sum=0.0;
        double denergy=energies[jjE+1]-energies[jjE];
        int jjR=0.0;
		st.photon.ps.iterate([&](const SpaceIterator& itR) {
			sum += lumOut[jjE][jjR];
            lum += lumOut[jjE][jjR]*(denergy/planck);
            double r=itR.val(DIM_R);
            double d = sqrt(r*r+zgap*zgap)*rg;
            double redshifted = energies[jjE] / sqrt(-g_tt(zgap,pi/2.01));
            nph[jjE] += lumOut[jjE][jjR]/planck/(cLight*pi*d*d*redshifted);
            jjR++;
		},{0,-1,0});
		flux[jjE] = sum; // * 0.25 / pi / (distance*distance);
	}
	printf("lumtot = %5.5e\n",lum);

	for(jE=0;jE<nE;jE++) {
		double frecuency = energies[jE]/planck;
        file << "\t" << frecuency //energies[jE] / 1.6021e-12 / sqrt(-g_tt(zgap,pi/2.01))
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumSyTot[jE]*frecuency
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumBrTot[jE]*frecuency
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumInIC[jE]*frecuency
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumOutICvec[jE]*frecuency
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << flux[jE]*frecuency
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << nph[jE]
			 << std::endl;
	};
		
	file.close();
	show_message(msgEnd, Module_luminosities);

}