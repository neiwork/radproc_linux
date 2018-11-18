// Standard C++ libraries
#include <iomanip>

// Boost libraries
#include <boost/property_tree/ptree.hpp>
#include <boost/math/tools/roots.hpp>

// Numerical Recipes
extern "C" {
	#include <nrMath/nrutil.h>
}

// Standard user-made libraries
#include <fparameters/SpaceIterator.h>
#include <fparameters/parameters.h>
#include <fluminosities/compton.h>
#include <fluminosities/thermalSync.h>
#include <fluminosities/thermalBremss.h>
#include <fluminosities/luminosityHadronic.h>
#include <fluminosities/thermalIC.h>
#include <fluminosities/blackBody.h>
#include <fmath/mathFunctions.h>
#include <fmath/physics.h>
#include <fmath/fbisection.h>

// Project headers
#include "globalVariables.h"
#include "messages.h"
#include "metric.h"
#include "modelParameters.h"
#include "read.h"
#include "thermalLuminosities.h"
#include "torusFunctions.h"
#include "write.h"

double volFactor,areaFactor,dr,dtheta;

// Namespaces
using namespace std;

void localProcesses(const State& st, Matrix& lumOutSy, Matrix& lumOutBr, Matrix& lumOutpp, 
			Matrix scatt, Vector& energies, const int flags[])
{
	double lumSync1,lumSync2;
	if (flags[0])
		matrixInit(lumOutSy,nE,nR,0.0);
		lumSync1=lumSync2=0.0;
	if (flags[1])
		matrixInit(lumOutBr,nE,nR,0.0);
	if (flags[2])
		matrixInit(lumOutpp,nE,nR,0.0);
	
	size_t jE=0;
	st.photon.ps.iterate([&](const SpaceIterator& itE) {
		energies[jE] = itE.val(DIM_E);
		double frecuency=energies[jE]/planck;
		size_t jR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itER) {
			double r=itER.val(DIM_R);
			double rB1=r-dr/2.0;
			double rB2=r+dr/2.0;
			double area=rB2*rB2*areaFactor*dtheta;
			double fluxToLum=area;
			double vol=(rB2*rB2*rB2-rB1*rB1*rB1)*volFactor*dtheta;
			double emissToLum=vol*4.0*pi;
			double lumBr=0.0;
			double lumSy=0.0;
			double lumRJ=0.0;
			double lumpp=0.0;
			st.photon.ps.iterate([&](const SpaceIterator& itERT) {
				double temp=st.tempElectrons.get(itERT);
				double temp_i=st.tempIons.get(itERT);
				double magf=st.magf.get(itERT);
				double dens_i=st.denf_i.get(itERT);
				double dens_e=st.denf_e.get(itERT);
				if (flags[0]) {
					double jSy = jSync(energies[jE],temp,magf,dens_e);
					lumSy += jSy*emissToLum;
					lumRJ += fluxToLum*bb_RJ(frecuency,temp);
				}
				if (flags[1])
					lumBr += jBremss(energies[jE],temp,dens_i,dens_e)*emissToLum;
				if (flags[2] && temp_i > 1.0e11 && frecuency > 1.0e20 && frecuency < 1.0e26) {
					double jpp = luminosityHadronic(energies[jE],dens_i,temp_i);
					lumOutpp[jE][jR] += jpp*emissToLum;
				}
			},{itE.coord[DIM_E],itER.coord[DIM_R],-1});
            if (flags[0]) {
				if (jR>0)
					lumSync2 = min((1.0-scatt[jR-1][jR])*lumSync1+2.0*lumSy,2.0*lumRJ);
				else
					lumSync2 = 2.0*min(lumSy,lumRJ);
				lumSy= max(lumSync2-lumSync1,0.0);
				lumOutSy[jE][jR]=lumSy;
				//lumSy=lumSync2-lumSync1;
				lumSync1=lumSync2;
            }
			if (flags[1]) {
				lumBr *= 2.0;
				if (flags[0] && frecuency < 1.0e8)
					lumBr=max(lumBr,lumSy);
				lumOutBr[jE][jR]=lumBr;
			}
		jR++;
		},{itE.coord[DIM_E],-1,0});
	jE++;	
	},{-1,0,0});
}

double maxTempNorm(const State& st, const SpaceIterator& i)
{
	double temp1 = 0.0;
	double temp = 0.0;
	st.photon.ps.iterate([&](const SpaceIterator& itRT) {
		temp1 = st.tempElectrons.get(itRT);
		temp = (temp > temp1 ? temp : temp1);
	},{0,i.coord[DIM_R],-1});
	return temp;
}

void thermalCompton(const State& st, Matrix lumOutSy, Matrix lumOutBr, Matrix lumOutpp, 
			Matrix scatt, Matrix& lumOutIC, Vector energies)
{
	show_message(msgStart,Module_thermalCompton);
	
	matrixInit(lumOutIC,nE,nR,0.0);
	Matrix lumOutLocalConst,lumOutLocal;
	matrixInit(lumOutLocalConst,nE,nR,lumOutSy,lumOutBr,lumOutpp);
	matrixInit(lumOutLocal,nE,nR,lumOutLocalConst);
			
	double res;
	size_t it=1;
	Vector lumInIC(nE,0.0);
	do {
		cout << "Iteration number = " << it << endl;
		res=0.0;
		
		// SETEAMOS A 0 EN CADA ITERACIÓN
		fill(lumInIC.begin(),lumInIC.end(),0.0);
				
		if(it > 1)
		{
			for (size_t jE=0;jE<nE;jE++) {
				lumInIC[jE]=0.0;
				for (size_t jR=0;jR<nR;jR++) {
					lumOutIC[jE][jR]=0.0;
				}
			}
		}
		
		// PARA CADA CELDA
		size_t jR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itR) {
			// CALCULAMOS LinC PARA CADA E        			 function LinC
			for (size_t jE=0;jE<nE;jE++) {
				for (size_t jjR=0;jjR<nR;jjR++) {
					lumInIC[jE] += scatt[jjR][jR]*lumOutLocal[jE][jjR];
				}
			}
			
			double temp = maxTempNorm(st,itR);
			double normtemp = boltzmann*temp/(electronMass*cLight2);
			
			// CALCULAMOS LA COMPTONIZACION					 function compton()
			if (normtemp >= 0.1) 
			{
				for (size_t jjE=0;jjE<nE;jjE++) {   // para cada energía del espectro inicial
													// calculamos la contribución a la luminosidad
													// total a cada energía del espectro final
													// y sumamos.
					double frecuency=energies[jjE]/planck;
					if (lumInIC[jjE]*frecuency > 1.0e20) 
						compton2(lumOutIC,lumInIC[jjE],energies,temp,nE,jjE,jR);
				}
			}
			
			for (size_t jE=0;jE<nE;jE++) {
				if (lumOutLocal[jE][jR] > 0.0)
					res += log10(1+lumOutIC[jE][jR]/lumOutLocal[jE][jR]);
				lumOutLocal[jE][jR] = lumOutLocalConst[jE][jR] + lumOutIC[jE][jR];
			}
			jR++;
		},{0,-1,0});
		
		res /= (nE*nR);
		cout << "Residuo = " << res << endl;
	++it;	
	} while (res > 1.0e-2);
	show_message(msgEnd,Module_thermalCompton);
}

void gravRedshift(const State& st, Vector energies, Matrix lumSy, Matrix lumBr, Matrix lumpp,
			Matrix lumIC, Matrix& lumRed)
{
	matrixInit(lumRed,nE,nR,0.0);

	//DEFINO FRONTERAS PARA LAS CELDAS DE ENERGIA
	double eVar = pow(energies[nE-1]/energies[0],1.0/(nE-1));
	Vector benergies(nE+1,0.0);
	benergies[0] = energies[0]/sqrt(eVar);
	benergies[nE] = energies[nE-1]*sqrt(eVar);
	for (size_t jE=1;jE<nE;jE++){
		benergies[jE] = sqrt(energies[jE-1]*energies[jE]);
	}

	/////////////////////////////////////////////////////////////////
	// CALCULO LA REDISTRIBUCION EN ENERGIAS POR EFECTOS 
	// GRAVITACIONALES
	
	size_t jR=0;
	st.photon.ps.iterate([&](const SpaceIterator& itR) {
		double r = itR.val(DIM_R);
		double gfactor = redshiftFactor(r,0.0);
		size_t items[nE];
		for (size_t jE=0;jE<nE;jE++) {
			double renergy = energies[jE] * gfactor;
			size_t count=0;
			size_t jjE=0;
			while(jjE < nE && count == 0) {
				if (renergy > benergies[jjE] && renergy < benergies[jjE+1]) {
					items[jE] = jjE;
					count++;
				}
				jjE++;
			}
			if (count == 0) items[jE]=nE;
		}
		for (size_t jE=0;jE<nE;jE++) {
			if (items[jE] < nE) {
				double denergy_original = benergies[jE+1]-benergies[jE];
				double denergy_new = benergies[items[jE]+1]-benergies[items[jE]];
				double lumOut = lumSy[jE][jR]+lumBr[jE][jR]+lumpp[jE][jR]+lumIC[jE][jR];
				lumRed[items[jE]][jR] += lumOut * (gfactor*denergy_original/denergy_new);
			}
		}
		jR++;
	},{0,-1,0});
}

void binEmissivities(const State& st, Vector energies, Matrix lumSy, Matrix lumBr, Matrix lumpp,
				Matrix lumIC)
{
	ofstream filePar, fileLum;
	fileLum.open("emissivities.txt", ios::out);
	filePar.open("parSpace.txt", ios::out);
	FILE *fileEmissBin;
	fileEmissBin=fopen("emissivities.bin","wb");

	float *emissivities;
	emissivities=(float*)calloc(nE*nR*nTheta,sizeof(float));
	
	size_t jj=0;
	double eVar = pow(energies[nE-1]/energies[0],1.0/nE);
	for (size_t jE=0;jE<nE;jE++) {
		double sum=0.0;
        double denergy=energies[jE]*(eVar-1.0);
        size_t jR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itR) {
			double r = itR.val(DIM_R);
			double thetamax = fbisection([&](double x){return normalizedPotential(r,x);},0.1,pi/2.1,1.0e-3);
			double rB1=r-dr/2.0;
			double rB2=r+dr/2.0;
			double vol=(rB2*rB2*rB2-rB1*rB1*rB1)*volFactor;
			double emissToLum=vol*4.0*pi;
			double sumFactor=0.0;
			double lumOut = lumSy[jE][jR]+lumBr[jE][jR]+lumpp[jE][jR]+lumIC[jE][jR];
			st.photon.ps.iterate([&](const SpaceIterator& itTheta) {
				sumFactor += st.tempElectrons.get(itTheta);
			},{0,itR.coord[DIM_R],-1});
			st.photon.ps.iterate([&](const SpaceIterator& itTheta) {
				double lum = lumOut*st.tempElectrons.get(itTheta) / sumFactor / 2.0;
				fileLum << jE << "\t" << jR << "\t" 
						  << lum << endl;
				emissivities[jj++] = (float)(lum/emissToLum);
			},{0,itR.coord[DIM_R],-1});
            jR++;
		},{0,-1,0});
	}
	fwrite(emissivities,sizeof(float),jj,fileEmissBin);
	
	float deltaLogE = log10(energies[1]/energies[0]);
	filePar << nE << "\t" << nR << "\t" << nTheta << endl;
	filePar << cuspRadius << "\t" << dr << "\t" << minPolarAngle << "\t"
			<< dtheta << "\t" << log10(energies[0]) << "\t" << deltaLogE << endl;
	filePar << gravRadius << "\t" << blackHoleSpinPar << "\t" 
			<< specificAngMom << "\t" << torusCenterRadius << endl;
	
	fclose(fileEmissBin);
	fileLum.close();
	filePar.close();
}

void thermalLuminosities(State& st, const string& filename, Matrix &scatt, Vector& esc)
{
	show_message(msgStart,Module_thermalLuminosities);

	Matrix lumOutSy,lumOutBr,lumOutpp,lumOutIC,lumOut,lumOutRed;
	Vector energies(nE,0.0);
	
	dr = (edgeRadius-cuspRadius)/nR;
	dtheta = (maxPolarAngle-minPolarAngle)/nTheta;
	areaFactor=2.0*pi*gravRadius*gravRadius;
	volFactor=areaFactor*gravRadius/3.0;
	int processesFlags[numProcesses];
	
	readThermalProcesses(processesFlags);
	if (processesFlags[0] || processesFlags[1] || processesFlags[2]) {
		localProcesses(st,lumOutSy,lumOutBr,lumOutpp,scatt,energies,processesFlags);
		if (processesFlags[3])
			thermalCompton(st,lumOutSy,lumOutBr,lumOutpp,scatt,lumOutIC,energies);
		gravRedshift(st,energies,lumOutSy,lumOutBr,lumOutpp,lumOutIC,lumOutRed);
		binEmissivities(st,energies,lumOutSy,lumOutBr,lumOutpp,lumOutIC);
	}
	
/*	
	ofstream file;
	file.open(filename.c_str(), ios::out);
	for(size_t jE=0;jE<nE;jE++) {
		double frecuency = energies[jE]/planck;
		double eVenergy = energies[jE]/1.6e-12;
        file << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << eVenergy
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumSyTot[jE]*frecuency
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumBrTot[jE]*frecuency
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumICinTot[jE]*frecuency
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumICout_vec[jE]*frecuency
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumppTot[jE]*frecuency
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumOutTot[jE]*frecuency
			 << std::endl;
	};
		
	file.close(); */
	
	show_message(msgEnd,Module_thermalLuminosities);
}