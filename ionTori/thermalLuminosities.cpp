// Standard C++ libraries
#include <iomanip>

// Boost libraries
#include <boost/property_tree/ptree.hpp>
#include <boost/math/tools/roots.hpp>


#include <nrMath/nrutil.h>


// Standard user-made libraries
#include <fparameters/SpaceIterator.h>
#include <fparameters/parameters.h>
#include <fluminosities/compton.h>
#include <fluminosities/thermalCompton.h>
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
			Matrix scatt, Vector& energies, const int flags[], Matrix& lumOut)
{
	double lumSync1,lumSync2;
	matrixInit4(lumOutSy,lumOutBr,lumOutpp,lumOut,nE,nR,0.0);
	
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
					//lumSy += jSy*emissToLum;
					lumSy += 2.0*min(jSy*emissToLum,bb_RJ(frecuency,temp)*fluxToLum);
					//lumRJ += fluxToLum*bb_RJ(frecuency,temp);
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
				else { lumSync2 = 2.0*min(lumSy,lumRJ);}
				//lumSy= max(lumSync2-lumSync1,0.0);
				//lumSy=lumSync2-lumSync1;
				if (lumSy > 0.0)
					lumOutSy[jE][jR]=lumSy;
				lumSync1=lumSync2;
            }
			if (flags[1]) {
				lumBr *= 2.0;
				if (flags[0] && frecuency < 1.0e8)
					lumBr=max(lumBr,lumSy);
				lumOutBr[jE][jR]=lumBr;
			}
			lumOut[jE][jR] = lumOutSy[jE][jR]+lumOutBr[jE][jR]+lumOutpp[jE][jR];
		jR++;
		},{itE.coord[DIM_E],-1,0});
	jE++;	
	},{-1,0,0});
}

double maxTemp(const State& st, const SpaceIterator& i)
{
	double temp1 = 0.0;
	double temp = 0.0;
	st.photon.ps.iterate([&](const SpaceIterator& itRT) {
		temp1 = st.tempElectrons.get(itRT);
		temp = (temp > temp1 ? temp : temp1);
	},{0,i.coord[DIM_R],-1});
	return temp;
}

void thermalCompton(const State& st, Matrix& lumOut, Matrix scatt,
				Matrix& lumOutIC, Vector energies)
{
	show_message(msgStart,Module_thermalCompton);
	
	Matrix lumOutLocal;
	Vector lumInIC(nE,0.0);
	matrixInit(lumOutIC,nE,nR,0.0);
	matrixInitCopy(lumOutLocal,nE,nR,lumOut);
			
	double res;
	size_t it=1;
	do {
		cout << "Iteration number = " << it << endl;
		res=0.0;
		
		// SETEAMOS A 0 EN CADA ITERACIÓN
		fill(lumInIC.begin(),lumInIC.end(),0.0);
				
		if(it > 1)
		{
			for (size_t jE=0;jE<nE;jE++) {
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
					lumInIC[jE] += scatt[jjR][jR]*lumOut[jE][jjR];
				}
			}
			
			double temp = maxTemp(st,itR);
			double normtemp = boltzmann*temp/(electronMass*cLight2);
			
			// CALCULAMOS LA COMPTONIZACION					 function compton()
			if (normtemp >= 0.001) 
			{
				for (size_t jjE=0;jjE<nE;jjE++) {   // para cada energía del espectro inicial
													// calculamos la contribución a la luminosidad
													// total a cada energía del espectro final
													// y sumamos.
					double frecuency=energies[jjE]/planck;
					if (lumInIC[jjE]*frecuency > 1.0e20 && energies[jjE] > 1.0e-18)
						comptonNew(lumOutIC,lumInIC[jjE],energies,temp,nE,jjE,jR);
				}
				for (size_t jE=0;jE<nE;jE++) {
					if (lumOutIC[jE][jR] > 0.0) {
						double lum = lumOutLocal[jE][jR] + lumOutIC[jE][jR];
				
						if (lumOut[jE][jR] > 0.0)
							res += abs(log10(lum)-log10(lumOut[jE][jR]));
						lumOut[jE][jR] = lum;
					}
				}
			}
			jR++;
		},{0,-1,0});
		
		res /= (nE*nR);
		cout << "Residuo = " << res << endl;
	++it;	
	} while (res > 1.0e-1);
	show_message(msgEnd,Module_thermalCompton);
}

void thermalCompton2(const State& st, Matrix& lumOut, Matrix scatt,
				Matrix& lumOutIC, Vector energies)
{
	show_message(msgStart,Module_thermalCompton);
	
	size_t nOm = 100;
	size_t nG = 10;
	Vector p(nR*nG*nE*nOm,0.0);
	
	size_t jR = 0;
	st.photon.ps.iterate([&](const SpaceIterator& itR) {
		double normtemp = boltzmann*maxTemp(st,itR)/(electronMass*cLight2);
		if (normtemp >= 0.01) {
			comptonRedistribution(p,nG,nE,nOm,jR,normtemp,energies,lumOut);
			
			double aux=0.;
		}
		jR++;
	},{0,-1,0});
	
	Matrix lumOutLocal;
	Vector lumInIC(nE,0.0);
	matrixInit(lumOutIC,nE,nR,0.0);
	matrixInitCopy(lumOutLocal,nE,nR,lumOut);
			
	double res;
	size_t it=1;
	do {
		cout << "Iteration number = " << it << endl;
		res=0.0;

		fill(lumInIC.begin(),lumInIC.end(),0.0);
				
		if(it > 1) {
			for (size_t jE=0;jE<nE;jE++) {
				for (size_t jR=0;jR<nR;jR++)
					lumOutIC[jE][jR]=0.0;
			}
		}
		
		// PARA CADA CELDA
		size_t jR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itR) {
			// CALCULAMOS LinC PARA CADA E        			 function LinC
			for (size_t jE=0;jE<nE;jE++) {
				for (size_t jjR=0;jjR<nR;jjR++)
					lumInIC[jE] += scatt[jjR][jR]*lumOut[jE][jjR];
			}
			double normtemp = boltzmann*maxTemp(st,itR)/(electronMass*cLight2);
			if (normtemp >= 0.01) 
			{
				for (size_t jjE=0;jjE<nE;jjE++) {
					double frecuency=energies[jjE]/planck;
					if (lumInIC[jjE]*frecuency > 1.0e20)
						comptonNew2(lumOutIC,lumInIC[jjE],energies,nG,nE,nOm,normtemp,p,jjE,jR);
				}
				for (size_t jE=0;jE<nE;jE++) {
					if (lumOutIC[jE][jR] > 0.0) {
						double lum = lumOutLocal[jE][jR] + lumOutIC[jE][jR];
						if (lumOut[jE][jR] > 0.0)
							res += abs(log10(lum/lumOut[jE][jR]));
						lumOut[jE][jR] = lum;
					}
				}
			}
			jR++;
		},{0,-1,0});
		
		res /= (nE*nR);
		cout << "Residuo = " << res << endl;
	++it;	
	} while (res > 1.0e-5);
	show_message(msgEnd,Module_thermalCompton);
}


void gravRedshift(const State& st, Vector energies, Matrix lum, Matrix& lumRed)
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
				lumRed[items[jE]][jR] += lum[jE][jR] * (gfactor*denergy_original/denergy_new);
			}
		}
		jR++;
	},{0,-1,0});
}

void binEmissivities(const State& st, Vector esc, Vector energies, Matrix lumOut)
{
	ofstream filePar, fileLum;
	FILE *fileEmissBin;
	
	fileLum.open("emissivities.txt", ios::out);
	filePar.open("parSpace.txt", ios::out);
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
			double rB1=r-dr/2.0;
			double rB2=r+dr/2.0;
			double vol=(rB2*rB2*rB2-rB1*rB1*rB1)*volFactor;
			double emissToLum=vol*4.0*pi;
			double sumFactor=0.0;
			st.photon.ps.iterate([&](const SpaceIterator& itTheta) {
				sumFactor += st.tempElectrons.get(itTheta);
			},{0,itR.coord[DIM_R],-1});
			st.photon.ps.iterate([&](const SpaceIterator& itTheta) {
				double lum = lumOut[jE][jR]*esc[jR]*st.tempElectrons.get(itTheta)/(sumFactor*2.0);
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

void writeLuminosities(State& st, Vector energies, Matrix lumOutSy, Matrix lumOutBr,
						Matrix lumOutpp, Matrix lumOutIC, Matrix lumOut, 
						Matrix lumOutRed, Vector esc, const string& filename)
{
	ofstream file1,file2;
	file1.open(filename.c_str(),ios::out);
	file2.open("lumRadius.txt",ios::out);
	
	double lumSy,lumBr,lumIC,lumpp,lum,lumRed;
	for (size_t jE=0;jE<nE;jE++) {
		double frecuency = energies[jE]/planck;
		double energyEV = energies[jE]/1.6e-12;
		lumSy = lumBr = lumIC = lumpp = lum = lumRed = 0.0;
		for (size_t jR=0;jR<nR;jR++) {
			lumSy += lumOutSy[jE][jR] * esc[jR];
			lumBr += lumOutBr[jE][jR] * esc[jR];
			lumIC += lumOutIC[jE][jR] * esc[jR];
			lumpp += lumOutpp[jE][jR] * esc[jR];
			lum += lumOut[jE][jR] * esc[jR];
			lumRed += lumOutRed[jE][jR] * esc[jR];
		}
        file1
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << energyEV
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumSy*frecuency
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumBr*frecuency
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumIC*frecuency
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumpp*frecuency
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lum*frecuency
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumRed*frecuency
			<< endl;
	};
	
	double eVar = pow(energies[nE-1]/energies[0],1.0/nE);
	size_t jR=0;
	st.photon.ps.iterate([&](const SpaceIterator& itR) {
		double r = itR.val(DIM_R);
		lumSy = lumBr = lumIC = lumpp = lum = 0.0;
		for (size_t jE=0;jE<nE;jE++)  {
			double dfrecuency = energies[jE]/planck * (eVar-1.0);
			lumSy += lumOutSy[jE][jR]*dfrecuency;
			lumBr += lumOutBr[jE][jR]*dfrecuency;
			lumIC += lumOutIC[jE][jR]*dfrecuency;
			lumpp += lumOutpp[jE][jR]*dfrecuency;
			lum += lumOut[jE][jR];
		}
		
		file2
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << r
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumSy
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumBr
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumIC
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumpp
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lum
			<< endl;
		jR++;
	},{0,-1,0});
	
	file1.close();
	file2.close();
}

void photonDensity(State& st, Vector energies, Matrix lumOut, Vector esc)
{
	double zGap = GlobalConfig.get<double>("zGap");
	
	ofstream file;
	file.open("photonDensity.txt",ios::out);
	file << "zGap = " << zGap << endl;
	file << "energy [eV] \t n_ph [cm^-3 erg^-1]" << endl << endl;

	size_t jE=0;
	st.photon.ps.iterate([&](const SpaceIterator& itE) {
		double nPh = 0.0;
		double energy = itE.val(DIM_E);
		size_t jR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itER) {
			double r = itER.val(DIM_R);
			double dist2 = (zGap*zGap + r*r) * gravRadius*gravRadius;
			nPh += lumOut[jE][jR]*esc[jR] / (pi*dist2*cLight*energy*planck);
			jR++;
		},{itE.coord[DIM_E],-1,0});
		file << energy/1.6e-12 << "\t" << nPh << endl;
		jE++;
	},{-1,0,0});
	file.close();
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
		localProcesses(st,lumOutSy,lumOutBr,lumOutpp,scatt,energies,processesFlags,lumOut);
		if (processesFlags[3])
			thermalCompton2(st,lumOut,scatt,lumOutIC,energies);
		gravRedshift(st,energies,lumOut,lumOutRed);
		binEmissivities(st,esc,energies,lumOut);
		photonDensity(st,energies,lumOut,esc);
		writeLuminosities(st,energies,lumOutSy,lumOutBr,lumOutpp,
							lumOutIC,lumOut,lumOutRed,esc,filename);
	}
	show_message(msgEnd,Module_thermalLuminosities);
}