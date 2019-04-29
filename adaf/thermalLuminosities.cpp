// Standard C++ libraries
#include <iomanip>

// Boost libraries
#include <boost/property_tree/ptree.hpp>
#include <boost/math/tools/roots.hpp>

// Standard user-made libraries
#include <fparameters/SpaceIterator.h>
#include <fparameters/parameters.h>
#include <fluminosities/thermalCompton.h>
#include <fluminosities/thermalSync.h>
#include <fluminosities/thermalBremss.h>
#include <fluminosities/luminosityHadronic.h>
#include <fluminosities/blackBody.h>
#include <fluminosities/reflection.h>
#include <fmath/mathFunctions.h>
#include <fmath/physics.h>
#include <fmath/fbisection.h>

// Project headers
#include "globalVariables.h"
#include "messages.h"
#include "read.h"
#include "thermalLuminosities.h"
#include "adafFunctions.h"
#include "write.h"

// Namespaces
using namespace std;

void localProcesses(const State& st, Matrix& lumOutSy, Matrix& lumOutBr, Matrix& lumOutpp, 
			Matrix scatt, Vector& energies, const int flags[], Matrix& lumOut)
{
	double lumSync1,lumSync2;
	lumSync1 = lumSync2 = 0.0;
	matrixInit4(lumOutSy,lumOutBr,lumOutpp,lumOut,nE,nR,0.0);
	
	size_t jE=0;
	st.photon.ps.iterate([&](const SpaceIterator& itE) {
		energies[jE] = itE.val(DIM_E);
		double frecuency=energies[jE]/planck;
		size_t jR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itER) {
			double r=itER.val(DIM_R);
			double thetaH = st.thetaH.get(itER);
			double rB1=r/sqrt(paso_r);
			double rB2=r*sqrt(paso_r);
			double area=2.0*pi*rB2*rB2*(2.0*cos(thetaH)+sin(thetaH)*sin(thetaH));
			double fluxToLum=area;
			double vol=(rB2*rB2*rB2-rB1*rB1*rB1)*(4.0/3.0)*pi*cos(thetaH);
			double emissToLum=vol*4.0*pi;
			double lumBr=0.0;
			double lumSy=0.0;
			double lumRJ=0.0;

			double temp=st.tempElectrons.get(itER);
			double temp_i=st.tempIons.get(itER);
			double magf=st.magf.get(itER);
			double dens_i=st.denf_i.get(itER);
			double dens_e=st.denf_e.get(itER);
			if (flags[0]) {
				double jSy = jSync(energies[jE],temp,magf,dens_e);
				lumRJ = bb_RJ(frecuency,temp) * fluxToLum;
				lumSy = jSy*emissToLum;
			}
			if (flags[1])
				lumBr = jBremss(energies[jE],temp,dens_i,dens_e)*emissToLum;
			if (flags[2] && temp_i > 1.0e11 && frecuency > 1.0e20 && frecuency < 1.0e26) {
				double jpp = luminosityHadronic(energies[jE],dens_i,temp_i);
				lumOutpp[jE][jR] = jpp*emissToLum;
			}
            if (flags[0]) {
				if (jR > 0)
					lumSync2 = min((1.0-scatt[jR-1][jR])*lumSync1+lumSy,lumRJ);
				else { lumSync2 = min(lumSy,lumRJ);}
				if (lumSync2>lumSync1)
					lumOutSy[jE][jR]=lumSync2-lumSync1;
				lumSync1=lumSync2;
            }
			if (flags[1]) {
				if (flags[0] && frecuency < 1.0e8)
					lumBr=min(lumBr,lumOutSy[jE][jR]);
				lumOutBr[jE][jR]=lumBr;
			}
			lumOut[jE][jR] = lumOutSy[jE][jR]+lumOutBr[jE][jR]+lumOutpp[jE][jR];
			jR++;
		},{itE.coord[DIM_E],-1,0});
		jE++;	
	},{-1,0,0});
}

/*double maxTemp(const State& st, const SpaceIterator& i)
{
	double temp1 = 0.0;
	double temp = 0.0;
	st.photon.ps.iterate([&](const SpaceIterator& itRT) {
		temp1 = st.tempElectrons.get(itRT);
		temp = (temp > temp1 ? temp : temp1);
	},{0,i.coord[DIM_R],-1});
	return temp;
}*/

void reflectedSpectrum(Matrix lumOut, Matrix absCD, Matrix& lumOutRefl, Vector energies,
						double pasoE)
{
	size_t numInt = 30;
	for (size_t jE=nE/2;jE<nE;jE++) {
		double freq = energies[jE]/planck;
		double x = energies[jE]/(electronMass*cLight2);
		double xMax = 1.0e2;
		for (size_t jRcd=0;jRcd<nRcd;jRcd++) {
			double freq0,pasoFreq,aux1;
			int logicalInt = 1;
			if (x > 0.03 && x < xMax) {
				freq0 = energies[jE]/planck;
				pasoFreq = pow(xMax/x,1.0/(numInt+1));
				aux1 = 0.0;
			} else if (x <= 0.03) {
				freq0 = 0.03*electronMass*cLight2/planck;
				pasoFreq = pow(xMax/0.03,1.0/(numInt+1));
				aux1 = greenDeltaFunc(x);
			} else { logicalInt = 0; aux1 = 0.0; }
			
			double lum = 0.0;
			if (logicalInt == 1) {
				for (size_t jjFreq=0;jjFreq<numInt;jjFreq++) {
					double dfreq0 = freq0*(pasoFreq-1.0);
					double lumInc = lumShellInc(freq0,jE,jRcd,lumOut,absCD,energies,nR);
					double green = greenFuncRefl(freq,freq0);
					lum += green*lumInc*dfreq0;//*(freq/freq0);
					freq0 *= pasoFreq;
				}
			}
			lumOutRefl[jE][jRcd] =
					lum + aux1*lumShellInc(freq,jE,jRcd,lumOut,absCD,energies,nR);
		}
	}
}

void coldDiskLuminosity(const State& st, Matrix absCD, Matrix lumOut, Matrix& lumOutRefl, 
						Matrix& lumOutCD, Vector energies)
{
	matrixInit(lumOutCD,nE,nRcd,0.0);
	matrixInit(lumOutRefl,nE,nRcd,0.0);
	
	size_t jRcd=0;
	st.photon.ps.iterate([&](const SpaceIterator& itRcd) {
		double rCd = itRcd.val(DIM_Rcd);
		double lj = rCd/sqrt(paso_rCD);
		double lj1 = rCd*sqrt(paso_rCD);
		double area = 2.0*pi*(lj1*lj1-lj*lj);
		double aux1 = auxCD(rCd);
		double aux2 = 0.0;
		
		double pasoE = pow(energies[nE-1]/energies[0],1.0/(nE+1));
		reflectedSpectrum(lumOut,absCD,lumOutRefl,energies,pasoE);

		for (size_t jjE=0;jjE<nE;jjE++) {
			double frequency = energies[jjE]/planck;
			double dfreq = frequency*(pasoE-1.0);
			double lumInc = 0.0;
			for (size_t kR=0;kR<nR;kR++) {
				lumInc += lumOut[jjE][kR]*absCD[kR][jRcd];
			}
			aux2 += (lumInc - lumOutRefl[jjE][jRcd])*dfreq;
		}
		double aux = aux1+aux2;
		double temp = pow(aux/area/ stefanBoltzmann,0.25);
		cout << temp/1.0e6 << endl;

		for(size_t jE=0;jE<nE;jE++) {
			double frecuency=energies[jE]/planck;
			lumOutCD[jE][jRcd] = area * bb(frecuency,temp);
		};

		jRcd++;
	},{0,0,-1});
}

void thermalCompton2(const State& st, Matrix& lumOut, Matrix scattADAF, Matrix scattCD, 
					Matrix absCD, Matrix& lumOutIC, Vector energies)
{
	show_message(msgStart,Module_thermalCompton);
	
	size_t nOm = 30;
	size_t nG = 10;
	Vector p(nR*nG*nE*nOm,0.0);
	
	size_t jR = 0;
	st.photon.ps.iterate([&](const SpaceIterator& itR) {
		double temp = st.tempElectrons.get(itR);
		double normtemp = boltzmann*temp/(electronMass*cLight2);
		if (normtemp >= 0.01) {
			comptonRedistribution(p,nG,nE,nOm,jR,normtemp,energies,lumOut);
		}
		jR++;
	},{0,-1,0});
	
	Matrix lumOutLocal,lumCD,lumOutRefl;
	Vector lumInIC(nE,0.0);
	matrixInit(lumOutIC,nE,nR,0.0);
	matrixInitCopy(lumOutLocal,nE,nR,lumOut);
			
	double res;
	size_t it=1;
	do {
		cout << "Iteration number = " << it << endl;
		res=0.0;

		fill(lumInIC.begin(),lumInIC.end(),0.0);
		matrixInit(lumOutIC,nE,nR,0.0);
		
		coldDiskLuminosity(st,absCD,lumOut,lumOutRefl,lumCD,energies);
		
		// PARA CADA CELDA
		size_t jR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itR) {
			// CALCULAMOS LinC PARA CADA E        			 function LinC
			for (size_t jE=0;jE<nE;jE++) {
				for (size_t jjR=0;jjR<nR;jjR++)
					lumInIC[jE] += scattADAF[jjR][jR]*lumOut[jE][jjR];
				for (size_t jjRcd=0;jjRcd<nRcd;jjRcd++)
					lumInIC[jE] += scattCD[jjRcd][jR]*lumCD[jE][jjRcd];
			}
			double temp = st.tempElectrons.get(itR);
			double normtemp = boltzmann*temp/(electronMass*cLight2);
			if (normtemp >= 0.01) 
			{
				for (size_t jjE=0;jjE<nE;jjE++) {
					double frecuency=energies[jjE]/planck;
					if (lumInIC[jjE]*frecuency > 1.0e15)
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


/*void gravRedshift(const State& st, Vector energies, Matrix lum, Matrix& lumRed)
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
}*/

/*void binEmissivities(const State& st, Vector esc, Vector energies, Matrix lumOut)
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
}*/

void writeLuminosities(State& st, Vector energies, Matrix lumOutSy, Matrix lumOutBr,
						Matrix lumOutpp, Matrix lumOutIC, Matrix lumOut, 
						Matrix lumOutCD, Matrix lumOutRefl,
						Vector esc, Vector escCD, const string& filename)
{
	ofstream file1,file2;
	file1.open(filename.c_str(),ios::out);
	file2.open("lumRadius.txt",ios::out);
	
	double lumSy,lumBr,lumIC,lumpp,lum,lumCD,lumRefl;
	for (size_t jE=0;jE<nE;jE++) {
		double frecuency = energies[jE]/planck;
		double energyEV = energies[jE]/1.6e-12;
		lumSy = lumBr = lumIC = lumpp = lum = lumCD = lumRefl = 0.0;
		for (size_t jR=0;jR<nR;jR++) {
			lumSy += lumOutSy[jE][jR] * esc[jR];
			lumBr += lumOutBr[jE][jR] * esc[jR];
			lumIC += lumOutIC[jE][jR] * esc[jR];
			lumpp += lumOutpp[jE][jR] * esc[jR];
			lum += lumOut[jE][jR] * esc[jR];
		}
		for (size_t jRcd=0;jRcd<nRcd;jRcd++) {
			lumCD += lumOutCD[jE][jRcd] * escCD[jRcd];
			lumRefl += lumOutRefl[jE][jRcd] * escCD[jRcd];
		}
        file1
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << energyEV
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumSy*frecuency
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumBr*frecuency
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumIC*frecuency
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumpp*frecuency
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lum*frecuency
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumCD*frecuency
			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumRefl*frecuency
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

/*void photonDensity(State& st, Vector energies, Matrix lumOut, Vector esc)
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
}*/

void thermalLuminosities(State& st, const string& filename, Matrix &scattADAF, 
							Matrix& scattCD, Matrix& absCD, Vector& esc, Vector& escCD)
{
	show_message(msgStart,Module_thermalLuminosities);

	Matrix lumOutSy,lumOutBr,lumOutpp,lumOutIC,lumOut,lumOutCD,lumOutRefl;
	matrixInit(lumOutIC,nE,nR,0.0);
	Vector energies(nE,0.0);

	int processesFlags[numProcesses];
	
	readThermalProcesses(processesFlags);
	if (processesFlags[0] || processesFlags[1] || processesFlags[2] || processesFlags[3]) {
		if (processesFlags[0] || processesFlags[1] || processesFlags[2])
			localProcesses(st,lumOutSy,lumOutBr,lumOutpp,scattADAF,energies,processesFlags,lumOut);
		if (processesFlags[4] && processesFlags[3])
			thermalCompton2(st,lumOut,scattADAF,scattCD,absCD,lumOutIC,energies);
		if (processesFlags[3])
			coldDiskLuminosity(st,absCD,lumOut,lumOutRefl,lumOutCD,energies);
		//gravRedshift(st,energies,lumOut,lumOutRed);
		//binEmissivities(st,esc,energies,lumOut);
		//photonDensity(st,energies,lumOut,esc);
		writeLuminosities(st,energies,lumOutSy,lumOutBr,lumOutpp,
							lumOutIC,lumOut,lumOutCD,lumOutRefl,esc,escCD,filename);
	}
	show_message(msgEnd,Module_thermalLuminosities);
}