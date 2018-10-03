#include "modelParameters.h"
#include "write.h"
#include "luminosities.h"
//#include "functions.h"
#include "messages.h"
#include <fparameters/SpaceIterator.h>
#include <fparameters/parameters.h>
#include <boost/property_tree/ptree.hpp>
#include <fluminosities/compton.h>
#include <fluminosities/thermalSync.h>
#include <fluminosities/thermalBremss.h>
#include <fluminosities/luminosityHadronic.h>
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

void luminosities(State& st, const string& filename, Matrix &prob) {
   
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
	Vector lumOutpp(nE,0.0);
	Vector energies(nE,0.0);
				
	Matrix lumOut1;
	Matrix lumOutSy;
    Matrix lumOutBr;
	matrixInit(lumOut1,nE,nR,0.0);
	matrixInit(lumOutSy,nE,nR,0.0);
    matrixInit(lumOutBr,nE,nR,0.0);

	int jE=0;
	double areaFactor=2.0*pi*rg*rg;
	double volFactor=areaFactor*rg/3.0;
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
			double lumBr=0.0;
			double lumSy=0.0;
			double lumRJ=0.0;
			double lumpp=0.0;
			st.photon.ps.iterate([&](const SpaceIterator& itERT) {
				double temp=st.tempElectrons.get(itERT);
				double magf=st.magf.get(itERT);
				double dens_i=st.denf_i.get(itERT);
				double dens_e=st.denf_e.get(itERT);
				lumBr += jBremss(energies[jE],temp,dens_i,dens_e)*emissToLum;
				double jSy=0.0;
				if (temp >= 1.0e8) {
					jSy = jSync(energies[jE],temp,magf,dens_e);
				}
				double jpp = luminosityHadronic(energies[jE],st.proton,dens_i,itERT);
				lumpp += jpp*emissToLum;
				lumSy += jSy*emissToLum;
				lumRJ += bb(frecuency,temp)*fluxToLum;
			},{itE.coord[DIM_E],itER.coord[DIM_R],-1});
            if (jR>0) {
				double a=1.0-prob[jR-1][jR];
				double b=a*lumSync1;
				double c=b+2.0*lumSy;
                lumSync2 = min(c,2.0*lumRJ);
				
            } else {
                lumSync2 = 2.0*min(lumSy,lumRJ);
            }
            lumBr *= 2.0;
			lumSy= max(lumSync2-lumSync1,0.0);
			//lumSy=lumSync2-lumSync1;
			lumSync1=lumSync2;
            if (lumBr > lumSy && frecuency < 1.0e8) lumBr=lumSy;
			lumOutSy[jE][jR]=lumSy;
            lumOutBr[jE][jR]=lumBr;
			lumOutpp[jE] += lumpp;
			lumOut1[jE][jR]=lumSy+lumBr+lumpp;
			lumSyTot[jE] += lumSy;
			lumBrTot[jE] += lumBr;
		jR++;
		},{itE.coord[DIM_E],-1,0});
	jE++;	
	},{-1,0,0});
	
	Vector lumICinCell(nE,0.0);
	Vector lumICinTot(nE,0.0);
	Vector lumICout_vec(nE,0.0);
	Matrix lumICout;
	Matrix lumOut;
	matrixInit(lumICout,nE,nR,0.0);
	matrixInit(lumOut,nE,nR,0.0);
	cout << "Starting IC" << endl;
	
	double res;
	int it=1;     // ITERACIONES PARA CONVERGENCIA DEL COMPTON
	do {
		cout << "Iteration number = " << it << endl;
		
		res=0.0;
		
		// SETEAMOS A 0 EN CADA ITERACIÓN
		if(it>1)
		{
			for (int jE=0;jE<nE;jE++) {
				lumICout_vec[jE]=0.0;
				lumICinTot[jE]=0.0;
				for (int jR=0;jR<nR;jR++) {
					lumICout[jE][jR]=0.0;
				}
			}
		}
		
		// PARA CADA CELDA
		int jR=0;
        st.photon.ps.iterate([&](const SpaceIterator& itR) {
			
			cout << "jR = " << jR << endl;
			// CALCULAMOS LinC PARA CADA E
            for (int jE=0;jE<nE;jE++) {
				double sumTot=0.0;
                for (int jjR=0;jjR<nR;jjR++) {
					sumTot += prob[jjR][jR]*lumOut1[jE][jjR];
                }
				lumICinCell[jE] = sumTot;
				lumICinTot[jE] += sumTot;
            }
			
			// CALCULAMOS LA MAYOR TEMPERATURA DE LA CELDA
            double temp=0.0;
            double temp1=0.0;
            st.photon.ps.iterate([&](const SpaceIterator& itRT) {
                temp1 = st.tempElectrons.get(itRT);
                temp = (temp > temp1 ? temp : temp1);
            },{0,itR.coord[DIM_R],-1});
			
			// CALCULAMOS LA COMPTONIZACION
			if (temp > 1.0e6) 
			{
				for (int jjE=0;jjE<nE;jjE++) {      // para cada energía del espectro inicial
													// calculamos la contribución a la luminosidad
													// total a cada energía del espectro final
													// y sumamos.
					double frecuency=energies[jjE]/planck;
					if (lumICinTot[jjE]*frecuency > 1.0e10) 
						compton2(lumICout,lumICinCell[jjE],energies,temp,nE,jjE,jR);
				}
			}
			
			for (int jE=0;jE<nE;jE++) {
				lumICout_vec[jE] += lumICout[jE][jR];
				double lumNew = lumOutSy[jE][jR] + lumOutBr[jE][jR] + lumICout[jE][jR];
				if (lumNew > 0.0 && lumOut1[jE][jR] > 0.0)
					res += P2(log10(lumNew) - log10(lumOut1[jE][jR]));
				lumOut[jE][jR] = lumNew;
			}
            
/*			if (temp > 1.0e6) {
			for (int jE=0;jE<nE;jE++) {
				lumICout[jE][jR] = compton(lumICinCell,energies,temp,nE,jE);
				lumICout_vec[jE] += lumICout[jE][jR];
				double lumNew = lumOutSy[jE][jR] + lumOutBr[jE][jR] + lumICout[jE][jR];
				if (lumNew > 0.0 && lumOut1[jE][jR] > 0.0)
					res += P2(log10(lumNew) - log10(lumOut1[jE][jR]));
				lumOut[jE][jR] = lumNew;
			}
			}*/
			
			jR++;
        },{0,-1,0});
		
		// COPIO LA MATRIZ A LA MATRIZ PREVIA
		for (int jE=0;jE<nE;jE++) {
			for (int jR=0;jR<nR;jR++) {
				lumOut1[jE][jR]=lumOut[jE][jR];
			}
		}
		
		res = sqrt(res)/(nE*nR);
		cout << "Residuo = " << res << endl;
	
	++it;	
	} while (res > 1.0e-3);
			
/*		st.photon.ps.iterate([&](const SpaceIterator& itR) {
			st.photon.ps.iterate([&](const SpaceIterator& itTh) {
				double dens_e=st.denf_e.get(itTh);
				double temp=st.tempElectrons.get(itTh);
					
				int jE=0;
				st.photon.ps.iterate([&](const SpaceIterator& itE) {
					double energy = itE.val(DIM_E);
					double frec=energy/planck;
				
					for (int jjR=0;jjR<nR;jjR++) {
						for (int jjE=0;jjE<nE;jjE++) {
							lumIC_in[jjE] += prob[jjR][jR]*lumOut[jR][jjE];
						}
					}
					double comp=compton(lumIC_in,energies,frec,temp,nE,energies[nE-1],energies[0]);
					lumIC_out=dens_e*comp;
					
					lumOut[jR][jE] = lumOut1[jR][jE] + lumIC_out;
					jE++;
				},{-1,itTh.coord[DIM_R],itTh.coord[DIM_THETA]}); //iterate en energia
			},{0,itR.coord[DIM_R],-1}); //iterate en theta
			jR++;
		},{0,-1,0}); ////iterate en r
	}
*/
	
	double distance = 100*rg;
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
		flux[jjE] = sum;// * 0.25 / pi / (distance*distance);
	}
	printf("lumtot = %5.5e\n",lum);
	
/*	Vector flux(nE,0.0);		
	jR=0;
	st.photon.ps.iterate([&](const SpaceIterator& itR) {                                 // loop en radios
		st.photon.ps.iterate([&](const SpaceIterator& itTh) {                           // loop en angulo polar
			double dens_e=st.denf_e.get(itTh);
			double temp=st.tempElectrons.get(itTh);
			Vector lumIC_in(nE,0.0);
			Vector lumIC_out(nE,0.0);

			int jE=0;
			st.photon.ps.iterate([&](const SpaceIterator& itE) {
				                    
				double energy=itE.val(DIM_E);
				double frec=energy/planck;
				for (int jjE=0;jjE<nE;jjE++) {
					for (int jjR=0;jjR<nR;jjR++) {
						lumIC_in[jjE]+=prob[jjR][jR]*lumOut[jjR][jjE];
					}
				}
				double comp = compton(lumIC_in,energies,frec,temp,nE,energies[nE-1],energies[0]);
				lumIC_out[jE]=dens_e*comp;
				lumOut[jR][jE]+=lumIC_out[jE];
				double aux = lumOut[jR][jE];
				flux[jE]+= aux;
				jE++;
			},{-1,itR.coord[DIM_R],itTh.coord[DIM_THETA]});
		},{0,itR.coord[DIM_R],-1});
		jR++;
	},{0,-1,0});
*/	
	for(jE=0;jE<nE;jE++) {
		double frecuency = energies[jE]/planck;
        file << "\t" << frecuency //energies[jE] / 1.6021e-12 / sqrt(-g_tt(zgap,pi/2.01))
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumSyTot[jE]*frecuency
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumBrTot[jE]*frecuency
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumICinTot[jE]*frecuency
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumICout_vec[jE]*frecuency
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumOutpp[jE]*frecuency
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << flux[jE]*frecuency
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << nph[jE]
			 << std::endl;
	};
		
	file.close();
	show_message(msgEnd, Module_luminosities);

}