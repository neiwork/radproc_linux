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
	
    /* int nT = GlobalConfig.get<int>("model.particle.default.dim.theta.samples");
    int nR = GlobalConfig.get<int>("model.particle.default.dim.radius.samples");
    int nEnergy = GlobalConfig.get<int>("model.particle.default.dim.energy.samples");
    double thetaMin = GlobalConfig.get<double>("thetamin");
    double thetaMax = GlobalConfig.get<double>("thetamax");
    double rMin = GlobalConfig.get<double>("rmin");
    double rMax = GlobalConfig.get<double>("rmax");
    rBoundary[0] = rMin;
    for (size_t i=1; i <= nR; i++) {
        rBoundary[i]=rBoundary[i-1]+dr;
        r[i]=rBoundary[i-1]+dr/2.0;
    } */
    /*
    double prob[nR][nR];
    

        flux *= 0.25 / (pi * distance*distance);
        file << "\t" << log10(frecuency)
                           << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency * flux
                           << std::endl;
    }, {-1, 0, 0} );

     */

	double lumBr;
	double lumSync1, lumSync2;
	double lumSy;
	Vector lumSyTot(nE,0.0);
	Vector lumBrTot(nE,0.0);
	Vector energies(nE,0.0);
				
    //double **lumOut;//,**prob;
	Matrix lumOut;
	Matrix lumOutSy;
    Matrix lumOutBr;
	matrixInit(lumOut,nE,nR,0.0);
	matrixInit(lumOutSy,nE,nR,0.0);
    matrixInit(lumOutBr,nE,nR,0.0);
    //matrixInit(prob,nR,nR,0.0);
    //matrixRead(prob,nR,nR);

	int jE=0;
	st.photon.ps.iterate([&](const SpaceIterator& it1) { // LOOP EN ENERGÍAS
		energies[jE]=it1.val(DIM_E);
		double frecuency=energies[jE]/planck;
		int jR=0;
		st.photon.ps.iterate([&](const SpaceIterator& it2) { // LOOP EN RADIOS
			double r=it2.val(DIM_R);
			double rB1=r-dr/2.0;
			double rB2=r+dr/2.0;
			double area=2.0*pi*rB2*rB2*rg*rg*dtheta;
			double fluxToLum=area;
			double vol=2.0*pi/3.0 * (rB2*rB2*rB2-rB1*rB1*rB1)*rg*rg*rg * dtheta;
			double emissToLum=vol*4.0*pi;
			double lumBr=0.0;
			double lumSy=0.0;
			double lumRJ=0.0;
			st.photon.ps.iterate([&](const SpaceIterator& it3) {
				double temp=st.tempElectrons.get(it3);
				double magf=st.magf.get(it3);
				double dens_i=st.denf_i.get(it3);
				double dens_e=st.denf_e.get(it3);
				lumBr += jBremss(energies[jE],temp,dens_i,dens_e)*emissToLum;
				double jSy;
				if (temp >= 5.0e7) {
					jSy = jSync(energies[jE],temp,magf,dens_e);
				} else {
					jSy = 0.0;
				}
				lumSy += jSy*emissToLum;
				lumRJ += bb_RJ(frecuency,temp)*fluxToLum;
			},{it1.coord[DIM_E],it2.coord[DIM_R],-1});
            if (jR>0) {
                lumSync2 = min((1.0-prob[jR-1][jR])*lumSync1+2.0*lumSy,2.0*lumRJ);
            } else {
                lumSync2 = lumSy;
            }
            lumBr *= 2.0;
			lumSy= max(lumSync2-lumSync1,0.0);
            //lumSy=lumSync2-lumSync1;
			lumSync1=lumSync2;
            if (lumBr > lumSy && frecuency < 1.0e8) lumBr=lumSy;
			lumOutSy[jE][jR]=lumSy;
            lumOutBr[jE][jR]=lumBr;
			lumOut[jE][jR]=lumOutSy[jE][jR]+lumOutBr[jE][jR];
			lumSyTot[jE] += lumSy;
			lumBrTot[jE] += lumBr;
		jR++;
		},{it1.coord[DIM_E],-1,0});
	jE++;	
	},{-1,0,0});

/*	int jR=0;
	st.photon.ps.iterate([&](const SpaceIterator& itR) {      // loop en radios
		double rB1=itR.val(DIM_R);
		if (itR.its[DIM_R].canPeek(1)) {
			double rB2=itR.its[DIM_R].peek(1);
			double area=2.0*pi*rB2*rB2*rg*rg;
			double vol=2.0*pi/3.0*(rB2*rB2*rB2-rB1*rB1*rB1)*rg*rg*rg;
			st.photon.ps.iterate([&](const SpaceIterator& itTh) {    // loop en angulo polar
				if (itTh.its[DIM_THETA].canPeek(1)) {
					double theta=itTh.val(DIM_THETA);
					double dth=abs(itTh.its[DIM_THETA].peek(1)-theta);
					area *= dth;
					vol *= dth;
					double fluxToLum=area;
					double emissToLum=4.0*pi*vol;
					double dens_e=st.denf_e.get(itTh);
					double dens_i=st.denf_i.get(itTh);
					double temp=st.tempElectrons.get(itTh);
					double magf=st.magf.get(itTh);
			
					int jE=0;
					st.photon.ps.iterate([&](const SpaceIterator& itE) {
						double energy = itE.val(DIM_E);
						if ( jR==1 ) {
							energies[jE] = energy;
						}
						double frec=energy/planck;
						double jBr=jBremss(energy,temp,dens_i,dens_e);
						lumBr = jBr*emissToLum;
						double jSy = (temp >= 5.0e8 ? jSync(energy,temp,magf,dens_e) : 0.0);
						auxSync=lumSync;
						double lumRJ=bb_RJ(frec,temp)*fluxToLum;
						lumSync=min((1.0-prob[jR][jR+1])*auxSync+jSy*emissToLum,lumRJ);
						lumSy = max(lumSync-auxSync,0.0);
						//
						lumOut1[jR][jE] += lumSy;
						lumOut[jR][jE] = lumOut1[jR][jE];
						jE++;
					},{-1,itTh.coord[DIM_R],itTh.coord[DIM_THETA]}); //iterate en energia
				}
			},{0,itR.coord[DIM_R],-1}); //iterate en theta
			jR++;
		}
	},{0,-1,0}); ////iterate en r
*/	
	
    Vector lumICinBr(nE,0.0);
	double lumICoutSy, lumICoutBr;
	Vector lumICoutvecSy(nE,0.0);
    Vector lumICoutvecBr(nE,0.0);
	
	for (int it=1;it<=2;it++) {                    // ITERACIONES PARA CONVERGENCIA DEL COMPTON
        int jR=0;
        Vector lumICinSy(nE,0.0);
        st.photon.ps.iterate([&](const SpaceIterator& it1) {
            for (int jjE=0;jjE<nE;jjE++) {
                double sumSy=0.0;
                double sumBr=0.0;
                double sumIC=0.0;
                for (int jjR=0;jjR<nR;jjR++) {
                    sumSy += prob[jjR][jR]*lumOut[jjE][jjR];
                    //sumBr += prob[jjR][jR]*lumOutBr[jjE][jjR];
                }
                lumICinSy[jjE] += sumSy;
                //lumICinBr[jjE] += sumBr;
            }
            double temp=0.0;
            double temp1=0.0;
            st.photon.ps.iterate([&](const SpaceIterator& it2) {
                temp1 = st.tempElectrons.get(it2);
                temp = (temp > temp1 ? temp : temp1);
            },{0,it1.coord[DIM_R],-1});
/*            for (int jE=0;jE<nE;jE++) {
                double om=energies[jE]/(electronMass*cLight2);
                lumICout1 = compton(lumICin1,energies,temp,nE,jE);
                lumICoutvec1[jE] += lumICout1;
                lumICout2 = compton(lumICin2,energies,temp,nE,jE);
                lumICoutvec2[jE] += lumICout2;
                lumOut[jE][jR] += lumICout1+lumICout2;
            } */
            
            if ( temp > 1.0e3) {
                for (int jE=0;jE<nE;jE++) {
                    if ( lumICinSy[jE] > 1.0e10 ) compton2(lumICoutvecSy,lumICinSy[jE],energies,temp,nE,jE);
                    //if (lumICinBr[jE] > 1.0e10) compton2(lumICoutvecBr,lumICinBr[jE],energies,temp,nE,jE);
                    lumOut[jE][jR] += lumICoutvecSy[jE];//+lumICoutvecBr[jE];
                }
            }
                
            jR++;
        },{0,-1,0});
    }
			
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
	
	double distance = 100.0*rg;
	Vector flux(nE,0.0);
    Vector nph(nE,0.0);
    double lum=0.0;
	for (int jjE=0;jjE<nE-1;jjE++) {
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
		flux[jjE] = sum * 0.25 / pi / (distance*distance);
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
//			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumICinSy[jE]*frecuency
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumICinBr[jE]*frecuency
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumICoutvecSy[jE]*frecuency
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumICoutvecBr[jE]*frecuency
			 << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency*nph[jE]
			 << std::endl;
	};
		
	file.close();
	show_message(msgEnd, Module_luminosities);

}