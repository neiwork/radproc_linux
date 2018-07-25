
//#include "probCompton.h"
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
	file << setw(10) << "frecuency [Hz]"
          << setw(10) << "Bremss"
		  << setw(10) << "Synchr"
		  << setw(10) << "IC_Bremss"
          << setw(10) << "IC_Synchr"
          << std::endl;
	double rg=GlobalConfig.get<double>("rg");
    
	int nE=GlobalConfig.get<int>("model.particle.default.dim.energy.samples");
    int nR=GlobalConfig.get<int>("model.particle.default.dim.radius.samples");
    //double rmax=GlobalConfig.get<double>("rEdge");
    //double rmin=GlobalConfig.get<double>("rCusp");
	
	
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
     
    // RECORRIENDO PRIMERO LOS RADIOS

	Vector lumBr(nE+1,0.0);
	Vector lumSync(nE+1,0.0);
	Vector lumSy(nE+1,0.0);
	Vector energies(nE+1,0.0);
				
    double **lumOut;//,**prob;
    //prob=dmatrix(0,nR+1,0,nR+1);
    lumOut=dmatrix(0,nR+1,0,nE+1);
    int jR=0;
   // probCompton(prob,nR);
    st.photon.ps.iterate([&](const SpaceIterator& itR) {                                // loop en radios
        double rB1=itR.val(DIM_R);
        if (itR.its[DIM_R].canPeek(1)) {
            double rB2=itR.its[DIM_R].peek(1);
            st.photon.ps.iterate([&](const SpaceIterator& itTh) {                           // loop en angulo polar
                double theta=itTh.val(DIM_THETA);
                double dth=itTh.its[DIM_THETA].peek(1)-theta;
                double area=2.0*pi*rB2*rB2*dth*rg*rg;
                double vol=2.0*pi/3.0*dth*(rB2*rB2*rB2-rB1*rB1*rB1)*rg*rg*rg;
                double fluxToLum=area;
                double emissToLum=4.0*pi*vol;
                double dens_e=st.denf_e.get(itTh);
                double dens_i=st.denf_i.get(itTh);
                double temp=st.tempElectrons.get(itTh);
                double magf=st.magf.get(itTh);
               
                //Vector energies(nE+1,0.0);
                int jE=0;
                st.photon.ps.iterate([&](const SpaceIterator& itE) {

					double energy = itE.val(DIM_E);
					energies[jE] = energy;
                    double frec=energy/planck;
                    double jBr=jBremss(energy,temp,dens_i,dens_e);
                    lumBr[jE]+=jBr*emissToLum;
                    double jSy=(temp >= 5.0e8 ? jSync(energy,temp,magf,dens_e) : 0.0);
                    double auxSync=lumSync[jE];
                    double lumRJ=bb_RJ(frec,temp)*fluxToLum;
                    lumSync[jE]=min((1.0-prob[jR][jR+1])*jSy*emissToLum+auxSync,lumRJ);
                    lumSy[jE]=lumSync[jE]-auxSync; //esto da 0 ?
                    lumOut[jR][jE]=lumBr[jE]+lumSy[jE];
                    jE++;
                },{-1,itR.coord[DIM_R],itTh.coord[DIM_THETA]}); //iterate en energia
            },{0,itR.coord[DIM_R],-1}); //iterate en theta
            jR++;
        }
    },{0,-1,0}); ////iterate en r
    
	
	
	Vector flux(nE+1,0.0);
    jR=0;
    st.photon.ps.iterate([&](const SpaceIterator& itR) {                                 // loop en radios
        st.photon.ps.iterate([&](const SpaceIterator& itTh) {                           // loop en angulo polar
            double dens_e=st.denf_e.get(itTh);
            double temp=st.tempElectrons.get(itTh);
            Vector lumIC_in(nE+1,0.0);
            Vector lumIC_out(nE+1,0.0);

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
	
	
	int jE = 0;
	st.photon.ps.iterate([&](const SpaceIterator& itE) {
	 //flux *= 0.25 / (pi * distance*distance);
        file << "\t" << log10(itE.val(DIM_E)/1.6e-12)
                           << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumBr[jE]
						   << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumSync[jE]
						   << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << flux[jE]
                           << std::endl;
		jE++;
	},{-1,0,0});
	
	
	file.close();
	show_message(msgEnd, Module_luminosities);

}