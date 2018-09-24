///#include <fluminosities/compton.h>
//#include <fluminosities/probMatrix.h>
#include "modelParameters.h"
#include "luminosities.h"
#include "messages.h"
#include "write.h"
#include <fparameters/SpaceIterator.h>
#include <fparameters/parameters.h>
#include <fparameters/nrutil.h>
#include <boost/property_tree/ptree.hpp>
#include <boost/math/special_functions/bessel.hpp>
#include <fluminosities/thermalSync.h>
#include <fluminosities/thermalBremss.h>
#include <fluminosities/blackBody.h>
#include <fmath/mathFunctions.h>
#include <fmath/physics.h>
#include <ssADAF/ssfunctions.h>
#include <boost/math/tools/roots.hpp>
#include <iomanip>

using namespace std;

void luminosities(State& st, const string& filename) {
   
	show_message(msgStart, Module_luminosities);
    std::ofstream file;
	file.open(filename.c_str(), ios::out);
	file  << setw(10) << "frecuency [Hz]"
          << setw(10) << "Bremss"
		  << setw(10) << "Synchr"
		  << setw(10) << "IC_Bremss"
          << setw(10) << "IC_Synchr"
          << std::endl;
	double rg=GlobalConfig.get<double>("rg");
    int nE=GlobalConfig.get<int>("model.particle.default.dim.energy.samples");
    int nR=GlobalConfig.get<int>("model.particle.default.dim.radius.samples");
	int nTheta=GlobalConfig.get<int>("mode.particle.default.dim.theta.samples");
    double rmax=GlobalConfig.get<double>("rMax");
    double rmin=GlobalConfig.get<double>("rMin");
	double thetamax=GlobalConfig.get<double>("thetaMax");
	double thetamin=GlobalConfig.get<double>("thetaMin");

	Vector energies(nE+1,0.0);
    double **lumOut,**prob;
    prob=dmatrix(0,nR,0,nR);
	torusSampling(prob,rmax,rmin,thetamax,thetamin,nR,nTheta);
    lumOut=dmatrix(0,nR,0,nE+1);
    int jR=0;
    st.photon.ps.iterate([&](const SpaceIterator& itR) {       // loop en radios
        r=itR.val(DIM_R);
        if (itR.its[DIM_R].canPeek(1) && itR.its[DIM_R].canPeek(-1)) {
			double rB1=sqrt(itR.its[DIM_R].peek(-1)*rad);
            double rB2=sqrt(itR.its[DIM_R].peek(1)*rad);
			double costhetaH=sqrt(pi/2.0*sqrdsoundvel())/(angularvel()*r*rg)*
								erff(angularvel()*r*rg/sqrt(2*sqrdsoundvel()));
			double sinthetaH2=1.0-costhetaH*costhetaH;
			double area=2.0*pi*rB2*rB2*rg*rg*(2.0*costhetaH+sinthetaH2);
			double vol=4.0*pi/3.0*(rB2*rB2*rB2-rB1*rB1*rB1)*rg*rg*rg*costhetaH;
					double fluxToLum=area;
					double emissToLum=4.0*pi*vol;
					double dens_e=st.denf_e.get(itTh);
					double dens_i=st.denf_i.get(itTh);
					double temp=st.tempElectrons.get(itTh);
					double magf=st.magfield.get(itTh);
					Vector lumBr(nE+1,0.0);
					Vector lumSync(nE+1,0.0);
					Vector lumSy(nE+1,0.0);
					Vector energies(nE+1,0.0);
					int jE=0;
					st.photon.ps.iterate([&](const SpaceIterator& itE) {
						energies[jE]=itE.val(DIM_E);
						double frec=energies[jE]/planck;
						double jBr=jBremss(energies[jE],temp,dens_i,dens_e);
						lumBr[jE]+=jBr*emissToLum;
						double jSy=(temp >= 5.0e8 ? jSync(energies[jE],temp,magf,dens_e) : 0.0);
						double auxSync=lumSync[jE];
						double lumRJ=bb_RJ(frec,temp)*fluxToLum;
						lumSync[jE]=min((1.0-prob[jR][jR+1])*jSy*emissToLum+auxSync,lumRJ);
						lumSy[jE]=lumSync[jE]-auxSync;
						lumOut[jR][jE]=lumBr[jE]+lumSy[jE];
						jE++;
					},{-1,itR.coord[DIM_R],itTh.coord[DIM_THETA]});
				},{0,itR.coord[DIM_R],-1});
				jR++;
			}
		},{0,-1,0});
    
		jR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itR) {                                 // loop en radios
			st.photon.ps.iterate([&](const SpaceIterator& itTh) {                           // loop en angulo polar
				double dens_e=st.denf_e.get(itTh);
				double temp=st.tempElectrons.get(itTh);
				Vector lumIC_in(nE+1,0.0);
				Vector lumIC_out(nE+1,0.0);
				Vector flux(nE+1,0.0);
				int jE=0;
				st.photon.ps.iterate([&](const SpaceIterator& itE) {
					double energy=itE.val(DIM_E);
					double frec=energy/planck;
					for (int jjE=0;jjE<nE;jjE++) {
						for (int jjR=0;jjR<nR;jjR++) {
							lumIC_in[jjE]+=prob[jjR][jR]*lumOut[jjR][jjE];
						}
					}
					lumIC_out[jE]=dens_e*compton(lumIC_in,energies,frec,temp,nE,energies[nE-1],energies[0]);
					lumOut[jR][jE]+=lumIC_out[jE];
					flux[jE]+=lumOut[jR][jE];
					jE++;
				},{-1,itR.coord[DIM_R],itTh.coord[DIM_THETA]});
			},{0,itR.coord[DIM_R],-1});
			jR++;
		},{0,-1,0});
	}
}