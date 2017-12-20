#include "luminosities.h"

#include "modelParameters.h"
#include "write.h"

#include "functions.h"

#include <fparameters/SpaceIterator.h>
#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>

#include <fluminosities/thermalSync.h>
#include <fluminosities/thermalBremss.h>
#include <fluminosities/thermalIC.h>
#include <fmath/mathFunctions.h>
#include <fmath/physics.h>

#include <iomanip>
using namespace std;

void luminosities2(State& st, const std::string& filename) {
    
    std::ofstream file;
	file.open(filename.c_str(), std::ios::out);

	file << setw(10) << "frecuency [Hz]"
		//  << setw(10) << "Synchr"
        //<< setw(10) << "Brems"
		//  << setw(10) << "IC_Bremss"
        //  << setw(10) << "IC_Sync"
          << setw(10) << "Total"
          << std::endl;

    double nT = GlobalConfig.get<double>("model.particle.default.dim.theta.samples");
    double nR = GlobalConfig.get<double>("model.particle.default.dim.radius.samples");
    
    //nT = 10.0;
    //nR = 10.0;
    
    double thetaMin = 0.0;
    double thetaMax = pi / 4.0;
    double dtheta = (thetaMax - thetaMin) / nT;
    
    double rMin = GlobalConfig.get<double>("rCusp")*1.1;
    double rMax = GlobalConfig.get<double>("rCenter")*3.0;
    
    rMin = rMin * gravitationalConstant * 1.0e6 * solarMass / cLight2;
    rMax = rMax * gravitationalConstant * 1.0e6 * solarMass / cLight2;
    
    double dr = (rMax - rMin) / nR;

    st.photon.ps.iterate([&](const SpaceIterator& i) {
        
        cout << rMin << "\n" << rMax << "\n" << dr << "\n" << thetaMin << "\n" << thetaMax << "\n" << dtheta << std::endl;
        
        double energy = i.val(DIM_E);
        double frecuency = energy / planck;
        double luminosity = 0.0;

        st.photon.ps.iterate([&](const SpaceIterator& j) {
                
            // double fmtE = log10(energy / 1.6e-12);

            double r = j.val(DIM_R);
            
            r = r * gravitationalConstant * 1.0e6 * solarMass / cLight2;
            double theta = j.val(DIM_THETA);
        
            double vol = 2.0*pi*dtheta* r*r * dr;
        
            double norm_temp = boltzmann * st.tempElectrons.get(j) / electronMass / cLight2;
		
            double jBr = st.tpf1.get(j) * energy*energy;
            double jSy = st.tpf2.get(j);
                        
            double jIC1 = jIC_Bremss(energy, norm_temp, r, theta, j, st.denf_e, jBr, rMax);
            double jIC2 = jIC_Sync(energy, norm_temp, r, theta, j, st.denf_e, jSy, rMax);
        
            double jTotal =  jBr + jSy + jIC1 + jIC2;
        
            double emissToLum = 4.0 * pi * vol;
        
            //double luminosityBr = jBr * emissToLum;
            //double luminositySync = jSy * emissToLum;
            //double luminosityIC1 = jIC1 * emissToLum;
            //double luminosityIC2 = jIC2 * emissToLum;
            luminosity = luminosity +  jTotal * emissToLum;
            
                            
        }, { i.coord[DIM_E], -1, -1 } );
        
        file << "\t" << frecuency << "\t"
        //                        << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency * luminositySync
        //                        << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency * luminosityBr
        //                        << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency * luminosityIC1
        //                        << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency * luminosityIC2
                                    << frecuency * luminosity
                                    << std::endl;
    
    }, { -1, 0, 0 } );
    
     file.close();
    
}