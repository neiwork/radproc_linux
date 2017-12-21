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

    double rg = gravitationalConstant * 1.0e6 * solarMass / cLight2;
    
    double nT = GlobalConfig.get<double>("model.particle.default.dim.theta.samples");
    double nR = GlobalConfig.get<double>("model.particle.default.dim.radius.samples");
    
    double thetaMin = GlobalConfig.get<double>("thetamin");
    double thetaMax = GlobalConfig.get<double>("thetamax");
    double dtheta = (thetaMax - thetaMin) / nT;
    
    double rMin = GlobalConfig.get<double>("rmin") * rg;
    double rMax = GlobalConfig.get<double>("rmax") * rg;

    double dr = (rMax - rMin) / nR;

    st.photon.ps.iterate([&](const SpaceIterator& i) {
                
        double energy = i.val(DIM_E);
        
        double lumTot=0.0;
        double lumBr=0.0;
        double lumSy=0.0;
        double lumIC_Br=0.0;
        double lumIC_Sy=0.0;

        st.photon.ps.iterate([&](const SpaceIterator& j) {
                
            // double fmtE = log10(energy / 1.6e-12);

            double r = j.val(DIM_R) * rg;
            double theta = j.val(DIM_THETA);
        
            double vol = 2.0*pi*dtheta* r*r * dr;
        
            double norm_temp = boltzmann * st.tempElectrons.get(j) / electronMass / cLight2;
		
            double jBr = st.tpf1.get(j) * 0.25 * energy*energy / pi;
            double jSy = st.tpf2.get(j);
                        
            double jIC1 = jIC_Bremss(energy, norm_temp, r, theta, j, st.denf_e, jBr);
            double jIC2 = jIC_Sync(energy, norm_temp, r, theta, j, st.denf_e, jSy);
        
            double jTotal =  jBr + jSy + jIC1 + jIC2;
        
            double emissToLum = 4.0 * pi * vol;
        
            lumBr += jBr * emissToLum;
            lumSy += jSy * emissToLum;
            lumIC_Br = jIC1 * emissToLum;
            lumIC_Sy = jIC2 * emissToLum;
            lumTot +=  jTotal * emissToLum;
            
                            
        }, { i.coord[DIM_E], -1, -1 } );
        
        double frecuency = energy / planck;
        
        file << "\t" << frecuency << "\t"
                                << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency * lumBr
                                << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency * lumSy
                                << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency * lumIC_Br
                                << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency * lumIC_Sy
                                << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency * lumTot
                                << std::endl;
    
    }, { -1, 0, 0 } );
    
     file.close();
    
}