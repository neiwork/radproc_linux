#include "luminosities.h"

#include "thermalCompton.h"

#include "modelParameters.h"
#include "write.h"
#include <fmath/bisection.h>

#include "functions.h"

#include <fparameters/SpaceIterator.h>
#include <fparameters/parameters.h>

#include <boost/property_tree/ptree.hpp>
#include <boost/math/special_functions/bessel.hpp>

#include <fluminosities/thermalSync.h>
#include <fluminosities/thermalBremss.h>
#include <fluminosities/thermalIC.h>
#include <fluminosities/blackBody.h>
#include <fmath/mathFunctions.h>
#include <fmath/physics.h>
#include <boost/math/tools/roots.hpp>

#include <iomanip>
using namespace std;

/*
template <class T = double>
struct aux_function {
    aux_function(T const& a, T const& b, T const& c, T const& d) : 
                                        r(a), denf(b), magfield(c), norm_temp(d) {}
    T operator()(T const& xc) {
        T fx = exp(1.8899* pow(xc, 1.0/3) ) * pow(xc, 5.0/3) - 3.13e-9*denf*r / magfield / 
            (norm_temp*norm_temp*norm_temp) / boost::math::cyl_bessel_k(2, 1.0/norm_temp) *
            ( sqrt(xc) + 0.4/pow(xc,1.0/4) + 0.5316 );
        return fx;
    }
private:
    T r;
    T denf;
    T magfield;
    T norm_temp;
};

template <class T = double>
T xc_root(T r, T a, T b, T c) {   
    using namespace std;
    using namespace boost::math::tools;

   //T guess = 1.e-20;
   //T factor = 1.1;                     // How big steps to take when searching.
   
   T min = 1.e-5;
   T max = 1.e3;

   const boost::uintmax_t maxit = 5000;  // Limit to maximum iterations.
   boost::uintmax_t it = maxit;             // Initally our chosen max iterations, but updated with actual.
   // bool is_rising = true;
   // Define a termination condition, stop when nearly all digits are correct, but allow for
   // the fact that we are returning a range, and must have some inaccuracy in the elliptic integral:
   eps_tolerance<T> tol(std::numeric_limits<T>::digits - 2);   
   // Call bracket_and_solve_root to find the solution, note that this is a rising function:
   //pair<T, T> x = bracket_and_solve_root(aux_function<T>(r, a, b, c), guess, factor, is_rising, tol, it);
   pair<T,T> x = bisect(aux_function<T>(r,a,b,c), min, max, tol, it);
   // Result is midway between the endpoints of the range:
   return x.first + (x.second - x.first) / 2;
}
*/

void luminosities2(State& st, const string& filename) {
   

    ofstream file;
	file.open(filename.c_str(), ios::out);

	file << setw(10) << "frecuency [Hz]"
          << setw(10) << "Brems"
		  << setw(10) << "Synchr"
		  << setw(10) << "BB"
          << setw(10) << "BB_RJ"
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

  /*  st.photon.ps.iterate([&](const SpaceIterator& i) {
                
        double energy = i.val(DIM_E);
        
        double lumBr=0.0;
        double lumSy=0.0;
        double lumBB=0.0;
        //double lumIC_Br=0.0;
        //double lumIC_Sy=0.0;
        //double lum_icBr = 0.0;
        //double lum_icSy = 0.0;

        st.photon.ps.iterate([&](const SpaceIterator& j) {

            double r = j.val(DIM_R) * rg;
            double theta = j.val(DIM_THETA);
        
            double vol = 2.0 * 2.0*pi * dtheta * (r*r*dr) + dr*dr*dr/12.0;
            double area = 2.0 * 2.0*pi * dtheta * (r+dr/2.0)*(r+dr/2.0);
			
			double temp = st.tempElectrons.get(j);
            double norm_temp = boltzmann * temp / electronMass / cLight2;
		
            double jBr = st.tpf1.get(j) * 0.25 * energy*energy / pi;
            double jSy = st.tpf2.get(j)* 0.25 * energy*energy / pi;
			
			double fluxBB_RJ = bb(frecuency, temp);
            
            double emissToLum = 4.0 * pi * vol;
            double fluxToLum = area;
            double rf = redshiftFactor(r/rg, theta);
            
            //if (fluxBB_RJ > jSy *4.0 * pi * r / 3.0) {
                //jSy *= rf;
                //lumSy += (jSy * emissToLum);
            //} else {
                //fluxBB_RJ *= rf;
                //lumSy += (fluxBB_RJ * fluxToLum);
            //}
            jSy *= rf;
            lumSy += (jSy*emissToLum);
            
            fluxBB_RJ *= rf;
            lumBB += (fluxBB_RJ * fluxToLum);
            
            //double jIC1(0.0), jIC2(0.0);
			
            //double icBr = thCompton(energy, st.electron, j, st.tpf1, st.photon.emin());
            //double icSy = thCompton(energy, st.electron, j, st.tpf2, st.photon.emin());
            
            //const double magfield{ st.magf.get(j) };
            //const double denf{ st.denf_e.get(j) };
            //double xc = xc_root(r, denf, magfield, norm_temp);
                
            //jIC1 = jIC_Bremss(energy, norm_temp, r, theta, j, st.denf_e, jBr, xc);
            //jIC2 = jIC_Sync(energy, norm_temp, r, theta, j, st.denf_e, jSy, xc);
            
            //jIC1 *= rf;
            //jIC2 *= rf;

            jBr *= rf;
            lumBr += (jBr * emissToLum);
            //lumIC_Br += jIC1 * emissToLum;
            //lumIC_Sy += jIC2 * emissToLum;
            //lumTot +=  jTotal * emissToLum;
            //lum_icBr += icBr * emissToLum;
            //lum_icSy += icSy * emissToLum;
            
                            
        }, { i.coord[DIM_E], -1, -1 } );
        
        double frecuency = energy / planck;
        
        file      //<< "\t" << log10(energy/1.6e-12) 
				   << "\t" << log10(frecuency)
                                << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency * lumBr
                                << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency * lumSy
                                << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency * lumBB
                  //          << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency * lumIC_Br
                  //          << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency * lumIC_Sy
                  //          << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency * lumTot
		          //			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lum_icBr
		          //			<< setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lum_icSy
                                << std::endl;
    
    }, { -1, 0, 0 } ); */
    
//INTERCAMBIANDO EL ORDEN DE LAS ITERACIONES

    int energyDim = 50;
    vector<double> frecuency(energyDim, 0.0);
    vector<double> lumSy(energyDim, 0.0);
    vector<double> lumBr(energyDim, 0.0);
    vector<double> lumBB(energyDim, 0.0);
    vector<double> lumBB_RJ(energyDim, 0.0);
    
   st.photon.ps.iterate([&](const SpaceIterator& i) {
        
        double r = i.val(DIM_R) * rg;
        double theta = i.val(DIM_THETA);
        
        double vol = 2.0 * 2.0*pi * dtheta * (r*r*dr + dr*dr*dr/12.0);
        double area = 2.0 * 2.0*pi * dtheta * (r+dr/2.0)*(r+dr/2.0);
        
        double emissToLum = 4.0 * pi * vol;
        double fluxToLum = area;
        
        //double norm_temp = boltzmann * st.tempElectrons.get(i) / electronMass / cLight2;
        
        //const double magfield{ st.magf.get(i) };
        //const double denf{ st.denf_e.get(i) };
        const double temp{ st.tempElectrons.get(i) };
        //double xc = xc_root(r, denf, magfield, norm_temp);
        
        double rf = redshiftFactor(r/rg, theta);

        int k=0;
        st.photon.ps.iterate([&](const SpaceIterator& j) {
            
            double energy = j.val(DIM_E);
            frecuency[k] = energy / planck;
            double jBr = rf * st.tpf1.get(j) * 0.25 * energy*energy / pi;
            double jSy = rf * st.tpf2.get(j)* 0.25 * energy*energy / pi;
            double fluxBB_RJ = rf * bb_RJ(frecuency[k], temp);
            
            if (fluxBB_RJ > jSy *4.0 * pi * r / 3.0) {
                lumSy[k] += (jSy * emissToLum);
            } else {
                lumSy[k] += (fluxBB_RJ * fluxToLum);
            }
           
            lumBr[k] += jBr*emissToLum;

            k++;
                            
        }, { -1, i.coord[DIM_R], i.coord[DIM_THETA] } );
    }, { 0, -1, -1 } );
    
    for (int k=0; k<energyDim; k++) {
        file << "\t" << log10(frecuency[k]) 
                           << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency[k] * lumBr[k]
                           << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency[k] * lumSy[k]
                           << std::endl;
    };
	
	// PRUEBAS PARA SYNCROTHRON //////////////////////////
	/*
	std::ofstream file2;
	file2.open("sync.txt", std::ios::out);

	file2 << setw(10) << "frecuency [Hz]"
		  << setw(10) << "Synchr"
		  << setw(10) << "BB_RJ"
          << setw(10) << "BB"
          << std::endl;
	
	st.photon.ps.iterate([&](const SpaceIterator& i) {
		
		double rg = gravitationalConstant * 1.0e6 * solarMass / cLight2;
		double r = i.val(DIM_R) * rg;
		double temp = st.tempElectrons.get(i);
		
		st.photon.ps.iterate([&](const SpaceIterator& j) {
			double energy = j.val(DIM_E);
			double frecuency = energy / planck;
			
            double jSy = st.tpf2.get(j)* 0.25 * energy*energy / pi;
			
			double fluxBB_RJ = bb_RJ(frecuency, temp);
			double fluxBB = bb(frecuency, temp);
			
			double lumSy = jSy * 4.0 *pi / 3.0 * r*r*r;
			double lumBB_RJ = fluxBB_RJ * 4.0 * pi * r*r;
			double lumBB = fluxBB * 4.0 * pi * r*r;
			
			
			file2 << "\t" << log10(frecuency)
						  << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumSy
						  << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumBB_RJ
						  << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << lumBB
                          << std::endl;
			
			
		}, {-1, i.coord[DIM_R], i.coord[DIM_THETA]} );
	}, {0,-1,-1} );
	
    file2.close(); */
	
	///////////////////////////////////////////////////////////
	
	file.close();
    
}