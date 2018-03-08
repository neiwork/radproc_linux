#include "luminosities.h"

#include "thermalCompton.h"

#include "modelParameters.h"
#include "write.h"
#include <fmath/bisection.h>

#include "functions.h"

#include "messages.h"

//#include "coppiBlandford.h"
//#include "photonEscape.h"

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
   
	show_message(msgStart, Module_luminosities);
	
    ofstream file;
	file.open(filename.c_str(), ios::out);

	file << setw(10) << "frecuency [Hz]"
          << setw(10) << "Bremss"
		  << setw(10) << "Synchr"
		  << setw(10) << "IC_Bremss"
          << setw(10) << "IC_Synchr"
          << std::endl;

	double rg = GlobalConfig.get<double>("rg");
    
    //double nT = GlobalConfig.get<double>("model.particle.default.dim.theta.samples");
    double nR = GlobalConfig.get<double>("model.particle.default.dim.radius.samples");
    //int nEnergy = GlobalConfig.get<int>("model.particle.default.dim.energy.samples");
    
    //double thetaMin = GlobalConfig.get<double>("thetamin");
    //double thetaMax = GlobalConfig.get<double>("thetamax");
    
    double rMin = GlobalConfig.get<double>("rmin");
    double rMax = GlobalConfig.get<double>("rmax");
    
    /////////////////////////////////////
    // VECTOR DE FRONTERAS DE LAS CELDAS //
    /////////////////////////////////////
    
    vector<double> rBoundary(nR+1, 0.0);
    double rVar = pow(rMax/rMin, 1.0/(rBoundary.size()-1.0));
    rBoundary[0] = rMin;
    for (size_t i=1; i < rBoundary.size(); i++) {
        rBoundary[i] = rBoundary[i-1]*rVar;
    }
    /////////////////////////////////////
    /////////////////////////////////////
    
    st.photon.ps.iterate([&](const SpaceIterator& k0) {   // loop en las energías
        
        double energy = k0.val(DIM_E);
        double frecuency = energy / planck;
        
        vector<double> lumOut(nR+1, 0.0);
        vector<double> lumBr(nR+1, 0.0);
        vector<double> lumSy(nR+1, 0.0);
        vector<double> lumSync(nR+1, 0.0);
        vector<double> lumIC_in(nR+1, 0.0);
        vector<double> lumIC_out(nR+1, 0.0);
        double flux = 0.0;
        
        int jR = 1;
        st.photon.ps.iterate([&](const SpaceIterator& k1) {   // loop en radios
            
            double l0 = rBoundary[jR-1] * rg;
            double l1 = rBoundary[jR] * rg;
            
            double lumRJ = 0.0;
            
            st.photon.ps.iterate([&](const SpaceIterator& k2) {   // loop en ángulo theta
                
                double theta = k2.val(DIM_THETA);
                double dtheta = k2.its[DIM_THETA].peek(1) - theta;
                
                double temp = st.tempElectrons.get(k2);
                double dens_i = st.denf_i.get(k2);
                double dens_e = st.denf_e.get(k2);
                double magfield = st.magf.get(k2);
                
                // TENGO QUE VER COMO ARMAR LOS VOLÚMENES
                double area = 2.0 * 2.0*pi * l1*l1 * dtheta;
                double vol = 2.0 * 2.0*pi / 3.0 * dtheta * (l1*l1*l1 - l0*l0*l0);
                double fluxToLum = area;
                double emissToLum = 4.0*pi*vol;
                /////////////////////////
                
                double jBr = jBremss(energy, temp, dens_i, dens_e);
                
                double jSy;
                if (temp >= 5.0e8) {
                    jSy = jSync(energy, temp, magfield, dens_e);
                } else { jSy = 0.0; };
                
                lumBr[jR] += jBr * emissToLum;
                lumSync[jR] += jSy * emissToLum;     // temporalmente
                lumRJ += bb_RJ(frecuency, temp) * fluxToLum;
                
            }, {0, k1.coord[DIM_R], -1} );
            
            //lumSync[jR] = min( (1.0- prob[jR, jR+1])*lumSync[jR-1] + lumSync[jR], lumRJ);
            lumSync[jR] = min(lumSync[jR-1] + lumSync[jR], lumRJ);
            lumSy[jR] = lumSync[jR] - lumSync[jR-1];
            
            //lumIC_out[j] = compton(lumIC_in);
            lumOut[jR] = lumBr[jR] + lumSy[jR] + lumIC_out[jR];
            
      //      for (int kR=0; kR<nR; kR++) {
       //         lumIC_in[jR] += prob[kR, jR] * lumOut[kR];
       //     }
            
            
            //flux += escape[jR]*lumOut[jR]; 
            flux += lumOut[jR];
            
            jR++;
            
        }, {k0.coord[DIM_E], -1, 0} );
         
        double distance = 1.0; 
        flux *= 0.25 / (pi * distance*distance);
        
        file << "\t" << log10(frecuency)
                           << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency * flux
                           << std::endl;
        
    }, {-1, 0, 0} );
    
 //   double ftheta = pow(thetaMax/thetaMin, 1.0 / double(nT));
 //  double fr = pow(rMax/rMin, 1.0 / double(nR));

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
    
//vector<double> r(nR-1, 0.0);
//vector<double> thetamax(nR-1, 0.0);

/*
int h = 0;
double l0, l3;
l3 = rMin;

st.photon.ps.iterate([&](const SpaceIterator& i) {
    
    if (i.its[DIM_R].canPeek(1)) {
        l0 = l3;
        l3 = i.its[DIM_R].peek(1);
        r[h] = sqrt( l0 * l3 );
        double radius = r[h];
        thetamax[h] = bisection(0.0, thetaMax, [&radius](double theta){return w(radius, theta);});
        h++;
    }
    
}, {0, -1, 0} ); */

/*    
//INTERCAMBIANDO EL ORDEN DE LAS ITERACIONES

    vector<double> frecuency(nEnergy, 0.0);
    vector<double> frecuency2(nEnergy, 0.0);
    vector<double> lumSy(nEnergy, 0.0);
    vector<double> lumBr(nEnergy, 0.0);
    vector<double> lumBB(nEnergy, 0.0);
    vector<double> lumBB_RJ(nEnergy, 0.0);
    vector<double> lumICB(nEnergy, 0.0);
    vector<double> lumICS(nEnergy, 0.0);
	
   st.photon.ps.iterate([&](const SpaceIterator& i) {

        double l1 = i.val(DIM_R) * rg;
        double l2 = i.its[DIM_R].peek(1) * rg;
        // Acá lo que quiero es seleccionar el punto intermedio entre los dos próximos valores de r (equiespaciado
        // logaritmicamente, de manera que r1 = sqrt(l1*l2). Si promediamos en theta para compton sólo habrá que
        // recorrer el espacio de parámetros en r
        
        double r = sqrt(l1*l2);
        double theta = i.val(DIM_THETA);
        
        double dtheta = theta*(ftheta-1.0);
        //double dr1 = r*(1.0-1.0/fr);
        double dr1 = r-l1;
        //double dr2 = r*(fr-1.0);
        double dr2 = l2-r;
        double vol = 2.0 * 2.0*pi * dtheta * 0.5*( r*r*(dr1+dr2) + 0.5*r*(dr2*dr2-dr1*dr1) +
                                                                                                                (dr1*dr1*dr1 + dr2*dr2*dr2)/12.0 );
        //double area = 2.0 * 2.0*pi * dtheta * (l2/2.0)*(l2/2.0);
        
        double emissToLum = 4.0 * pi * vol;
        //double fluxToLum = area;
        
        const double dens_e{ st.denf_e.get(i) };
        //const double dens_i{ st.denf_i.get(i) };
        const double magfield{ st.magf.get(i) };
        const double temp{ st.tempElectrons.get(i) };
        const double normTemp = boltzmann * temp / electronMass / cLight2;
        
        double A = 1.0 + 4.0 * normTemp + 16.0 * normTemp*normTemp;

        int k=0;
        st.photon.ps.iterate([&](const SpaceIterator& j) {
            
            double energy = j.val(DIM_E);
            frecuency[k] = energy / planck;
                        
            double jBr = st.tpf1.get(j) * 0.25 * energy*energy / pi;
            double jSy = st.tpf2.get(j) * 0.25 * energy*energy / pi;
                        
            lumBr[k] += jBr*emissToLum;
            lumSy[k] += (jSy * emissToLum);
            
            // INVERSE COMPTON /////////////////////////////////////////////////
            
            double jICB = 0.0;
            double jICS = 0.0;
            double jBrSource, jSySource;
            
            double auxEnergy = energy;
            for (int s = 1; s < 6; s++) {
                auxEnergy = auxEnergy / A;
                
				if (temp > 5.e8) {
					jSySource = jSync(auxEnergy, temp, magfield, dens_e);
				} else {
					jSySource = 0.0;
				}
				
				jBrSource = st.tpf1.interpolate({ {DIM_E, auxEnergy}  }, &(const SpaceCoord&) j)* 0.25*auxEnergy*auxEnergy / pi;
            
                if ( bb_RJ(auxEnergy/planck, temp) < jSySource * 4.0*pi*r/3.0 || auxEnergy > 3.0*boltzmann*temp) {
                    jSySource = 0.0;
                    jBrSource = 0.0;
                }
                
                jICB += jIC(jBrSource, normTemp, r/rg, rMax, theta, j, st.denf_e, s);
                jICS += jIC(jSySource, normTemp, r/rg, rMax, theta, j, st.denf_e, s);
            }

            lumICB[k] += jICB * emissToLum;
            lumICS[k] += jICS * emissToLum;
            
            ///////////////////////////////////////////////////////////////////
			
			//jIC_Syn[k] += thComptonCoppi(energy, st.electron, j, st.tpf1, st.photon.emin())*energy*fluxToLum*cLight*
			//					escapePhoton(energy, st.electron, st.denf_e, i);

            k++;
                            
        }, { -1, i.coord[DIM_R], i.coord[DIM_THETA] } );
    }, { 0, -1, -1 } );
    
    for (int k=0; k<nEnergy; k++) {
        file << "\t" << log10(frecuency[k])
                           << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency[k] * lumBr[k]
                           << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency[k] * lumSy[k]
                           << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency[k] * lumICB[k]
                           << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << frecuency[k] * lumICS[k]

						  // << setw(10) << setiosflags(ios::fixed) << scientific << setprecision(2) << jIC_Syn[k]

                           << std::endl;
    };
     */
     
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
    
	show_message(msgEnd, Module_luminosities);
	
}