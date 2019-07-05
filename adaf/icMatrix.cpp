#include "icMatrix.h"
#include <stdio.h>
#include "modelParameters.h"
#include "globalVariables.h"
#include <fparameters/SpaceIterator.h>
#include "adafFunctions.h"
#include <boost/property_tree/ptree.hpp>
#include <fparameters/Dimension.h>
#include <fparameters/parameters.h>

#include <fmath/constants.h>
#include <fmath/RungeKutta.h>
#include <stdio.h>
#include "write.h"

extern "C" {
	#include <nrMath/random.h>
}

#define RANDOM_GENERATOR gsl_rng_r250 /* R250 random number generator from GNU Scientific Library. */

#include <boost/property_tree/ptree.hpp>
#include <fparameters/Dimension.h>
#include <fparameters/parameters.h>

#define new_max(x,y) ((x) >= (y)) ? (x) : (y)
#define new_min(x,y) ((x) <= (y)) ? (x) : (y)
#define new_abs(x,y) ((x) >= (y)) ?(x-y) : (y-x)

void icMatrixWrite()
{
	matrixWrite("scattAA.dat",scattAA,nR,nR);
	matrixWrite("scattDA.dat",scattDA,nRcd,nR);
	matrixWrite("reachAD.dat",reachAD,nR,nRcd);
	matrixWrite("reachDA.dat",reachDA,nRcd,nR);
	matrixWrite("reachAA.dat",reachAA,nR,nR);
	vectorWrite("escapeAi.dat",escapeAi,nR);
	vectorWrite("escapeDi.dat",escapeDi,nRcd);
}

void icMatrix(State& st)
{
	size_t nPhot = GlobalConfig.get<size_t>("nPhotMatrix");
	size_t nTheta = GlobalConfig.get<size_t>("nTheta");

	Vector rCellsBoundaries(nR+1,0.0), rCellsBoundariesCD(nRcd+1,0.0);
	rCellsBoundaries[0] = st.electron.ps[DIM_R][0]/sqrt(paso_r);
	rCellsBoundariesCD[0]=rTr;
    for(size_t iR=1;iR<=nR;iR++) rCellsBoundaries[iR]=rCellsBoundaries[iR-1]*paso_r;
	for(size_t iRcd=1;iRcd<=nRcd;iRcd++) rCellsBoundariesCD[iRcd]=rCellsBoundariesCD[iRcd-1]*paso_rCD;

	matrixInit(scattAA,nR,nR,0.0);
	matrixInit(scattDA,nRcd,nR,0.0);
	matrixInit(reachAD,nR,nRcd,0.0);
	matrixInit(reachAA,nR,nR,0.0);
	matrixInit(reachDA,nRcd,nR,0.0);
    escapeAi.resize(nR,0.0);
	escapeDi.resize(nRcd,0.0);
    
    InitialiseRandom(RANDOM_GENERATOR);
	size_t iR=0;
	//ADAF SCATTERING MATRIX
	
	double pasoprim = pow(rCellsBoundaries[nR]/rCellsBoundaries[0],1.0/(nR*10.0));
	double pasoprimmin = pow(rCellsBoundaries[nR]/rCellsBoundaries[0],1.0/(nR*100.0));
	double pasoprimmax = pow(rCellsBoundaries[nR]/rCellsBoundaries[0],1.0/(0.5*nR));
	
    InitialiseRandom(RANDOM_GENERATOR);
	st.photon.ps.iterate([&](const SpaceIterator& it1) {
		double r0=it1.val(DIM_R);
		double thetaMin=st.thetaH.get(it1);
		double dyaux=(1.0-sin(thetaMin))/nTheta;
        for(size_t kTh=1;kTh<=nTheta;kTh++) {
            double yaux=sin(thetaMin)+kTh*dyaux;  // Theta distributed uniformly in sin(theta).
            double theta0=asin(yaux);
			double y0=r0*sin(theta0);
            double z0=r0*cos(theta0);
			for(size_t jPh=1;jPh<=nPhot;jPh++) {
                double random_number = gsl_rng_uniform(RandomNumberGenerator);
                double phiprim = 2.0*pi*random_number;
                random_number = gsl_rng_uniform(RandomNumberGenerator);
                double thetaprim = acos(1.0-2.0*random_number);   // Photon directions
																  // distributed
																	// isotropically.
				double pescap = 1.0;
                double ne;
                double r1=r0;
				double theta1;
				double z1,z1ant;
				z1ant = z0;
				double drprim=r0*(pasoprimmin-1.0);       // Minimum step.
				double rprim=drprim;
				vector<size_t> countAA(nR,0);
                do {
					double xprim=rprim*sin(thetaprim)*cos(phiprim);
					double yprim=rprim*sin(thetaprim)*sin(phiprim);
					double zprim=rprim*cos(thetaprim);
					double x1=xprim;
					double y1=y0+yprim;
					z1=z0+zprim;
					r1=sqrt(x1*x1+y1*y1+z1*z1);
					theta1=atan(sqrt(x1*x1+y1*y1)/abs(z1));
					ne = electronDensityTheta(r1,theta1);
					double exptau = exp(-ne*thomson*drprim);
                    double psc = 1.0-exptau;    // Probability of scattering.
					if (r1 > rTr && z1*z1ant < 0.0) {
						for(size_t jRcd=0;jRcd<=nRcd;jRcd++) {
							if(r1 > rCellsBoundariesCD[jRcd] && r1 < rCellsBoundariesCD[jRcd+1]) {
								reachAD[iR][jRcd] += pescap; 		// Add the probability of
																	// absorption by the ring
																	// jRcd of the thin disk.
								goto LOOP;
							}
						}
                    }
					z1ant = z1;
                    for(size_t jR=0;jR<=nR;jR++) {
                        if(r1 > rCellsBoundaries[jR] && r1 < rCellsBoundaries[jR+1]) {
							scattAA[iR][jR] += psc*pescap;
							if (countAA[jR] == 0) {
								reachAA[iR][jR] += pescap;
								countAA[jR]++;
							}
						}
                    }
					pescap *= exptau;
                    drprim=r1*(pasoprimmin-1.0);
                    rprim += drprim;
                } while( r1 < rCellsBoundaries[nR] );     // Escape from the region.
				LOOP:;
            }
        }
        for(size_t jR=0;jR<nR;jR++) {
            scattAA[iR][jR] /= (nPhot*nTheta);     // Dividing by the number of photons launched.
			reachAA[iR][jR] /= (nPhot*nTheta);
        }
		for (size_t jRcd=0;jRcd<nRcd;jRcd++) {
			reachAD[iR][jRcd] /= (nPhot*nTheta);
		}
		iR++;
    },{0,-1,0});
	
	//////////////////////////////////////////////////////////////////////////////////////
	
	iR=0;
	// ADAF ESCAPE PHOTONS
    st.photon.ps.iterate([&](const SpaceIterator& it2) {
		double r0=it2.val(DIM_R);
		double thetaMin=st.thetaH.get(it2);
		double dyaux=(1.0-sin(thetaMin))/nTheta;
        for(size_t kTh=1;kTh<=nTheta;kTh++) {
            double yaux=sin(thetaMin)+kTh*dyaux;  // Theta distributed uniformly in sin(theta).
            double theta0 = (yaux < 1.0 ? asin(yaux) : pi/2.0);
			double phi0 = 0.0;
			double dPhi = 2.0*pi/nPhot;
			double x0=r0*sin(theta0)*cos(phi0);
			double y0=r0*sin(theta0)*sin(phi0);
            double z0=r0*cos(theta0);
			for(size_t jPh=1;jPh<=nPhot;jPh++) {
                double thetaprim = inclination*(pi/180.0);
				double u = sin(thetaprim);
				double w = cos(thetaprim);
				size_t Nk = 100;
				double rMax = rCellsBoundaries[nR];
				double paso_r = pow(rMax/r0,1.0/Nk);
				double rNow = r0;
				double xNow = x0; double yNow = y0; double zNow = z0;
				double thetaNow = theta0;
				size_t k = 0;
				double tau = 0.0;
				while (k<Nk && (electronDensityTheta(rNow,thetaNow)>1.0e-3)){
					double drNow = rNow*(paso_r-1.0);
					xNow += u*drNow;
					zNow += w*drNow;
					rNow = sqrt(xNow*xNow+yNow*yNow+zNow*zNow);
					thetaNow = asin(abs(zNow)/rNow);
					tau += thomson*electronDensityTheta(rNow,thetaNow)*drNow;
					k++;
				}
				escapeAi[iR] += exp(-tau);
				phi0 += dPhi;
            }
        }
		escapeAi[iR] /= (nPhot*nTheta);     // Dividing by the number of photons launched.
		iR++;
    },{0,-1,0});
	
	//////////////////////////////////////////////////////////////////////////////////////

	// COLD DISK SCATTERING MATRIX
	pasoprim = pow(rCellsBoundariesCD[nRcd]/rCellsBoundariesCD[0],1.0/(nRcd*10.0));
	pasoprimmin = pow(rCellsBoundariesCD[nRcd]/rCellsBoundariesCD[0],1.0/(nRcd*100.0));
	pasoprimmax = pow(rCellsBoundariesCD[nRcd]/rCellsBoundariesCD[0],1.0/(0.5*nRcd));
	
	size_t iRcd=0;
	st.photon.ps.iterate([&](const SpaceIterator& it2) {
		double r0cd=it2.val(DIM_Rcd);
		double drprim=r0cd*(pasoprim-1.0);             // Initial step for the photon path.
		double drprimmin=r0cd*(pasoprimmin-1.0);       // Minimum step.
		for(size_t jPh=1;jPh<=nPhot;jPh++) {
			double random_number = gsl_rng_uniform(RandomNumberGenerator);
			double phiprim = 2.0*pi*random_number;
			random_number = gsl_rng_uniform(RandomNumberGenerator);
			double thetaprim = acos(1.0-2.0*random_number);   // Photon directions
															  // distributed
			double rprim=drprim; 	                          // isotropically.
			double rprimant=0.0;
			double pescap = 1.0;
			double ne;
			double r1=r0cd;
			double theta1;
			double neant = electronDensity(r0cd);
			vector<size_t> countDA(nR,0);
			do {
				size_t count=0;
				double control;
				do {
					if(count >= 1) {              // Control of the step.
						pasoprim = sqrt(pasoprim);
						drprim= new_max(r1*(pasoprim-1.0),drprimmin);
						rprim = rprimant+drprim;
					}
					double xprim=rprim*sin(thetaprim)*cos(phiprim);
					double yprim=rprim*sin(thetaprim)*sin(phiprim);
					double zprim=rprim*cos(thetaprim);
					double x1=xprim;
					double y1=r0cd+yprim;
					double z1=zprim;
					r1=sqrt(x1*x1+y1*y1+z1*z1);
					theta1=atan(sqrt(x1*x1+y1*y1)/abs(z1));
					double thetaMinLocal = acos(costhetaH(r1));
					ne = (theta1 > thetaMinLocal && theta1 < (pi-thetaMinLocal)) 
												? electronDensity(r1) : 0.0;
					control = abs(ne-neant)/(0.5*(ne+neant));       // dn/n
					count++;
				} while( (control > 1.0e-2) && (drprim-drprimmin > 1.0e-9) );
				double exptau = exp(-ne*thomson*drprim);
				double psc = 1.0-exptau;    // Probability of scattering.
				for(size_t jR=0;jR<=nR;jR++) {
					if(r1 > rCellsBoundaries[jR] && r1 < rCellsBoundaries[jR+1]) {
						scattDA[iRcd][jR] += psc*pescap;
						if (countDA[jR] == 0) {
							reachDA[iRcd][jR] += pescap;
							countDA[jR]++;
						}
					}
				}
				neant=ne;
				rprimant=rprim;
				pasoprim=new_min(pasoprim*pasoprim,pasoprimmax);
				drprim=r1*(pasoprim-1.0);
				rprim += drprim;
				pescap *= exptau;                // Probability that a photon reaches 
												 // the previous position.
			} while( r1 < rCellsBoundaries[nR] );     // Escape from the region.
        }
        for(size_t jR=0;jR<nR;jR++) {
            scattDA[iRcd][jR] /= nPhot;     // Dividing by the number of photons launched.
			reachDA[iRcd][jR] /= nPhot;
        }
		iRcd++;
    },{0,0,-1});
	
	//COLD DISK ESCAPE PHOTONS
	pasoprim = pow(rCellsBoundariesCD[nRcd]/rCellsBoundariesCD[0],1.0/(nRcd*10.0));
	iRcd=0;
	st.photon.ps.iterate([&](const SpaceIterator& it2) {
		double r0cd=it2.val(DIM_Rcd);
		double drprim=r0cd*(pasoprim-1.0);             // Initial step for the photon path.
		double drprimmin=r0cd*(pasoprimmin-1.0);       // Minimum step.
		for(size_t jPh=1;jPh<=nPhot;jPh++) {
			double random_number = gsl_rng_uniform(RandomNumberGenerator);
			double phiprim = 2.0*pi*random_number;
			double thetaprim = inclination*(pi/180.0);
			
			double rprim=drprim; 	                          
			double rprimant=0.0;
			double pescap = 1.0;
			double ne;
			double r1=r0cd;
			double theta1;
			double neant = electronDensity(r0cd); 
			do {
				size_t count=0;
				double control;
				do {
					if(count >= 1) {              // Control of the step.
						pasoprim = sqrt(pasoprim);
						drprim= new_max(r1*(pasoprim-1.0),drprimmin);
						rprim = rprimant+drprim;
					}
					double xprim=rprim*sin(thetaprim)*cos(phiprim);
					double yprim=rprim*sin(thetaprim)*sin(phiprim);
					double zprim=rprim*cos(thetaprim);
					double x1=xprim;
					double y1=r0cd+yprim;
					double z1=zprim;
					r1=sqrt(x1*x1+y1*y1+z1*z1);
					theta1=atan(sqrt(x1*x1+y1*y1)/abs(z1));
					double thetaMinLocal = acos(costhetaH(r1));
					ne = (theta1 > thetaMinLocal && theta1 < (pi-thetaMinLocal)) 
												? electronDensity(r1) : 0.0;
					control = abs(ne-neant)/(0.5*(ne+neant));       // dn/n
					count++;
				} while( (control > 1.0e-2) && (drprim-drprimmin > 1.0e-9) );
				pescap *= exp(-ne*thomson*drprim);
				neant=ne;
				rprimant=rprim;
				pasoprim=new_min(pasoprim*pasoprim,pasoprimmax);
				drprim=r1*(pasoprim-1.0);
				rprim += drprim;
			} while( r1 < rCellsBoundaries[nR] );     // Escape from the region.
			escapeDi[iRcd] += pescap;
        }
        escapeDi[iRcd] /= nPhot;
		iRcd++;
    },{0,0,-1});
    FinaliseRandom();
	icMatrixWrite();
}

void icMatrixRead(State& st)
{
	matrixInit(scattAA,nR,nR,0.0);
	matrixInit(scattDA,nRcd,nR,0.0);
	matrixInit(reachAD,nR,nRcd,0.0);
	matrixInit(reachAA,nR,nR,0.0);
	matrixInit(reachDA,nRcd,nR,0.0);
    escapeAi.resize(nR,0.0);
	escapeDi.resize(nRcd,0.0);
	
	matrixRead("scattAA.dat",scattAA,nR,nR);
	matrixRead("scattDA.dat",scattDA,nRcd,nR);
	matrixRead("reachAD.dat",reachAD,nR,nRcd);
	matrixRead("reachAA.dat",reachAA,nR,nRcd);
	matrixRead("reachDA.dat",reachDA,nR,nRcd);
	vectorRead("escapeAi.dat",escapeAi,nR);
	vectorRead("escapeDi.dat",escapeDi,nRcd);
}