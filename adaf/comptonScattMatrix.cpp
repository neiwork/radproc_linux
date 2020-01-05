#include "comptonScattMatrix.h"
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

void comptonScattMatrixWrite()
{
	matrixWrite("scattAA.dat",scattAA,nR,nR);
	matrixWrite("scattDA.dat",scattDA,nRcd,nR);
	matrixWrite("reachAD.dat",reachAD,nR,nRcd);
	matrixWrite("reachDA.dat",reachDA,nRcd,nR);
	matrixWrite("reachAA.dat",reachAA,nR,nR);
	vectorWrite("escapeAi.dat",escapeAi,nR);
	vectorWrite("escapeDi.dat",escapeDi,nRcd);
}

void comptonScattMatrix(State& st)
{
	size_t nPhot = GlobalConfig.get<size_t>("scatt.nRandomPhot");
	size_t nTheta = GlobalConfig.get<size_t>("scatt.nTheta");
    
	double captureRadius = 0.5*sqrt(27.0)*schwRadius;

	Vector rCellsBoundaries(nR+1,0.0), rCellsBoundariesCD(nRcd+1,0.0);
	rCellsBoundaries[0] = st.denf_e.ps[DIM_R][0]/sqrt(paso_r);
	rCellsBoundariesCD[0] = rTr;
    for(size_t iR=1;iR<=nR;iR++) rCellsBoundaries[iR]=rCellsBoundaries[iR-1]*paso_r;
	for(size_t iRcd=1;iRcd<=nRcd;iRcd++) rCellsBoundariesCD[iRcd]=rCellsBoundariesCD[iRcd-1]*paso_rCD;

	double rBound = max(rCellsBoundaries[nR],rCellsBoundariesCD[nRcd]);
	matrixInit(scattAA,nR,nR,0.0);
	matrixInit(scattDA,nRcd,nR,0.0);
	matrixInit(reachAD,nR,nRcd,0.0);
	matrixInit(reachAA,nR,nR,0.0);
	matrixInit(reachDA,nRcd,nR,0.0);
    escapeAi.resize(nR,0.0);
	escapeDi.resize(nRcd,0.0);
    
    InitialiseRandom(RANDOM_GENERATOR);
	//ADAF SCATTERING MATRIX
	
	double pasoprim = pow(rCellsBoundaries[1]/rCellsBoundaries[0],1.0/5.0);
	
	#pragma omp parallel for
	for (int iR=0;iR<nR;iR++) {
		double r0 = st.denf_e.ps[DIM_R][iR];
		double thetaMin;
		if (height_method == 0)
			thetaMin = st.thetaH.get({0,iR,0});
		else
			thetaMin = atan(r0/height_fun(r0));
		double dyaux=(1.0-sin(thetaMin))/nTheta;
        for(size_t kTh=1;kTh<=nTheta;kTh++) {
            double yaux=sin(thetaMin)+kTh*dyaux;  // Theta distributed uniformly in sin(theta).
            double theta0=asin(yaux);
			double y0,z0;
			if (height_method == 0) {
				y0=r0*sin(theta0);
				z0=r0*cos(theta0);
			} else {
				y0 = r0;
				z0 = r0/tan(theta0);
			}
			for(size_t jPh=1;jPh<=nPhot;jPh++) {
				double random_number = gsl_rng_uniform(RandomNumberGenerator);
				double phiprim = 2.0*pi*random_number;
				random_number = gsl_rng_uniform(RandomNumberGenerator);
				double thetaprim = acos(1.0-2.0*random_number);   // Photon directions
																  // distributed
																  // isotropically.
				double pescap = 1.0;
				double z1ant = z0;
				double drprim=r0*(pasoprim-1.0);       // Step.
				double rprim=drprim;
				vector<size_t> countAA(nR,0);
				double r1 = r0;
				double z1 = z0;
				double r1aux = r1;
                do {
					double xprim=rprim*sin(thetaprim)*cos(phiprim);
					double yprim=rprim*sin(thetaprim)*sin(phiprim);
					double zprim=rprim*cos(thetaprim);
					double x1=xprim;
					double y1=y0+yprim;
					z1=z0+zprim;
					r1=sqrt(x1*x1+y1*y1+z1*z1);
					r1aux = (height_method == 0) ? r1 : sqrt(x1*x1+y1*y1);
					double theta1=atan(sqrt(x1*x1+y1*y1)/abs(z1));
					double ne = electronDensityTheta(r1,theta1);
					double exptau = exp(-ne*thomson*drprim);
                    double psc = 1.0-exptau;    // Probability of scattering.
					if (r1aux > rTr && z1*z1ant < 0.0) {
						for(size_t jRcd=0;jRcd<=nRcd;jRcd++) {
							if(r1aux > rCellsBoundariesCD[jRcd] && r1aux < rCellsBoundariesCD[jRcd+1]) {
								reachAD[iR][jRcd] += pescap; 		// Add the probability of
																	// absorption by the ring
																	// jRcd of the thin disk.
								goto LOOP;
							}
						}
                    }
					z1ant = z1;
					for(size_t jR=0;jR<nR;jR++) {
						
						if(r1aux > rCellsBoundaries[jR] && r1aux < rCellsBoundaries[jR+1]) {
							
							scattAA[iR][jR] += psc*pescap;
							if (countAA[jR] == 0) {
								reachAA[iR][jR] += pescap;
								countAA[jR]++;
							}
						}
                    }
					pescap *= exptau;
                    drprim=r1*(pasoprim-1.0);
                    rprim += drprim;
                } while(r1aux < rBound && z1 < rBound && r1 > schwRadius);
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
	}
	
	//////////////////////////////////////////////////////////////////////////////////////
	
	// ADAF ESCAPE PHOTONS
	#pragma omp parallel for
	for (int iR=0;iR<nR;iR++) {
		double r0 = st.thetaH.ps[DIM_R][iR];
		double drprim=r0*(pasoprim-1.0);
		double thetaMin;
		if (height_method == 0)
			thetaMin = st.thetaH.get({0,iR,0});
		else
			thetaMin = atan(r0/height_fun(r0));
		double dyaux=(1.0-sin(thetaMin))/nTheta;
        for(size_t kTh=1;kTh<=nTheta;kTh++) {
            double yaux=sin(thetaMin)+kTh*dyaux;  // Theta distributed uniformly in sin(theta).
            double theta0 = (yaux < 1.0 ? asin(yaux) : pi/2.0);
			double phi0 = 0.0;
			double dPhi = 2.0*pi/nPhot;
			
			double x0,y0,z0;
			if (height_method == 0) {
				x0 = r0*sin(theta0)*cos(phi0);
				y0 = r0*sin(theta0)*cos(phi0);
				z0 = r0*cos(theta0);
			} else {
				x0 = r0*sin(phi0);
				y0 = r0*cos(phi0);
				z0 = r0/tan(theta0);
			}
			for(size_t jPh=1;jPh<=nPhot;jPh++) {
				double thetaprim = inclination*(pi/180.0);
				double rprim=drprim;
				double r1 = r0;
				double z1 = z0;
				double r1aux = r1;
				double pescap = 1.0;
				do {
					double xprim=rprim*sin(thetaprim);
					double zprim=rprim*cos(thetaprim);
					double x1=x0+xprim;
					z1=z0+zprim;
					r1=sqrt(x1*x1+y0*y0+z1*z1);
					double theta1=atan(sqrt(x1*x1+y0*y0)/abs(z1));
					r1aux = (height_method == 0) ? r1 : sqrt(x1*x1+y0*y0);
					double ne = electronDensityTheta(r1,theta1);
					pescap *= exp(-ne*thomson*drprim);
					drprim=r1*(pasoprim-1.0);
					rprim += drprim;
				} while (r1aux < rCellsBoundaries[nR] && z1 < rCellsBoundaries[nR]);
				escapeAi[iR] += pescap;
				phi0 += dPhi;
            }
        }
		escapeAi[iR] /= (nPhot*nTheta);     // Dividing by the number of photons launched.
	}

	//////////////////////////////////////////////////////////////////////////////////////

	// COLD DISK SCATTERING MATRIX
	#pragma omp parallel for
	for (int iRcd=0;iRcd<nRcd;iRcd++) {
		double r0cd = st.denf_e.ps[DIM_Rcd][iRcd];
		double drprim = schwRadius*(pasoprim-1.0);             // Step for the photon path.
		for(size_t jPh=1;jPh<=nPhot;jPh++) {
			double random_number = gsl_rng_uniform(RandomNumberGenerator);
			double phiprim = 2.0*pi*random_number;
			random_number = gsl_rng_uniform(RandomNumberGenerator);
			double thetaprim = acos(1.0-2.0*random_number);   // Photon directions
															  // distributed
			double rprim=drprim; 	                          // isotropically.
			double pescap = 1.0;
			double r1 = r0cd;
			double r1aux = r1;
			double z1 = 0.0;
			vector<size_t> countDA(nR,0);
			do {
				double xprim=rprim*sin(thetaprim)*cos(phiprim);
				double yprim=rprim*sin(thetaprim)*sin(phiprim);
				double zprim=rprim*cos(thetaprim);
				double x1=xprim;
				double y1=r0cd+yprim;
				z1=zprim;
				r1=sqrt(x1*x1+y1*y1+z1*z1);
				double theta1=atan(sqrt(x1*x1+y1*y1)/abs(z1));
				r1aux = (height_method == 0) ? r1 : sqrt(x1*x1+y1*y1);
				double exptau = 1.0;
				if (r1aux < rCellsBoundaries[nR]) {
					double ne = electronDensityTheta(r1,theta1);
					exptau = exp(-ne*thomson*drprim);
					double psc = 1.0-exptau;    // Probability of scattering.
					for(size_t jR=0;jR<=nR;jR++) {
						if(r1aux > rCellsBoundaries[jR] && r1aux < rCellsBoundaries[jR+1]) {
							scattDA[iRcd][jR] += psc*pescap;
							if (countDA[jR] == 0) {
								reachDA[iRcd][jR] += pescap;
								countDA[jR]++;
							}
						}
					}
				}
				//drprim=r1*(pasoprim-1.0);
				rprim += drprim;
				pescap *= exptau;                // Probability that a photon reaches 
												 // the previous position.
			} while(r1aux < rBound && z1 < rBound && r1 > schwRadius);
        }
        for(size_t jR=0;jR<nR;jR++) {
            scattDA[iRcd][jR] /= nPhot;     // Dividing by the number of photons launched.
			reachDA[iRcd][jR] /= nPhot;
        }
	}

	//COLD DISK ESCAPE PHOTONS
	#pragma omp parallel for
	for (int iRcd=0;iRcd<nRcd;iRcd++) {
		double r0cd = st.denf_e.ps[DIM_Rcd][iRcd];
		double drprim=r0cd*(pasoprim-1.0);
		for(size_t jPh=1;jPh<=nPhot;jPh++) {
			double random_number = gsl_rng_uniform(RandomNumberGenerator);
			double phiprim = 2.0*pi*random_number;
			double thetaprim = inclination*(pi/180.0);
			double rprim=drprim;
			double r1 = r0cd;
			double z1 = 0.0;
			double r1aux = r1;
			double pescap = 1.0;
			do {
				double xprim=rprim*sin(thetaprim)*cos(phiprim);
				double yprim=rprim*sin(thetaprim)*sin(phiprim);
				double zprim=rprim*cos(thetaprim);
				double x1=xprim;
				double y1=r0cd+yprim;
				z1=zprim;
				r1=sqrt(x1*x1+y1*y1+z1*z1);
				double theta1=atan(sqrt(x1*x1+y1*y1)/abs(z1));
				r1aux = (height_method == 0) ? r1 : sqrt(x1*x1+y1*y1);
				double ne = electronDensityTheta(r1,theta1);
				pescap *= exp(-ne*thomson*drprim);
				drprim=r1*(pasoprim-1.0);
				rprim += drprim;
			} while(r1aux < rBound && z1 < rBound);     // Escape from the region.
			escapeDi[iRcd] += pescap;
        }
        escapeDi[iRcd] /= nPhot;
	}

    FinaliseRandom();
	comptonScattMatrixWrite();
}

void comptonScattMatrixRead(State& st)
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
	matrixRead("reachAA.dat",reachAA,nR,nR);
	matrixRead("reachDA.dat",reachDA,nRcd,nR);
	vectorRead("escapeAi.dat",escapeAi,nR);
	vectorRead("escapeDi.dat",escapeDi,nRcd);
}