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
#include <stdio.h>

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

void icMatrix(State& st, Matrix& scatt, Vector& escape)
{
	size_t nPhot = GlobalConfig.get<size_t>("nPhotMatrix");
	size_t nTheta = GlobalConfig.get<size_t>("nTheta");

	Vector rCellsBoundaries(nR+1,0.0);					// Cells' boundaries.
	rCellsBoundaries[0]=schwRadius;
    for(size_t iR=1;iR<=nR;iR++)
        {rCellsBoundaries[iR]=rCellsBoundaries[iR-1]*paso_r;}

	matrixInit(scatt,nR,nR,0.0);
    escape.resize(nR,0.0);
    
    InitialiseRandom(RANDOM_GENERATOR);
	size_t iR=0;
	
	double pasoprim=pow(rCellsBoundaries[nR]/rCellsBoundaries[0],1.0/(nR*10.0));
	double pasoprimmin=pow(rCellsBoundaries[nR]/rCellsBoundaries[0],1.0/(nR*100.0));
	double pasoprimmax=pow(rCellsBoundaries[nR]/rCellsBoundaries[0],1.0/(0.5*nR));
			
    InitialiseRandom(RANDOM_GENERATOR);
    st.photon.ps.iterate([&](const SpaceIterator& it1) {
		double r=it1.val(DIM_R);
		double thetaMin=st.thetaH.get(it1);
		double dyaux=(1.0-sin(thetaMin))/nTheta;
        for(size_t k=1;k<=nTheta;k++) {
            double yaux=sin(thetaMin)+k*dyaux;    // Theta distributed uniformly in sin(theta).
            double theta0=asin(yaux);
			double y0=r*sin(theta0);
            double z0=r*cos(theta0);
            double drprim=r*(pasoprim-1.0);             // Initial step for the photon path.
			double drprimmin=r*(pasoprimmin-1.0);       // Minimum step.
			for(size_t jPh=1;jPh<=nPhot;jPh++) {
                double random_number;
                random_number=gsl_rng_uniform(RandomNumberGenerator);
                double phiprim=2.0*pi*random_number;
                random_number=gsl_rng_uniform(RandomNumberGenerator);
                double thetaprim=acos(1.0-2.0*random_number);        // Photon directions distributed
                double rprim=drprim;                                 // isotropically.
                double rprimant=0.0;
				double pescap = 1.0;
                double accumulatedp=0.0;
                double ne;
                double r1,theta1;
                
				double neant = electronDensity(r); 
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
                        double y1=y0+yprim;
                        double z1=z0+zprim;
                        r1=sqrt(x1*x1+y1*y1+z1*z1);
                        theta1=atan(sqrt(x1*x1+y1*y1)/abs(z1));
						if (theta1 > thetaMin) {
							ne = electronDensity(r1);
						} else ne = 0.0;
						
                        control = abs(ne-neant)/(0.5*(ne+neant));       // dn/n
                        count++;
                    } while((control > 1.0e-2) && (drprim-drprimmin > 1.0e-9));
                    double psc=ne*thomson*drprim;    // Probability of scattering.
                    pescap=1.0-accumulatedp;       // Probability that a photon reaches the previous position.
                    for(size_t jR=0;jR<=nR;jR++) {
                        if(r1 > rCellsBoundaries[jR] && r1 < rCellsBoundaries[jR+1]) 
							scatt[iR][jR] += psc*pescap; 		// Add the probability of
															// interaction in the matrix
															// element.
                    }
                    accumulatedp += psc;  // The product of (1-p_k) gives the probability for a photon to reach the
										  // previous position. To first order, this is 1-sum(p_k). Hence, we accumulate
										  // the sum of the p_k.
                    neant=ne;
                    rprimant=rprim;
					pasoprim=new_min(pasoprim*pasoprim,pasoprimmax);
                    drprim=r1*(pasoprim-1.0);
                    rprim += drprim;
                } while(ne > 1.0e1);               // Escape from the torus.
				pescap = 1.0-accumulatedp;
				escape[iR] += pescap;
            }
        }
		escape[iR] /= (nPhot*nTheta);
        for(size_t jR=0;jR<nR;jR++) {
            scatt[iR][jR] /= (nPhot*nTheta);         // Dividing by the number of photons launched.
        }
		iR++;
    },{0,-1,0});
    FinaliseRandom();
}