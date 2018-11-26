#include "modelParameters.h"
#include "torusSampling.h"
#include "globalVariables.h"
#include <fparameters/SpaceIterator.h>
#include "torusFunctions.h"
#include <boost/property_tree/ptree.hpp>
#include <fparameters/Dimension.h>
#include <fparameters/parameters.h>

extern "C" {
	#include <nrMath/nrutil.h>
	#include <nrMath/random.h>
}

#include <fmath/constants.h>

#include <stdio.h>

#define RANDOM_GENERATOR gsl_rng_r250 /* R250 random number generator from GNU Scientific Library. */

#define new_max(x,y) ((x) >= (y)) ? (x) : (y)
#define new_min(x,y) ((x) <= (y)) ? (x) : (y)

void torusSampling(State& st, Matrix& scatt, Vector& escape)
{
    size_t nPhot=GlobalConfig.get<size_t>("nPhotMatrix");
    
    double dr=(edgeRadius-cuspRadius)/nR;      	        // Cells' radial size.
	Vector rCellsBoundaries(nR+1,0.0);					// Cells' boundaries.
	rCellsBoundaries[0]=cuspRadius;
    for(size_t i=1;i<=nR;i++)
        rCellsBoundaries[i]=rCellsBoundaries[i-1]+dr;

	matrixInit(scatt,nR,nR,0.0);
    escape.resize(nR,0.0);
    
    InitialiseRandom(RANDOM_GENERATOR);
	size_t iR=0;
	st.photon.ps.iterate([&](const SpaceIterator& it1) {
        double dyaux=(sin(maxPolarAngle)-sin(minPolarAngle))/nTheta;
		double r=it1.val(DIM_R);
        for(size_t k=1;k<=nTheta+1;k++) {
            double yaux=sin(minPolarAngle)+(k-1)*dyaux;   // Theta distributed uniformly in sin(theta).
            double theta0=asin(yaux);
            double y0=r*cos(theta0);
            double z0=r*sin(theta0);
            double drprim=dr/10.0;                        // Initial step for the photon path.
            double drprimmin=dr/1.0e6;                    // Minimum step.
            double drprimmax=dr/2.0;                      // Maximum step.
            for(size_t jPh=1;jPh<=nPhot;jPh++) {
                double random_number;
                random_number=gsl_rng_uniform(RandomNumberGenerator);
                double phiprim=2.0*pi*random_number;
                random_number=gsl_rng_uniform(RandomNumberGenerator);
                double thetaprim=acos((pi/2.0-random_number*pi)/(pi/2.0)); // Photon directions distributed
                double rprim=drprim;                                                                                       // isotropically.
                double rprimant=0.0;
                double accumulatedp=0.0;
                double neant, ne;
                double r1,theta1;
				neant=electronDensity(r,theta0);
				double pescap=1.0;
                do {
                    size_t count=0;
                    double control;
                    do {
                        if(count >= 1) {                                                             // Control of the step.
                            drprim = new_max(drprim/2.0,drprimmin);
                            rprim = rprimant+drprim;
                        }
                        double xprim=rprim*sin(thetaprim)*cos(phiprim);
                        double yprim=rprim*sin(thetaprim)*sin(phiprim);
                        double zprim=rprim*cos(thetaprim);
                        double x1=xprim;
                        double y1=y0+yprim;
                        double z1=z0+zprim;
                        r1=sqrt(x1*x1+y1*y1+z1*z1);
                        theta1=atan(abs(z1)/sqrt(x1*x1+y1*y1));
                        ne=electronDensity(r1,theta1);
                        control = abs(ne-neant)/(0.5*(ne+neant));                 // dn/n
                        count++;
                    } while((control > 1.0e-2) && (drprim-drprimmin > 1.0e-9));
                    double psc=ne*thomson*(drprim*gravRadius);    // Probability of scattering.
                    pescap=1.0-accumulatedp;       // Probability that a photon reaches the previous position.
                    for(size_t jR=0;jR<nR;jR++) {
                        if (r1 > rCellsBoundaries[jR] && r1 < rCellsBoundaries[jR+1]) 
							scatt[iR][jR] += psc*pescap;            // Add the probability
																    // interaction in the matrix
																	// element.
                    }
                    accumulatedp += psc;    // The product of (1-p_k) gives the probability for a photon to reach the
											// previous position. To first order, this is 1-sum(p_k). Hence, we accumulate
											// the sum of the p_k.
                    neant=ne;
                    rprimant=rprim;
                    drprim=new_min(drprim*2,drprimmax);
                    rprim+=drprim;
                } while (ne > 1.0e1);                          // Escape from the torus.
				pescap = 1.0-accumulatedp;
				escape[iR] += pescap;
            }
        }
		escape[iR] /= (nPhot*(nTheta+1));
        for(size_t jR=0;jR<nR;jR++) {
            scatt[iR][jR] /= (nPhot*(nTheta+1));      // Dividing by the number of photons launched.
        }
		iR++;
    },{0,-1,0});
    FinaliseRandom();
}