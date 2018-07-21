#include "torusSampling.h"
#include "nrutil.h"
#include "random.h"
#include <stdio.h>
#include <math.h>

#define RANDOM_GENERATOR gsl_rng_r250 /* R250 random number generator from GNU Scientific Library. */

#define new_max(x,y) ((x) >= (y)) ? (x) : (y)
#define new_min(x,y) ((x) <= (y)) ? (x) : (y)
#define new_abs(x,y) ((x) >= (y)) ?(x-y) : (y-x)

void torusSampling(double **prob)
{
	extern double rMax,rMin;
	extern int nR,nTheta;
	int nPhot=100;
	
	double pasoprim=pow(rMax/rMin,1.0/(nR*10.0));
	double pasoprimmin=pow(rMax/rMin,1.0/(nR*1000.0));
	double pasoprimmax=pow(rMax/rMin,1.0/(0.5*nR));
			
    InitialiseRandom(RANDOM_GENERATOR);
    for(int i=1;i<=nR;i++) {
        double dyaux=(sin(thetaMax)-sin(thetaMin))/nTheta;
        for(int k=1;k<=nTheta+1;k++) {
            double yaux=sin(thetaMin)+(k-1)*dyaux;                  // Theta distributed uniformly in sin(theta).
            //double theta0=asin(yaux);
            double theta0=0.0;
			double y0=r[i]*cos(theta0);
            double z0=r[i]*sin(theta0);

            double drprim=r[i]*(pasoprim-1.0);                    // Initial step for the photon path.
			double drprimmin=r[i]*(pasoprimmin-1.0);                 // Minimum step.
            double drprimmax=r[i]*(pasoprimmax-1.0);
			for(int j=1;j<=nPhot;j++) {
                double random_number;
                random_number=gsl_rng_uniform(RandomNumberGenerator);
                double phiprim=2.0*PI*random_number;
                random_number=gsl_rng_uniform(RandomNumberGenerator);
                double thetaprim=acos((PI/2.0-random_number*PI)/(PI/2.0));        // Photon directions distributed
                double rprim=drprim;                                                                                       // isotropically.
                double rprimant=0.0;
                double accumulatedp=0.0;
                double ne;
                double r1,theta1;
                double neant=eDensity(r[i],theta0);
                do {
                    int count=0;
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
                        theta1=atan(abs(z1)/sqrt(x1*x1+y1*y1));
                        ne=eDensity(r1,theta1);  //definir las funciones
                        control = new_abs(ne,neant)/(0.5*(ne+neant));                 // dn/n
                        count++;
                    } while((control > 1.0e-2) && (drprim-drprimmin > 1.0e-10));
                    double psc=ne*thomson*(drprim*RG);    // Probability of scattering.
                    double pescap=1.0-accumulatedp;       // Probability that a photon reaches the previous position.
                    for(int j=1;j<=nR;j++) {
                        if(r1 > rCells[j-1] && r1 < rCells[j]) prob[i][j] += psc*pescap; // Add the probability of
																				         // interaction in the matrix
																						 // element.
                    }
                    printf("r1 = %f, drprim=%6.3e, j = %d, i = %d\n",r1,drprim,j,i);
                    accumulatedp+=psc;    // The product of (1-p_k) gives the probability for a photon to reach the
										  // previous position. To first order, this is 1-sum(p_k). Hence, we accumulate
										  // the sum of the p_k.
                    neant=ne;
                    rprimant=rprim;
					pasoprim=new_min(pasoprim*pasoprim,pasoprimmax);
                    drprim=r1*(pasoprim-1.0);
                    rprim += drprim;
                } while(ne > 1.0);               // Escape from the torus.
            }
        }
        for(int j=1;j<=nR;j++) {
            prob[i][j] /= (nPhot*(nTheta+1));         // Dividing by the number of photons launched.
        }
    }
    FinaliseRandom();
}