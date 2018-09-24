#include "torusSampling2.h"



#include <stdio.h>
//#include <math.h>
//#include "torusSampling.h"
//#include "auxFunctions.h"
//#include "torusParameters.h"

#include <ssADAF/packageData.h>
#include <ssADAF/ssfunctions.h>

#include <boost/property_tree/ptree.hpp>
#include <fparameters/Dimension.h>
#include <fparameters/parameters.h>

extern "C" {
#include <nrMath/nrutil.h>
//#include <nrMath/nr.h>
#include <nrMath/random.h>
}

#include <fmath/constants.h>

#include <stdio.h>


//#define pi 3.14159265359
//#define thomson 6.652458734e-25
//#define RG 1.5e11
#define RANDOM_GENERATOR gsl_rng_r250 /* R250 random number generator from GNU Scientific Library. */

#define new_max(x,y) ((x) >= (y)) ? (x) : (y)
#define new_min(x,y) ((x) <= (y)) ? (x) : (y)
#define new_abs(x,y) ((x) >= (y)) ?(x-y) : (y-x)

void torusSampling2()
{
	FILE *fp2;
    fp2=fopen("prob.txt","w");
    
	//extern double *rCells, *r,paso,rMax,rMin;
	//extern int nR,nPhot,nTheta;
	
	static const double massBH=GlobalConfig.get<double>("massBH")*solarMass;
	static const double f=GlobalConfig.get<double>("f");
	static const double alphapar=GlobalConfig.get<double>("alpha");
	static const double gammapar=GlobalConfig.get<double>("gammapar");
	static const double mdot=GlobalConfig.get<double>("mdot");
	
	
	double RG= GlobalConfig.get<double>("rg");
	
	double rMin= GlobalConfig.get<double>("rMin"); //ver estos nombres
	double rMax=GlobalConfig.get<double>("rMax"); //1.0e4;
	int nR = GlobalConfig.get<double>("model.particle.default.dim.radius.samples");//100


    long nPhot=1000;
	
	
	double thetaMin=0.0;
    double thetaMax=pi/4.0;
	int nTheta= GlobalConfig.get<double>("model.particle.default.dim.theta.samples");//sale de aca o a mano?;
	
	double **prob;
	double *rCells,*r;
	
	prob=dmatrix(1,nR,1,nR);
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
                double phiprim=2.0*pi*random_number;
                random_number=gsl_rng_uniform(RandomNumberGenerator);
                double thetaprim=acos((pi/2.0-random_number*pi)/(pi/2.0));        // Photon directions distributed
                double rprim=drprim;                                                                                       // isotropically.
                double rprimant=0.0;
                double accumulatedp=0.0;
                double eDen;//este era ne pero lo reemplaze porque la funcion se llama ne
                double r1,theta1;
                
				double neant = ne(mdot, r[i], massBH, alphapar, gammapar, f);  //double neant=eDensity(r[i],theta0); 
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
                        eDen=ne(mdot, r1, massBH, alphapar, gammapar, f); //eDensity(r1,theta1);  
                        control = new_abs(eDen,neant)/(0.5*(eDen+neant));                 // dn/n
                        count++;
                    } while((control > 1.0e-2) && (drprim-drprimmin > 1.0e-10));
                    double psc=eDen*thomson*(drprim*RG);    // Probability of scattering.
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
                    neant=eDen;
                    rprimant=rprim;
					pasoprim=new_min(pasoprim*pasoprim,pasoprimmax);
                    drprim=r1*(pasoprim-1.0);
                    rprim += drprim;
                } while(eDen > 1.0);               // Escape from the torus.
            }
        }
        for(int j=1;j<=nR;j++) {
            prob[i][j] /= (nPhot*(nTheta+1));         // Dividing by the number of photons launched.
            fprintf(fp2,"%8.5e  ",prob[i][j]);
        }
        fprintf(fp2,"\n");
    }
    FinaliseRandom();
    fclose(fp2);
}