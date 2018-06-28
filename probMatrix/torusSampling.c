#include "torusSampling.h"
#include "auxFunctions.h"
#include "nrutil.h"
#include "torusParameters.h"
#include "random.h"
#include <stdio.h>

#define PI 3.14159265359
#define thomson 6.652458734e-25
#define RG 1.5e11
#define RANDOM_GENERATOR gsl_rng_r250 /* R250 random number generator from GNU Scientific Library. */

#define new_max(x,y) ((x) >= (y)) ? (x) : (y)
#define new_min(x,y) ((x) <= (y)) ? (x) : (y)

extern double rMin, rCenter, l_0;
extern double spinBH, lambda;

void torusSampling()
{
    //criticalRadii(&rMin, &rCenter);
	rMin=1.0;
    //l_0=specificAngularMom();
    //double rMax=edge();
	double rMax=1.0e4;
    int nR=100,nTheta=10;
    long nPhot=1000;
    FILE *fp1, *fp2;
    
    fp1=fopen("torus.txt","w");
    fp2=fopen("prob.txt","w");
    
    //double dr=(rMax-rMin)/nR;                                // Cells' radial size.
	double paso=pow(rMax/rMin,1.0/nR);
	double *rCells,*r,**prob;
    
    rCells=dvector(0,nR);                                          // Cells' boundaries.
    r=dvector(1,nR);                                                   // Cells' central position.
    prob=dmatrix(1,nR,1,nR);
    rCells[0]=rMin;
    /*for(int i=1;i<=nR;i++) {
        rCells[i]=rCells[i-1]+dr;
        r[i]=rCells[i-1]+dr/2.0;
    }*/
	for(int=1;i<=nR;i++) {
		rCells[i]=rMin*paso;
		r[i]=sqrt(rCells[i]*rCells[i-1]);
	}
    double thetaMin=0.0;
    double thetaMax=PI/4.0;
    
    InitialiseRandom(RANDOM_GENERATOR);
    for(int i=1;i<=nR;i++) {
        double dyaux=(sin(thetaMax)-sin(thetaMin))/nTheta;
        for(int k=1;k<=nTheta+1;k++) {
            double yaux=sin(thetaMin)+(k-1)*dyaux;                  // Theta distributed uniformly in sin(theta).
            double theta0=asin(yaux);
            double y0=r[i]*cos(theta0);
            double z0=r[i]*sin(theta0);
            /*if (w(r[i],theta0) > 0) {
                fprintf(fp1,"%f   %f\n",y0,z0);
                fprintf(fp1,"%f   %f\n",y0,-z0);
                fprintf(fp1,"%f   %f\n",-y0,z0);
                fprintf(fp1,"%f   %f\n",-y0,-z0);
            }*/
            double drprim=dr/100.0;                                   // Initial step for the photon path.
            double drprimmin=dr/1.0e6;                                // Minimum step.
            double drprimmax=dr/10.0;                                 // Maximum step.
            for(int j=1;j<=nPhot;j++) {
                double random_number;
                random_number=gsl_rng_uniform(RandomNumberGenerator);
                double phiprim=2.0*PI*random_number;
                random_number=gsl_rng_uniform(RandomNumberGenerator);
                double thetaprim=acos((PI/2.0-random_number*PI)/(PI/2.0));        // Photon directions distributed
                double rprim=drprim;                                                                                       // isotropically.
                double rprimant=0.0;
                double accumulatedp=0.0;
                double neant, ne;
                double r1,theta1;
                neant=electronDensity(r[i],theta0);
                do {
                    int count=0;
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
                    double psc=ne*thomson*(drprim*RG);    // Probability of scattering.
                    double pescap=1.0-accumulatedp;           // Probability that a photon reaches the previous position.
                    for(int j=1;j<=nR;j++) {
                        if(r1 > rCells[j-1] && r1 < rCells[j]) prob[i][j] += psc*pescap;            // Add the probability of
                                                                                                                                                      // interaction in the matrix
                                                                                                                                                      // element.
                    }
                    printf("r1 = %f, drprim=%6.3e, j = %d, i = %d\n",r1,drprim,j,i);
                    accumulatedp+=psc;    // The product of (1-p_k) gives the probability for a photon to reach the
                                                              // previous position. To first order, this is 1-sum(p_k). Hence, we accumulate
                                                              // the sum of the p_k.
                    neant=ne;
                    rprimant=rprim;
                    drprim=new_min(drprim*2,drprimmax);
                    rprim+=drprim;
                } while(ne > 1.0e-10);                          // Escape from the torus.
            }
        }
        for(int j=1;j<=nR;j++) {
            prob[i][j] /= (nPhot*(nTheta+1));         // Dividing by the number of photons launched.
            fprintf(fp2,"%8.5e  ",prob[i][j]);
        }
        fprintf(fp2,"\n");
    }
    FinaliseRandom();
    fclose(fp1);
    fclose(fp2);
}