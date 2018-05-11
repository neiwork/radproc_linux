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

extern double rMin, rCenter, l_0;
extern double spinBH, lambda;

void torusSampling()
{
    criticalRadii(&rMin, &rCenter);
    l_0 = specificAngularMom();
    double rMax=edge();
    int nR=10,nTheta=3;
    long nPhot=1000;
    FILE *fp1, *fp2;
    
    fp1=fopen("torus.txt","w");
    fp2=fopen("prob.txt","w");
    
    double dr=(rMax-rMin)/nR;
    double *rCells,*r,**prob;
    
    rCells=dvector(0,nR);
    r=dvector(1,nR);
    prob=dmatrix(1,nR,1,nR);
    rCells[0]=rMin;
    for(int i=1;i<=nR;i++) {
        rCells[i]=rCells[i-1]+dr;
        r[i]=rCells[i-1]+dr/2.0;
    }
    double thetaMin=0.0;
    double thetaMax=PI/6.0;
    
    InitialiseRandom(RANDOM_GENERATOR);
    for(int i=1;i<=nR;i++) {
        double dyaux=(sin(thetaMax)-sin(thetaMin))/nTheta;
        for(int k=1;k<=nTheta+1;k++) {
            double yaux=sin(thetaMin)+(k-1)*dyaux;
            double theta0=asin(yaux);
            double y0=r[i]*cos(theta0);
            double z0=r[i]*sin(theta0);
            if (w(r[i],theta0) > 0) {
                fprintf(fp1,"%f   %f\n",y0,z0);
                fprintf(fp1,"%f   %f\n",y0,-z0);
                fprintf(fp1,"%f   %f\n",-y0,z0);
                fprintf(fp1,"%f   %f\n",-y0,-z0);
            }
            double drprim=dr/100.0;
            for(int j=1;j<=nPhot;j++) {
                double random_number;
                random_number=gsl_rng_uniform(RandomNumberGenerator);
                double phiprim=2.0*PI*random_number;
                random_number=gsl_rng_uniform(RandomNumberGenerator);
                double thetaprim=acos((PI/2.0-random_number*PI)/(PI/2.0));
                double rprim=drprim;
                double accumulatedp=0.0;
                double neant, ne;
                double r1,theta1;
                neant=electronDensity(r[i],theta0);
                do {
                    int count=0;
                    do {
                        if(count >= 1) {
                            drprim /= 2.0;
                            rprim -= drprim;
                        }
                        double xprim=rprim*sin(thetaprim)*cos(phiprim);
                        double yprim=rprim*sin(thetaprim)*sin(phiprim);
                        double zprim=rprim*cos(thetaprim);
                        double x1=xprim;
                        double y1=y0+yprim;
                        double z1=z0+zprim;
                        r1=sqrt(x1*x1+y1*y1+z1*z1);
                        theta1=atan(z1/sqrt(x1*x1+y1*y1));
                        ne=electronDensity(r1,abs(theta1));
                        count++;
                    } while(0);//while(abs(ne-neant)/(0.5*(ne+neant)) > 1.0e-2 || ne < 1.0e-9 || neant < 1.0e-9);
                    double psc=ne*thomson*(drprim*RG);
                    double pescant=1.0-accumulatedp;
                    for(int j=1;j<=nR;j++) {
                        if(r1 > rCells[j-1] && r1 < rCells[j]) prob[i][j] += psc*pescant;
                    }
                    accumulatedp+=psc;
                    neant=ne;
                    rprim+=drprim;
                    //drprim *= 2;
                } while(ne > 1.0);
                printf("%d\n",j);
            }
        }
        for(int j=1;j<=nR;j++) {
            prob[i][j] /= (nPhot*(nTheta+1));
            fprintf(fp2,"%f  ",prob[i][j]);
        }
        fprintf(fp2,"\n");
    }
    FinaliseRandom();
    fclose(fp1);
    fclose(fp2);
}