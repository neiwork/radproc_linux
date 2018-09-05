#include "compton.h"
#include "probexact.h"
#include "probaprox.h"
#include <math.h>
#include <fmath/mathFunctions.h>
//#include <fmath/laguerre.h>
//#include <fparameters/nrutil.h>
extern "C" {
    #include <nrMath/nrutil.h>
    #include <nrMath/laguerre.h>
}
#include <fmath/physics.h>

/*double k(double om,double omprim,double normtemp,double ommin,double ommax)
{
    int N=6;
    double alf=0.0;
    void gaulag(double *x, double *w, int n, double alf);
    double *abscissas,*weights,*gamma,*c;
    abscissas=dvector(1,N);
    weights=dvector(1,N);
    gamma=dvector(1,N);
    c=dvector(1,N);
    double sum1=0.0,sum2=0.0;
    gaulag(abscissas, weights,N,alf);
    for (int i=1;i<=N;i++) {
        gamma[i]=abscissas[i]*normtemp+1.0;
        double beta=sqrt(1.0-1.0/(gamma[i]*gamma[i]));
        c[i]=weights[i]*gamma[i]*sqrt(gamma[i]*gamma[i]-1.0);
        if (om/omprim > (1.0-beta)/(1.0+2.0*beta*omprim/gamma[i]) && om/omprim < (1.0+beta)
        /(1.0-beta+2.0*omprim/gamma[i])) {
            //sum1+=c[i]*probaprox(om,omprim,gamma[i],ommin,ommax);//*rate(frecprim,gamma[i]);
            sum1 += c[i]*probexact(om,omprim,gamma[i]);
        } else {
            sum1 += 0.0;
        }
        sum2 += c[i];
    }
    return sum1/sum2;
} */

#include <iostream>
using namespace std;

double lumICin(Vector lumIn, Vector energies, double omprim, int nE) {
    double energy = electronMass*cLight2*omprim;
    int j=0;
    while (energy > energies[j] && j < nE-1) {
            j++;
    }
    if (j==nE-1) return 0.0;
    else {
        if (lumIn[j-1] > 0.0 && lumIn[j] > 0.0) {
            double alpha = log10(lumIn[j]/lumIn[j-1])*(log10(energy/energies[j])/log10(energies[j]/energies[j-1]));
            double lum = lumIn[j-1]*pow(10.0,alpha);
            return lum;
        } else {
            return 0.0;
        }
    }
}

double compton(Vector lumIn,Vector energies,double temp,int nE, int jE) {
    
    double om=energies[jE]/(electronMass*cLight2);
    double normtemp=boltzmann*temp/electronMass/cLight2;
    int nom=50;
    int N=6;
    double alf=0.0;
    void gaulag(double *x, double *w, int n, double alf);
    double *abscissas,*weights,*gamma,*c;
    abscissas=dvector(1,N);
    weights=dvector(1,N);
    gamma=dvector(1,N);
    c=dvector(1,N);
    double sum1,sum2;
    sum1=sum2=0.0;
    gaulag(abscissas, weights,N,alf);
    for (int i=1;i<=N;i++) {
        gamma[i]=abscissas[i]*normtemp+1.0;
        double beta=sqrt(1.0-1.0/(gamma[i]*gamma[i]));
        c[i]=weights[i]*gamma[i]*sqrt(gamma[i]*gamma[i]-1.0);
        double omprimmax = om/(1.0-beta*(1.0+2.0*om/gamma[i]));
        double omprimmin = om*(1.0-beta)/(1.0+beta-2.0*om/gamma[i]);
        if (omprimmin > 0.0 && omprimmax > omprimmin) {
            double omratio = omprimmax/omprimmin;
            //if (omratio < 1.0) return 0.0;
            double var_int=pow(omratio,1.0/nom);
            double omprim=omprimmin;//(omprimmin>0.0 ? omprimmin : 0.1);
            double ka=0.0;
            for (int j=1;j<nom;j++) {
                omprim=omprim*var_int;
                double domprim=omprim*(var_int-1.0);
                double lumICin_v = lumICin(lumIn,energies,omprim,nE);
                ka += lumICin_v*probexact(om,omprim,gamma[i])*domprim;
            }
            sum1 += c[i]*ka;
            sum2 += c[i];
        } else {
            sum1 = 0.0;
            sum2 = 1.0;
        }
    }
    return sum1/sum2;
}

void compton2(Vector& lumOut, double lumIn, Vector energies, double temp, int nE, int jE) {
    
    double omprim=energies[jE]/(electronMass*cLight2);
    double var_int=pow(energies[nE-1]/energies[0],1.0/nE);
    double dnuprim=energies[jE]/planck * (var_int-1.0);
    double normtemp=boltzmann*temp/(electronMass*cLight2);
    int nom=100;
    int N=6;
    double alf=0.0;
    void gaulag(double *x, double *w, int n, double alf);
    double *abscissas,*weights,*gamma,*c;
    abscissas=dvector(1,N);
    weights=dvector(1,N);
    gamma=dvector(1,N);
    c=dvector(1,N);
    double sum2=0.0;
    gaulag(abscissas, weights,N,alf);
    Vector lumOut1(nE,0.0);
    for (int i=1;i<=N;i++) {
        gamma[i]=abscissas[i]*normtemp+1.0;
        double beta=sqrt(1.0-1.0/(gamma[i]*gamma[i]));
        c[i]=weights[i]*gamma[i]*sqrt(gamma[i]*gamma[i]-1.0);
        double ommin = omprim*(1.0-beta)/(1.0+2.0*beta*omprim/gamma[i]);
        double eps = omprim/gamma[i];
        double aux = beta/(1.0+gamma[i]*(1.0+beta));
        double ommax = (eps < aux ? omprim*(1.0+beta)/(1.0-beta+2.0*omprim/gamma[i]) 
                                        : omprim+(gamma[i]-1.0));
        double var_int_om=pow(ommax/ommin,1.0/nom);
        double om=ommin;
        for (int jjE=1;jjE<nom;jjE++) {
            double lumOutAux=0.0;
            om = om*var_int_om;
            double dom=om*(var_int_om-1.0);
            double lim1=sqrt(energies[1]*energies[0]);
            double lim2;
            double energy = om*(electronMass*cLight2);
            int count=0;
            lumOutAux = lumIn * probexact(om,omprim,gamma[i]) * c[i];
            int jjjE=1;
            while(count == 0 && jjjE < nE-2) {
                lim2=sqrt(energies[jjjE]*energies[jjjE+1]);
                if ( energy < lim2 && energy > lim1) {
                    double dnu=energies[jjjE]/planck * (var_int-1.0);
                    lumOut1[jjjE] += lumOutAux * dom/dnu;
                    count=1;
                } else {
                    lim1=lim2;
                    jjjE++;
                }
            }
        }
        sum2 += c[i];
    }
    
    for (int jjE=0;jjE<nE;jjE++) {
        lumOut1[jjE] *= dnuprim/sum2;
        lumOut[jjE] += lumOut1[jjE];
    }
}


/* double compton2(Vector lumIn,Vector energyprim,double energy,double temp,int nE) {
    double omprim,domprim;
    double om=energy/(electronMass*cLight2);
    double normtemp=boltzmann*temp/electronMass/cLight2;
    double ommax=energyprim[nE-1]/(electronMass*cLight2);
    double ommin=energyprim[0]/(electronMass*cLight2);
    double var_int=pow(ommax/ommin,1.0/(nE-1));
    double sum=0.0;
    double omaux;
    omprim=ommin;
    for (int i=0;i<nE-1;i++) {
        omprim=energyprim[i]/(electronMass*cLight2);
        domprim=omprim*(var_int-1.0);
        double ka=k(om,omprim,normtemp,ommin,ommax);
        sum += lumIn[i]*ka*domprim;
        double sum2=0.0;
        for (int j=0;j<nE;j++) {
            omaux=energyprim[j]/(electronMass*cLight2);
            double dom=omaux*(var_int-1.0);
            sum2 += probability(omaux,omprim,3.0,ommin,ommax)*dom;
        }
        cout << sum2 << endl;
    }
    return sum;
} */