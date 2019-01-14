#include "compton.h"
#include "probexact.h"
#include "probaprox.h"

#include "../ionTori/maxwellJuttner.h"
#include <math.h>
#include <fmath/mathFunctions.h>
#include <fmath/RungeKutta.h>
//#include <fmath/laguerre.h>
//#include <fparameters/nrutil.h>
extern "C" {
    #include <nrMath/nrutil.h>
    #include <nrMath/laguerre.h>
}
#include <fmath/physics.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>

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

/*double lumICin(Vector lumIn, Vector energies, double omprim, int nE) {
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
}*/

/*double compton(Vector lumIn,Vector energies,double temp,int nE, int jE) {
    
    double om=energies[jE]/(electronMass*cLight2);
    double normtemp=boltzmann*temp/electronMass/cLight2;
    int nOm=100;
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
        double omPrimMax = om/(1.0-beta*(1.0+2.0*om/gamma[i]));
        double omPrimMin = om*(1.0-beta)/(1.0+beta-2.0*om/gamma[i]);
        if (omPrimMin > 0.0 && omPrimMax > omPrimMin) {
            double omRatio = omPrimMax/omPrimMin;
            //if (omratio < 1.0) return 0.0;
            double var_int=pow(omRatio,1.0/nOm);
            double omPrim=omPrimMin;//(omprimmin>0.0 ? omprimmin : 0.1);
            double ka=0.0;
            for (int jOm=1;jOm<nOm;jOm++) {
                omPrim=omPrim*var_int;
                double dOmPrim=omPrim*(var_int-1.0);
                double lumICin_v = lumICin(lumIn,energies,omPrim,nE);
                ka += lumICin_v*probexact(om,omPrim,gamma[i])*dOmPrim;
            }
            sum1 += c[i]*ka;
        }
		sum2 += c[i];
    }
    return sum1/sum2;
}*/

void compton2(Matrix& lumOut, double lumIn, Vector energies, 
				double temp, int nE, int jEprim, int jR) 
{    
    double omPrim=energies[jEprim]/(electronMass*cLight2);
    double var_int=pow(energies[nE-1]/energies[0],1.0/nE);
    double dnuprim=energies[jEprim]/planck * (var_int-1.0);
    double normtemp=boltzmann*temp/(electronMass*cLight2);
    int nOm=100;
    int N=10;
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
        double omMin = omPrim*(1.0-beta)/(1.0+2.0*beta*omPrim/gamma[i]);
        double eps = omPrim/gamma[i];
        double aux = beta/(1.0+gamma[i]*(1.0+beta));
        double omMax = (eps < aux ? omPrim*(1.0+beta)/(1.0-beta+2.0*omPrim/gamma[i]) 
                                        : omPrim+(gamma[i]-1.0));
		Vector probab(nOm,0.0);
		
		double probtotexact = RungeKuttaSimple(omMin,omMax,[&](double om)
						{return probexact(om,omPrim,gamma[i]);});    // DEBERÍAN DAR 1
		
/*		for (int jOm=1;jOm<nOm;jOm++) {
			om *= var_int_om;
			double dOm=om*(var_int_om-1.0);
			probab[jOm]=probexact(om,omprim,gamma[i])/probtotexact;
			//probab[jOm]=probaprox(om,omprim,gamma[i],omMin,omMax);
			if (probab[jOm] < 0.0) cout << "ERROR" << endl;
			//lumOutAux2 += probab[jOm] * dOm;
		}*/
		
		// POR AHORA PROBAMOS CON ESTO
		
		
		double var_int_om=pow(omMax/omMin,1.0/nOm);
		//double var_int_om_aprox=pow(omMax_aprox/omMin_aprox,1.0/nOm);
		
		if (probtotexact > 0.0) {
		
		double om=omMin;
		int jE=1;
		for (int jOm=1;jOm<nOm;jOm++) {
			om *= var_int_om;
			double dOm=om*(var_int_om-1.0);
			probab[jOm]=probexact(om,omPrim,gamma[i])/probtotexact;
			if (probab[jOm] < 0.0) cout << "ERROR" << endl;
			
			double lim1=sqrt(energies[1]*energies[0]);
			double lim2;
			double energy = om * (electronMass*cLight2);
			
			double lumOutAux = lumIn * probab[jOm] * (om/omPrim);
			int count=0;
			while (count == 0 && jE < nE-2) {
				lim2 = sqrt(energies[jE]*energies[jE+1]);
				if ( energy < lim2 && energy > lim1) {
					double dnu = (lim2-lim1)/planck;
					lumOut1[jE] += lumOutAux * c[i] * dOm/dnu;
					count++;
				} else {
					lim1=lim2;
					jE++;
				}
			}
		}
		
		}
		
		sum2 += c[i];
	}

/*		om=omMin;
        for (int jOm=1;jOm<nOm;jOm++) {
            double lumOutAux=0.0;
            om = om*var_int_om;
            double dOm=om*(var_int_om-1.0);
            double lim1=sqrt(energies[1]*energies[0]);
            double lim2;
            double energy = om*(electronMass*cLight2);
            int count=0;
            lumOutAux = lumIn * probab[jOm] * c[i] * dOm;
			if (lumOutAux2 > 0.0) lumOutAux /= lumOutAux2;
            int jE=1;
            while(count == 0 && jE < nE-2) {
                lim2=sqrt(energies[jE]*energies[jE+1]);
                if ( energy < lim2 && energy > lim1) {
                    double dnu=(lim2-lim1)/planck;
                    lumOut1[jE] += lumOutAux /dnu;
                    count++;
                } else {
                    lim1=lim2;
                    jE++;
                }
            }
        }
        sum2 += c[i];
    }*/
	
    for (int jE=0;jE<nE;jE++) {
        lumOut[jE][jR] += lumOut1[jE] * dnuprim /sum2;
    }
}

/*void compton3(Matrix& lumOut, double lumIn, Vector energies, 
				double temp, int nE, int jEprim, int jR) 
{    
    double omPrim=energies[jEprim]/(electronMass*cLight2);
    double var_int=pow(energies[nE-1]/energies[0],1.0/nE);
    double dnuprim=energies[jEprim]/planck * (var_int-1.0);
    double normtemp=boltzmann*temp/(electronMass*cLight2);
    int nOm=100;
    int N=10;
	double gamma_min = 1.0;
	double gamma_max = 10.0*normtemp;
	double var_int_g = pow(gamma_max/gamma_min,1.0/N);
	double gamma = gamma_min;
    Vector lumOut1(nE,0.0);
    for (int i=1;i<=N;i++) {
        gamma *= var_int_g;
        double beta=sqrt(1.0-1.0/(gamma*gamma));
		double ngamma = maxwellRel(gamma,normtemp,1.0);
        double omMin = omPrim*(1.0-beta)/(1.0+2.0*beta*omPrim/gamma);
        double eps = omPrim/gamma;
        double aux = beta/(1.0+gamma*(1.0+beta));
        double omMax = (eps < aux ? omPrim*(1.0+beta)/(1.0-beta+2.0*omPrim/gamma) 
                                        : omPrim+(gamma-1.0));
		Vector probab(nOm,0.0);
		
		double probtotexact = RungeKuttaSimple(omMin,omMax,[&](double om)
						{return probexact(om,omPrim,gamma);});    // DEBERÍAN DAR 1
		
		double var_int_om=pow(omMax/omMin,1.0/nOm);
		
		if (probtotexact > 0.0) {
		
		double om=omMin;
		int jE=1;
		for (int jOm=1;jOm<nOm;jOm++) {
			om *= var_int_om;
			double dOm=om*(var_int_om-1.0);
			probab[jOm]=probexact(om,omPrim,gamma)/probtotexact;
			if (probab[jOm] < 0.0) cout << "ERROR" << endl;
			
			double lim1=sqrt(energies[1]*energies[0]);
			double lim2;
			double energy = om * (electronMass*cLight2);
			double lumOutAux = lumIn  * probab[jOm];
			int count=0;
			while (count == 0 && jE < nE-2) {
				lim2 = sqrt(energies[jE]*energies[jE+1]);
				if ( energy < lim2 && energy > lim1) {
					double dnu = (lim2-lim1)/planck;
					lumOut1[jE] += lumOutAux * ngamma * dOm/dnu;
					count++;
				} else {
					lim1=lim2;
					jE++;
				}
			}
		}
		
		}
	}
	
    for (int jE=0;jE<nE;jE++) {
        lumOut[jE][jR] += lumOut1[jE] * dnuprim;
    }
}*/


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