#include "compton.h"
#include "coppiProbDist.h"

double compton(Vector lumIn,Vector energyprim,double frec,double temp,int nE,double Emax,double Emin) {
    double frecprim,dfrecprim;
    double normtemp=boltzmann*temp/electronMass/cLight2;
    double frecmax=Emax/planck;
    double frecmin=Emin/planck;
    double var_int=pow(frecmax/frecmin,1.0/(nE-1));
    double sum=0.0;
    for (int i=0;i<nE;i++) {
        frecprim=energyprim[i]/planck;
        dfrecprim=frecprim*(var_int-1.0);
        sum+=lumIn[i]*k(frec,frecprim,normtemp)*dfrecprim;
    }
    return sum;
}