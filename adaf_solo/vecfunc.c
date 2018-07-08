#include <math.h>
#include "vecfunc.h"
#include "ssfunctions.h"
#include "rates.h"

void vecfunc(int n, double x[], double fvec[])
{
    extern double r,betapar,f;
    double mdot=x[3];
    double tempi=x[1];
    double tempe=x[2];
    double rad=log10(r);
    double qie=qiefunc(tempi,tempe,mdot);
    double qemi=qem(tempe,mdot);
    double qplus=qp(mdot);
	
    fvec[1]=qie/qplus-(1.0-f);  //Eq 3.34
    fvec[2]=qie/qemi-1.0;  //Eq 3.35
	//fvec[2]=qemi/qplus-(1.0-f);

    fvec[3]=1700.14*tempi-(1039.44*betapar*c3()/r-tempe); //Eq. 2.16

    //fvec[3]=tempi-(1.08*tempe+6.66e12*betapar*c3()/r);
}
  

