#include <math.h>
#include "vecfunc2.h"
#include "ssfunctions2.h"
#include "rates2.h"

void vecfunc2(int n, double x[], double fvec[])
{
    extern double r,betapar;
    double f=x[3];
    double tempi=x[1];
    double tempe=x[2];
    double rad=log10(r);
    double qie=qiefunc2(tempi,tempe,f);
    double qemi=qem2(tempe,f);
    double qplus=qp2(f);
  
    fvec[1]=qie/qplus-(1.0-f);
    fvec[2]=qie/qemi-1.0;
    //fvec[3]=1700.14*tempi-(1039.44*betapar*c3(f)/r+tempe);
    fvec[3]=tempi-(1.08*tempe+6.66e12*betapar*c3(f)/r);
}