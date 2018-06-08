#include <math.h>
#include "vecfunc.h"
#include "ssfunctions.h"
#include "rates.h"

void vecfunc(int n, double x[], double fvec[])
{
    extern double r,betapar;
    double f=x[3];
    double tempi=x[1];
    double tempe=x[2];
    double rad=r;
    double qie=qiefunc(tempi,tempe,f);
    double qemi=qem(tempe,f);
    double qplus=qp(f);
  
    fvec[1]=qie/qplus-(1.0-f);
    fvec[2]=qie/qemi-1.0;
    fvec[3]=1700.14*tempi-(1039.44*betapar*c3(f)/r+tempe);
}
  

