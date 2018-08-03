#include "ssADAF.h"
#include "vecfunc.h"


#include <fmath/constants.h>
extern "C" {
	#include <nrMath/nr.h>
}

void adafSol(double *x, double r, dataADAF data)
{	
	int check;
	//broydn(x,3,&check,vecfunc);
	broydn(x,3,&check, [&r, &data](int n, double x[], double fvec[]) { return  vecfunc(n, x[], fvec[], r, data); });
}