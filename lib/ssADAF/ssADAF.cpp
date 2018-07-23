#include "ssADAF.h"
#include "vecfunc.h"
#include <fmath/constants.h>
extern "C" {
	#include <nrMath/nr.h>
}

void adafSol(double *x)
{	
	int check;
	broydn(x,3,&check,vecfunc);
}