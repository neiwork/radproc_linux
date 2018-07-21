//#include "ssADAF.h"


#include "rates.h"
#include "vecfunc.h"
#include "ssfunctions.h"
#include <fmath/constants.h>
extern "C" {
#include <nrMath/nrutil.h>
#include <nrMath/nr.h>
}
#include <stdio.h>
#include <math.h>

#define TOL 1.0e-20

/*#define SOLARMASS 1.998e33
#define PI 3.14159265359
#define BOLTZMANN 1.38e-16
#define ELECTRONMASS 9.109e-28
#define PROTONMASS 1.67262e-24
#define CLIGHT2 8.9875518e20*/

void adafSol(double *x)
{	
	broydn(x,N,&check,vecfunc);
}