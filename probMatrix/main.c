#include <stdio.h>
#include "functions.h"
#include "torusSampling2.h"

double *rCells, *r, *ne;
double paso;
double rMin=3.0;
double rMax=1.0e4;
int nR=10;
int nTheta=1;
long nPhot=10;

int main()
{
	initVectors();
	torusSampling2();
	return 0;
}
