#pragma once
#include <stdio.h>
#include <iostream>
#include "physics.h"
#include <math.h>
extern "C" {
	#include "../nrMath/nrutil.h"
//#include <nrMath/nr.h>
}

double fbisection(fun1,double,double,double);
