#pragma once
#include <fmath/physics.h>

double extrinf(double, void *);
double extrsup(double, void *);
double gslprobexact(double, void *);
double probexact (double,double,double);
double probexactNew (double,double,double);
double rate(double,double);
double rateThermal(double,double);
double probTemp(Vector,double,double,double);
double probInterpolated(Vector,double,double,double);
double probInterpolated2(Vector,double,double,double);