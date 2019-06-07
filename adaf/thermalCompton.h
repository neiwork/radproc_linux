#pragma once
#include <fmath/physics.h>

void comptonMatrix();
void comptonMatrix2();
double comptonNewNew(Vector,Vector,Vector,Vector,Vector,double,double,Vector,size_t);
double comptonNewNewNew(Vector,double,double,Vector,size_t);
double comptonNewLocal(Vector,Vector,double,size_t,Vector,size_t);