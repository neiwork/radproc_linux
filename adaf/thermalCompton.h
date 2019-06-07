#pragma once
#include <fmath/physics.h>

void comptonNew(Matrix&, double, Vector, double, size_t, size_t, size_t);
void comptonMatrix();
void comptonMatrix2();
void comptonNew2(Matrix&, double, Vector, size_t, size_t, size_t, double, Vector,size_t, size_t);
void comptonRedistribution(Vector&, size_t, size_t, size_t, size_t, double, Vector, Matrix);
double comptonNewNew(Vector, Vector, double, double,Vector,size_t);
double comptonNewLocal(Vector,Vector,double,size_t,Vector,size_t);