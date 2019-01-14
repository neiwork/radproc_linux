#pragma once
#include <fmath/physics.h>

void comptonNew(Matrix&, double, Vector, double, size_t, size_t, size_t);
void comptonNew2(Matrix&, double, Vector, size_t, size_t, size_t, double, Vector,size_t, size_t);
void comptonRedistribution(Vector&, size_t, size_t, size_t, size_t, double, Vector, Matrix);