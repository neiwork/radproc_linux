#pragma once

#include <iostream>
#include "mathematics.h"

void TriDiagSys(Vector a, Vector b, Vector c, Vector& d, int n);
void TriBlockDiagSys(Vector a, Vector ba, Vector bb, Vector bc, Vector c, Vector& d, size_t J, size_t M);
void GaussSeidel(Matrix& a, Vector& b, Vector& x, int n, int init, double &err);