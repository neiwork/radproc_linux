#pragma once
#include <fmath/physics.h>
#include <fstream>

double compton(Vector,Vector,double,int,int);
void compton2(Matrix&, double, Vector, double, int, int, int);
void compton3(Matrix&, double, Vector, double, int, int, int);