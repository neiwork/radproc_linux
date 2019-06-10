#pragma once
#include <fmath/physics.h>

void comptonMatrix();
void comptonMatrix2();
double comptonNewNew(Vector,Vector,Vector,Vector,Vector,double,double,Vector,size_t);
double comptonNewNewNew(Vector,double,double,Vector,size_t);
double comptonNewLocal(Vector,Vector,double,size_t,Vector,size_t);

double comptonNewNewPrueba(Vector p, Vector lumIn, size_t jTemp, Vector energies, size_t jE);
void comptonNewNewPruebaVector(Vector tempVec, Vector nuPrimVec, Vector nuVec, Vector probVec,
			size_t jTemp, Vector energies, size_t jE, Vector& p, double normtemp);
			
void comptonNewNewNewPruebaVector(size_t jR, Vector energies, size_t jE, Vector& p, double normtemp);