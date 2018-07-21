#pragma once



//double rCusp, rCenter, l_0;

/*double w(double, double);
double edge();
double electronDensity(double, double);*/

void torusParameters();

double energyDensity(double r, double theta);

double electronDensity(double r, double theta);

double pressureTot(double r, double theta);

double temp_e(double r, double theta);

double temp_i(double r, double theta);