#pragma once

void constants();

/*double c1();
double c2();
double c3();

double rvel();              // Radial velocity.
double epsp();
double angvel();            // Angular velocity.
double sqrdsoundvel();      // Squared sound speed.
double massdens(double);          // Mass density.
double press(double);             // Total pressure.
double magf(double);              // Magnetic Field.
double qp(double);                // Viscous dissipation per unit volume.
double ne(double);                // Electron density.
double ni(double);                // Ion density.
double taues(double);             // Scattering optical depth.
double height();*/

//double mdotcrit();


double epsp(double gammapar, double f);

//double g(double alphapar, double gammapar, double f) {

double c1(double alphapar, double gammapar, double f);

double c2(double alphapar, double gammapar, double f);

double c3(double alphapar, double gammapar, double f);

//double radialvel(double r, double alphapar, double gammapar, double f);  
//double rvel();  	 ver nombre

//double angularvel(double r, double massBH, double alphapar, double gammapar, double f);
//double angvel();  idem

double sqrdsoundvel(double r, double alphapar, double gammapar, double f);

double massdens(double mdot, double r, double massBH, double alphapar, double gammapar, double f);

double press(double mdot, double r, double massBH, double alphapar, double gammapar, double f);

double magf(double mdot, double r, double massBH, double betapar, double alphapar, double gammapar, double f);

double qp(double mdot, double r, double massBH, double alphapar, double gammapar, double f);

double ne(double mdot, double r, double massBH, double alphapar, double gammapar, double f);


double ni(double mdot, double r, double massBH, double alphapar, double gammapar, double f);

double taues(double mdot, double r, double alphapar, double gammapar, double f);

double height(double radius, double alphapar, double gammapar, double f);
