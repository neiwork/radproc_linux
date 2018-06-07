#include <math.h>
#include "ssfunctions.h"

void constants()
{
    extern double alphapar,m,mdot,RS,betapar,gammapar,rmin,rmax,step;
    extern int mgrid;
  
    m=1.0e8;                                        // [Solar masses]
    mdot=1.0e-6;                                // [Eddington luminosities]
    alphapar=0.3;
    RS=m*2.9675e5;                          // Schwarzschild radius
    betapar=0.5;
    gammapar=(32.0-24.0*betapar-3.0*betapar*betapar)/(24.0-21.0*betapar);
    mgrid=100;
    rmin=1.0;                                        // [Schwwarzschild radii]
    rmax=1.0e5;                                  //                     "
    step=pow(rmax/rmin,1.0/mgrid);
}

double epsp(double f) {
    extern double gammapar;
    return (5.0/3.0-gammapar)/(gammapar-1.0)/f;
}

double g(double f) {
    extern double alphapar;
    double e=epsp(f);
    return sqrt(1.0+18.0*alphapar*alphapar/((5.0+2.0*e)*(5.0+2.0*e)))-1.0;
}

double c1(double f) {
    extern double alphapar;
    return (5.0+2.0*epsp(f))/(3.0*alphapar*alphapar)*g(f);
}

double c2(double f) {
    extern double alphapar;
    double e=epsp(f);
    return sqrt(2.0*e*(5.0+2.0*e)*g(f)/(9.0*alphapar*alphapar));
}

double c3(double f) {
    extern double alphapar;
    return 2.0*(5.0+2.0*epsp(f))*g(f)/(9.0*alphapar*alphapar);
}

double rvel(double f) {
    extern double alphapar,r;
    return -2.12e10*alphapar*c1(f)/sqrt(r);          // [cm s^-1]
}

double angvel(double f) {
    extern double r,m;
    return 7.19e4*c2(f)/(m*pow(r,1.5));           // [rad s^-1]
}

double sqrdsoundvel(double f) {
    extern double r;
    return 4.5e20*c3(f)/r;                        // [cm^2 s^-2]
}

double massdens(double f) {
    extern double r,mdot,m,alphapar;
    return 3.79e-5*mdot/(m*alphapar*c1(f)*sqrt(c3(f))*pow(r,1.5));       // [g cm^-3]
}

double press(double f) {
    extern double r,mdot,m,alphapar;
    return 1.71e16*mdot*sqrt(c3(f))/(alphapar*c1(f)*m*pow(r,2.5));       // [g cm^-1 s^-2]
}

double magf(double f)   {
    extern double r,alphapar,betapar,mdot,m;
    return 6.55e8*sqrt((1.0-betapar)*sqrt(c3(f))*mdot/(c1(f)*m*alphapar))/pow(r,1.25);    // [G]
}

double qp(double f) {
    extern double r,mdot,m,m;
    return 1.84e21*epsp(f)*sqrt(c3(f))*mdot/(m*m*r*r*r*r);            // [erg cm^-3 s^-1]
}

double ne(double f) {
    extern double r,mdot,alphapar,m;
    return 2.0e19*mdot/(alphapar*c1(f)*m*sqrt(c3(f)*r*r*r));             // [cm^-3]
}

double ni(double f) {
    return 0.9268*ne(f);
}

double taues(double f) {
    extern double r,mdot,alphapar;
    return 12.4*mdot/(alphapar*c1(f)*sqrt(r));
}