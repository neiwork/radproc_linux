#include <math.h>
#include "ssfunctions.h"

void constants()
{
    extern double alphapar,m,f,RS,betapar,gammapar,rmin,rmax,step;
    extern int mgrid;
  
    m=1.0e8;                                        // [Solar masses]
    //mdot=1.0e-4;                                // [Eddington luminosities]
    f=0.9;
    alphapar=0.3;
    RS=m*2.9675e5;                          // Schwarzschild radius
    betapar=0.5;
    gammapar=(32.0-24.0*betapar-3.0*betapar*betapar)/(24.0-21.0*betapar);
    mgrid=100;
    rmin=1.0;                                        // [Schwwarzschild radii]
    rmax=1.0e6;                                  //                     "
    step=pow(rmax/rmin,1.0/mgrid);
}

//double epsp(double f) {
double epsp() {
    extern double gammapar,f;
    return (5.0/3.0-gammapar)/(gammapar-1.0)/f;
}

double g() {
    extern double alphapar;
    double e=epsp();
    return sqrt(1.0+18.0*alphapar*alphapar/((5.0+2.0*e)*(5.0+2.0*e)))-1.0;
}

double c1() {
    extern double alphapar;
    return (5.0+2.0*epsp())/(3.0*alphapar*alphapar)*g();
}

double c2() {
    extern double alphapar;
    double e=epsp();
    return sqrt(2.0*e*(5.0+2.0*e)*g()/(9.0*alphapar*alphapar));
}

double c3() {
    extern double alphapar;
    return 2.0*(5.0+2.0*epsp())*g()/(9.0*alphapar*alphapar);
}

double rvel() {
    extern double alphapar,r;
    return -2.12e10*alphapar*c1()/sqrt(r);          // [cm s^-1]
}

double angvel() {
    extern double r,m;
    return 7.19e4*c2()/(m*pow(r,1.5));           // [rad s^-1]
}

double sqrdsoundvel() {
    extern double r;
    return 4.5e20*c3()/r;                        // [cm^2 s^-2]
}

double massdens(double mdot) {
    extern double r,m,alphapar;
    return 3.79e-5*mdot/(m*alphapar*c1()*sqrt(c3())*pow(r,1.5));       // [g cm^-3]
}

double press(double mdot) {
    extern double r,m,alphapar;
    return 1.71e16*mdot*sqrt(c3())/(alphapar*c1()*m*pow(r,2.5));       // [g cm^-1 s^-2]
}

double magf(double mdot)   {
    extern double r,alphapar,betapar,m;
    return 6.55e8*sqrt((1.0-betapar)*sqrt(c3())*mdot/(c1()*m*alphapar))/pow(r,1.25);    // [G]
}

double qp(double mdot) {
    extern double r,m,m;
    return 1.84e21*epsp()*sqrt(c3())*mdot/(m*m*r*r*r*r);            // [erg cm^-3 s^-1]
}

double ne(double mdot) {
    extern double r,alphapar,m;
    return 2.0e19*mdot/(alphapar*c1()*m*sqrt(c3()*r*r*r));             // [cm^-3]
}

double ni(double mdot) {
    return 0.9268*ne(mdot);
}

double taues(double mdot) {
    extern double r,alphapar;
    return 12.4*mdot/(alphapar*c1()*sqrt(r));
}

double height() {
    extern double radius;
    return radius*sqrt(2.5*c3());
}