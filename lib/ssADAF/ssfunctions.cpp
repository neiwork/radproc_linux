
#include "ssfunctions.h"
#include <boost/property_tree/ptree.hpp>
#include <fparameters/parameters.h>
//#include <fmath/configure.h>
//#include <fparameters/parameters.h>
//#include <fmath/configure.h>
#include <math.h>

//static double alphapar,massBH,f,betapar,gammapar,rmin,rmax;
 
void constants()
{
    extern double alphapar,massBH,f,betapar,gammapar,rmin,rmax;
	
	massBH=GlobalConfig.get<double>("massBH");
	f=GlobalConfig.get<double>("f");
	alphapar=GlobalConfig.get<double>("alpha");
	betapar=GlobalConfig.get<double>("beta");
	rmin=GlobalConfig.get<double>("rmin");
	rmax=GlobalConfig.get<double>("rmax");
    gammapar=(32.0-24.0*betapar-3.0*betapar*betapar)/(24.0-21.0*betapar);
}

//double epsp(double f) {
double epsp() { //Eq 2.2
    extern double gammapar,f;
    return (5.0/3.0-gammapar)/(gammapar-1.0)/f;   
}

//Eq 2.2
double g() {
    extern double alphapar;
    double e=epsp();
    return sqrt(1.0+18.0*alphapar*alphapar/((5.0+2.0*e)*(5.0+2.0*e)))-1.0;  
}

//c1,c2,c3 Eq 2.1
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

//Eq 2.15
double radialvel() {
    extern double alphapar,r;
    return -2.12e10*alphapar*c1()/sqrt(r);          // [cm s^-1]
}

double angularvel() {  //Eq 2.15
    extern double r,massBH;
    return 7.19e4*c2()/(massBH*pow(r,1.5));           // [rad s^-1]
}

double sqrdsoundvel() {
    extern double r;
    return 4.5e20*c3()/r;                        // [cm^2 s^-2]
}

double massdens(double mdot) {
    extern double r,massBH,alphapar;
    return 3.79e-5*mdot/(m*alphapar*c1()*sqrt(c3())*pow(r,1.5));       // [g cm^-3]
}

double press(double mdot) {
    extern double r,massBH,alphapar;
    return 1.71e16*mdot*sqrt(c3())/(alphapar*c1()*massBH*pow(r,2.5));       // [g cm^-1 s^-2]
}

double magf(double mdot)   {
    extern double r,alphapar,betapar,massBH;
    return 6.55e8*sqrt((1.0-betapar)*sqrt(c3())*mdot/(c1()*m*alphapar))/pow(r,1.25);    // [G]
}

double qp(double mdot) {
    extern double r,massBH;
    return 1.84e21*epsp()*sqrt(c3())*mdot/(massBH*massBH*r*r*r*r);            // [erg cm^-3 s^-1]
}

double ne(double mdot) {
    extern double r,alphapar,massBH;
    return 2.0e19*mdot/(alphapar*c1()*massBH*sqrt(c3()*r*r*r));             // [cm^-3]
}

double ni(double mdot) {
    return 0.9268*ne(mdot);  //ni=ne*mu_e/mu_i
}

double taues(double mdot) {
    extern double r,alphapar;
    return 12.4*mdot/(alphapar*c1()*sqrt(r));
}
//hasta aca Eq 2.15

//Eq 2.4
double height() {
    extern double radius;
    return radius*sqrt(2.5*c3());
}