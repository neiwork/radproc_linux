#include <math.h>
#include "ssfunctions.h"
#include <fmath/constants.h>
//#include <fmath/configure.h>
#include <fparameters/parameters.h>
//#include <fmath/configure.h>
//#include <boost/property_tree/ptree.hpp>

/*void constants()
{
    extern double alphapar,massBH,f,betapar,gammapar,rmin,rmax;
	
	massBH=GlobalConfig.get<double>("massBH");
	f=GlobalConfig.get<double>("f");
	alphapar=GlobalConfig.get<double>("alpha");
	betapar=GlobalConfig.get<double>("beta");
	rmin=GlobalConfig.get<double>("rmin");
	rmax=GlobalConfig.get<double>("rmax");
    gammapar=(32.0-24.0*betapar-3.0*betapar*betapar)/(24.0-21.0*betapar);
}*/

//Eq 2.2
//double epsp() { 
//    extern double gammapar,f;
double epsp(double gammapar, double f) {
    return (5.0/3.0-gammapar)/(gammapar-1.0)/f;   
}

//Eq 2.2
//double g() {
//    extern double alphapar;
double g(double alphapar, double gammapar, double f) {
    double e=epsp(gammapar,f);
    return sqrt(1.0+18.0*alphapar*alphapar/((5.0+2.0*e)*(5.0+2.0*e)))-1.0;  
}

//c1,c2,c3 Eq 2.1
//double c1() {
//    extern double alphapar;
double c1(double alphapar, double gammapar, double f) {
    return (5.0+2.0*epsp(gammapar,f))/(3.0*alphapar*alphapar)*g(alphapar,gammapar,f);
}

//double c2() {
//   extern double alphapar;
double c2(double alphapar, double gammapar, double f) {
    double e=epsp(gammapar,f);
    return sqrt(2.0*e*(5.0+2.0*e)*g(alphapar,gammapar,f)/(9.0*alphapar*alphapar));
}

//double c3() {
//    extern double alphapar;
double c3(double alphapar, double gammapar, double f) {
    return 2.0*(5.0+2.0*epsp(gammapar,f))*g(alphapar,gammapar,f)/(9.0*alphapar*alphapar);
}

//Eq 2.15
//double radialvel() {
//    extern double alphapar,r;
double radialvel(double r, double alphapar, double gammapar, double f) {	
    return -2.12e10*alphapar*c1(alphapar,gammapar,f)/sqrt(r);          // [cm s^-1]  
}

//double angularvel() {  //Eq 2.15
//    extern double r,massBH;
double angularvel(double r, double massBH, double alphapar, double gammapar, double f) {
    return 7.19e4*c2(alphapar,gammapar,f)/(massBH*pow(r,1.5));           // [rad s^-1]
}

//double sqrdsoundvel() {
//    extern double r;
double sqrdsoundvel(double r, double alphapar, double gammapar, double f){
    return 4.5e20*c3(alphapar,gammapar,f)/r;                        // [cm^2 s^-2]
}

//double massdens(double mdot) {
//    extern double r,massBH,alphapar;
double massdens(double mdot, double r, double massBH, double alphapar, double gammapar, double f) {
	double m = massBH/solarMass;  //VER
    return 3.79e-5*mdot/(m*alphapar*c1(alphapar,gammapar,f)*sqrt(c3(alphapar,gammapar,f))*pow(r,1.5));       // [g cm^-3]
}

//double press(double mdot) {
//    extern double r,massBH,alphapar;
double press(double mdot, double r, double massBH, double alphapar, double gammapar, double f) {
    return 1.71e16*mdot*sqrt(c3(alphapar,gammapar,f))/(alphapar*c1(alphapar,gammapar,f)*massBH*pow(r,2.5));       // [g cm^-1 s^-2]
}

//double magf(double mdot)   {
//    extern double r,alphapar,betapar,massBH;
double magf(double mdot, double r, double massBH, double betapar, double alphapar, double gammapar, double f) {
    return 6.55e8*sqrt((1.0-betapar)*sqrt(c3(alphapar,gammapar,f))*mdot/(c1(alphapar,gammapar,f)*m*alphapar))/pow(r,1.25);    // [G]
}

//double qp(double mdot) {
//    extern double r,massBH;
double qp(double mdot, double r, double massBH, double alphapar, double gammapar, double f) {
    return 1.84e21*epsp(gammapar,f)*sqrt(c3(alphapar,gammapar,f))*mdot/(massBH*massBH*r*r*r*r);            // [erg cm^-3 s^-1]
}

//double ne(double mdot) {
//    extern double r,alphapar,massBH;
double ne(double mdot, double r, double massBH, double alphapar, double gammapar, double f) {
    return 2.0e19*mdot/(alphapar*c1(alphapar,gammapar,f)*massBH*sqrt(c3(alphapar,gammapar,f)*r*r*r));             // [cm^-3]
}

double ni(double mdot, double r, double massBH, double alphapar, double gammapar, double f) {
    return 0.9268*ne(mdot, r, massBH, alphapar, gammapar, f);  //ni=ne*mu_e/mu_i
}

//double taues(double mdot) {
//    extern double r,alphapar;
double taues(double mdot, double r, double alphapar, double gammapar, double f) {
    return 12.4*mdot/(alphapar*c1(alphapar,gammapar,f)*sqrt(r));
}
//hasta aca Eq 2.15

//Eq 2.4
//double height() {
//    extern double radius;
double height(double radius, double alphapar, double gammapar, double f) {
    return radius*sqrt(2.5*c3(alphapar,gammapar,f));
}