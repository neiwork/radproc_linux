#pragma once


#include <fparticle/Particle.h>

/*vel is the lateral jet expansion*/
double adiabaticLosses(double E, double z, double vel_lat, double gamma);  //en [erg/s]

/*Diffusion rate in the Bohm regime*/
double diffusionTimeTurbulence(double E, double h, Particle& p, double B);   //en [s]^-1
double diffusionTimeParallel(double E, double h, double B);   //en [s]^-1

double diffCoeff_r(double g, Particle& p, double height, double B);   // en cm^2 s^-1
double diffLength(double g, Particle& p, double r, double height, double B, double vR);   // en cm

double diffCoeff_p(double E, Particle& p, double height, double B, double rho);
double diffCoeff_g(double g, Particle& p, double height, double B, double rho);
double accelerationRate(double E, double magneticField); //en [s]^-1
double accelerationTimeSDA(double,Particle&,double,double,double);
double accelerationRateSDA(double E, Particle& p, double B, double height, double rho);
double escapeRate(double size, double vel);  //en [s]^-1