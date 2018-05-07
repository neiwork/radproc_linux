#include "auxFunctions.h"
#include <math.h>

// Keplerian specific angular momentum
double keplAngularMom(double r)
{
  return sqrt(massBH) * ( r*r - 2.0 * spinBH * sqrt(massBH*r) + spinBH*spinBH ) /
    (pow(r, 1.5) - 2.0 * massBH * sqrt(r) + spinBH* sqrt(massBH) );
}

void marginalOrbits(double *r_ms, double *r_mb)
{
  double z1 = 1.0 + pow( 1.0 - (spinBH/massBH)*(spinBH/massBH), 1.0/3.0) * 
    ( pow(1.0 + spinBH/massBH, 1.0/3.0) + pow(1.0 - spinBH/massBH, 1.0/3.0) );
  double z2 = sqrt( 3.0 * (spinBH/massBH)*(spinBH/massBH) + z1*z1);
  
  // marginally stable circular orbit
  *r_ms = massBH * (3.0 + z2 - sqrt( (3.0 - z1) * (3.0 + z1 + 2.0*z2) ) );
  // marginally bound circular orbit
  *r_mb = 2.0 * massBH - spinBH + 2.0 * sqrt(massBH) * sqrt(massBH-spinBH);
}
double specificAngularMom()
{
  double r_ms, r_mb;
  marginalOrbits(&r_ms, &r_mb);
  
  double l_ms = keplAngularMom(r_ms);             // Keplerian specific angular momentum at r = r_ms
  double l_mb = keplAngularMom(r_mb);             // Keplerian specific angular momentum at r = r_mb
  
  return (1.0 - lambda) * l_ms + lambda * l_mb;
}

// AUXILIARY FUNCTIONS /////////////////////////////
double modfKepl(double r)
{
  double l_0 = specificAngularMom();
  return keplAngularMom(r) - l_0;
}
void criticalRadii(double *rCusp, double *rCenter)
{
  double r_ms, r_mb;
  marginalOrbits(&r_ms, &r_mb);

  int maxmitr = 1000;
  double allerr = 1.0e-3;
  *rCusp = bisection(r_mb, r_ms, allerr, maxmitr, modfKepl);
  *rCenter = bisection(r_ms, 500.0, allerr, maxmitr, modfKepl);
}
