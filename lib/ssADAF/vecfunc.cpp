#include <math.h>
#include "vecfunc.h"
#include "ssfunctions.h"
#include "rates.h"


#include <boost/property_tree/ptree.hpp>


void vecfunc(int n, double x[], double fvec[] //)
		, double r, dataADAF data)
{
    //extern double r,betapar,f;
	
	double massBH = data.massBH;
	double f = data.f;
	double alphapar = data.alphapar;
	double betapar  = data.betapar;
	double rmin = data.rmin;
	double rmax = data.rmax;
	double eDensity = data.eDensity;
	double iDensity = data.iDensity;
	double magField = data.magField;
	double eTemp = data.eTemp; 
	
	
/*	static const double massBH=GlobalConfig.get<double>("massBH")*solarMass;
	static const double f=GlobalConfig.get<double>("f");
	static const double alphapar=GlobalConfig.get<double>("alpha");
	static const double betapar=GlobalConfig.get<double>("beta");
	static const double rmin=GlobalConfig.get<double>("rmin");
	static const double rmax=GlobalConfig.get<double>("rmax");*/
    double gammapar=(32.0-24.0*betapar-3.0*betapar*betapar)/(24.0-21.0*betapar);
	
	//falta r, eDensity, iDensity, magField, eTemp,
	
    double mdot=x[3];
    double tempi=x[1];
    double tempe=x[2];
    double rad=log10(r);
	
	double qie=qiefunc(tempi, tempe, mdot,r, massBH, alphapar, gammapar, f);
	
	double qemi=qem(tempe, mdot, eDensity, iDensity, magField, eTemp,
					r, massBH, betapar, alphapar, gammapar, f); 
	
	double qplus=qp(mdot, r, massBH, alphapar, gammapar, f);
	
	
   // double qie=qiefunc(tempi,tempe,mdot);
    //double qemi=qem(tempe,mdot);
    //double qplus=qp(mdot);
	
    fvec[1]=qie/qplus-(1.0-f);  //Eq 3.34
    fvec[2]=qie/qemi-1.0;  //Eq 3.35
	//fvec[2]=qemi/qplus-(1.0-f);
    fvec[3]=1700.14*tempi-(1039.44*betapar*c3(alphapar, gammapar, f)/r-tempe);  //Eq. 2.16
    //fvec[3]=tempi-(1.08*tempe+6.66e12*betapar*c3()/r);
}
  

	

