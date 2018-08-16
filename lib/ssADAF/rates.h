#pragma once

//double qem(double,double);
//void qsync(double *, double *, double);
//double qbremss(double);
//void qC(double, double, double, double, double *, double *,double);
//double qiefunc(double,double,double);


double qbremss(double tempe, double eDensity, double iDensity);

void qsync(double *qSy, double *nucrit, double tempe,
			double radius, double eDensity, double magField, double eTemp);
			
void qC(double qbr, double qsy, double nuc, double tempe, double *qCbr, double *qCsy, double mdot,
		double r, double alphapar, double gammapar, double f);
		
double qiefunc(double tempi, double tempe, double mdot,
		double r, double massBH, double alphapar, double gammapar, double f);
		
double qem(double tempe, double mdot,
		double eDensity, double iDensity, double magField, double eTemp,
		double r, double massBH, double betapar, double alphapar, double gammapar, double f); 