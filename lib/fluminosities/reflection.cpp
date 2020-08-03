#include "reflection.h"
#include <fmath/physics.h>

double lumShellInc(double freq0, size_t jE, size_t jRcd, Matrix lumOut, Matrix absCD,
					Vector energies, size_t nR, Matrix redshift_RIAF_to_CD)
{
	double lum = 0.0;
	for (size_t jjR=0;jjR<nR;jjR++) {
		size_t jjE = jE;
		while (freq0*planck > energies[jjE]) jjE++;
		double logLumFreq0;
		if (jjE > 0) {
			if (lumOut[jjE-1][jjR] > 0.0 && lumOut[jjE][jjR] > 0.0) {
				logLumFreq0 = log10(lumOut[jjE][jjR]/lumOut[jjE-1][jjR]) /
					log10(energies[jjE]/energies[jjE-1]) * log10(freq0*planck/energies[jjE-1])
					+ log10(lumOut[jjE-1][jjR]);
				lum += pow(10.0,logLumFreq0) * absCD[jjR][jRcd] * pow(redshift_RIAF_to_CD[jjR][jRcd],2);
			}
		}
	}
	return lum;
}

double sigmanu(double x)
{
	double c0,c1,c2;
	double e = x*electronMass*cLight2/1.6e-9;
	
	if (e > 0.03 && e < 0.1) 	     { c0 = 17.3;  c1 = 608.1;  c2 = -2150;}
	else if (e > 0.1   && e < 0.284) { c0 = 34.6;  c1 = 267.9;  c2 = -476.1;}
	else if (e > 0.284 && e < 0.4)   { c0 = 78.1;  c1 = 18.8;   c2 = 4.3;}
	else if (e > 0.4   && e < 0.532) { c0 = 71.4;  c1 = 66.8;   c2 = -51.4;}
	else if (e > 0.532 && e < 0.707) { c0 = 95.5;  c1 = 145.8;  c2 = -61.1;}
	else if (e > 0.707 && e < 0.867) { c0 = 308.9; c1 = -380.6; c2 = 294.0;}
	else if (e > 0.867 && e < 1.303) { c0 = 120.6; c1 = 169.3;  c2 = -47.7;}
	else if (e > 1.303 && e < 1.84)  { c0 = 141.3; c1 = 146.8;  c2 = -31.5;}
	else if (e > 1.84  && e < 2.471) { c0 = 202.7; c1 = 104.7;  c2 = -17.0;}
	else if (e > 2.471 && e < 3.21)  { c0 = 342.7; c1 = 18.7;   c2 = 0.0;}
	else if (e > 3.21  && e < 4.038) { c0 = 352.2; c1 = 18.7;   c2 = 0.0;}
	else if (e > 4.038 && e < 7.111) { c0 = 433.9; c1 = -2.4;   c2 = 0.75;}
	else if (e > 7.111 && e < 8.331) { c0 = 629.0; c1 = 30.9;   c2 = 0.0;}
	else if (e > 8.331 && e < 10.0)  { c0 = 701.2; c1 = 25.2;   c2 = 0.0;}
	else return 0.0;
	
	return (c0+c1*e+c2*e*e)/(e*e*e) * 1.0e-24;
}

double greenDeltaFunc(double x)
{
	double knu = sigmanu(x)/protonMass;
	double kes = 0.4;
	double epsilon = knu / (knu+kes);
	return (1.0-sqrt(epsilon))/(1.0+sqrt(epsilon)) * planck/(electronMass*cLight2);
}

double greenFuncRefl(double freq, double freq0)
{
	double x = planck*freq / (electronMass*cLight2);
	double x0 = planck*freq0 / (electronMass*cLight2);
	double y = 1.0/x;
	double y0 = 1.0/x0;
	double deltay = y-y0;
	double deltayc = 1.0e3-y0;
	double w = exp(0.25e-5*(y0*y0*y0*y0 - y*y*y*y));
	if (w > 0.0) {
		double Ges;
		double A = 0.56 + 1.12*pow(y0,-0.785)-0.34*pow(y0,-1.04);
		double alpha = -0.3*pow(y0,-0.51)+0.06*pow(y0,-0.824);
		double beta = 0.37-pow(y0,0.85);
		if (deltay < 2.0) {
			double B = (1.0-A*(2.0+(pow(0.5*deltayc,0.5+alpha)-1.0)/(0.5+alpha))/sqrt(deltayc)) / 
				(pow(y0,1.0-beta)*pow(y0+2.0,beta)*(pow(1.0+2.0*x0,1.0-beta)-1.0)/(1.0-beta));
			if (abs(alpha+0.5) < 1.0e-6) {
				B = (1.0-A*(2.0+log(0.5*deltay))/sqrt(deltayc)) /
				(pow(y0,1.0-beta)*pow(y0+2.0,beta)*(pow(1.0+2.0*x0,1.0-beta)-1.0)/(1.0-beta));
			}
			Ges = B * pow((y0+2.0)/y,beta);
		} else if ( deltay >= 2.0 && deltay < deltayc) {
			Ges = A*pow(deltay,-1.5)*pow(deltayc/deltay,alpha);
		} else {
			Ges = A*pow(deltay,-1.5);
		}
		Ges /= (x*x);
		return w*Ges*planck/(electronMass*cLight2);
	} else return 0.0;
}