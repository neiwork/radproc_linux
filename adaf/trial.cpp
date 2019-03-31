#include <fluminosities/probexact.h>
#include <gsl/gsl_math.h>
#include <nrMath/brent.h>
#include <fmath/physics.h>
#include <nrMath/integrators.h>
#include <iostream>
using namespace std;

void trial()
{
	gsl_function gsl_extrinf;
	gsl_function gsl_probexact;
	double g = 10.0;
	double beta = sqrt(1.0-1.0/(g*g));
	double zero = 1.0e-4;

	double domp = pow(10.0/0.01,1.0/10);
	double omp = 0.01;
	for (int i=0;i<10;i++) {
    
		double eps = omp/g;
		
		struct two_d_params extr_params = {eps,g};
		struct two_d_params probexact_params = {omp,g};
		gsl_function gsl_extrinf;
			gsl_extrinf.function = &extrinf;
			gsl_extrinf.params = &extr_params;
		gsl_function gsl_extrsup;
			gsl_extrsup.function = &extrsup;
			gsl_extrsup.params = &extr_params;
		gsl_function gsl_probexact;
			gsl_probexact.function = &gslprobexact;
			gsl_probexact.params = &probexact_params;
		
		int status,status1,status2;
		size_t numeval;
		double error;
		
		double omMin = omp*brent(&gsl_extrinf,0.0,1.0,&status1,&status2);
		double omMaxAbs = omp + (g - 1.0);
		double omMax = omp*(1.0+beta)/(1.0-beta+2.0*omp/g);
		double aux = beta/(1.0+g*(1.0+beta));
		omMax = (eps < aux) ? omMax : omMaxAbs;
		cout << omMin << "\t" << omMax << endl;
		double result = integrator_cquad(&gsl_probexact,omMin,omMax,1.0e-4,1.0e-4,100,&error,&status);
		cout << omp << "\t" << g << "\t" << "\t"<< status << "\t" << result << endl;
		omp *= domp;
	}
}