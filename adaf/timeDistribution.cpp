#include "timeDistribution.h"



#include "losses.h"
#include "injection.h"

//#include <fmath\interpolation.h>
//#include <fmath\elimiGaussiana.h>
//#include <fmath\matrixInit.h>
#include <fmath\physics.h>

#include <fparameters\Dimension.h>

#include <iostream>
#include <fstream>






double timeDistribution(double Ee, double r, double t, Particle& p, State& st, double Eeff, const SpaceCoord& i)
{ 
	//const double magneticField{ st.magf.get(i) };
	double Emax = p.emax();// eEmax(r, magneticField);  //ver si genero la funcion XXX
	double tmin = st.electron.ps[2][0]; //en teoria esto no es necesario, porque solo llamo a etsa funcion en el tmin
	//chequear que t == tmin

	double EpMin = Ee;   //Ep = Eprima
	double EpMax = Eeff;
	double Ep = EpMin;

	int nEp = 50;
	double Ep_int = pow((EpMax / EpMin), 1.0 / nEp);  //0.001

	double dEp(0.0), inj(0.0);

	double sum_Ep = 0.0;

	for (int l=0; l < nEp; ++l)	{

		dEp = Ep*(Ep_int-1.0);

		if (Ep > Emax){
			inj = 0.0;
		}
		else{
			inj = p.injection.interpolate({ { DIM_E, Ep }, { DIM_R, r }}, &i); //, { DIM_T, tmin} }); //paso los valores E,r,t en donde quiero evaluar Q
			//interpolDoble(Ep, Time[j]-tau, Ee, Time, injection);  //chequear que ande!
		}

		sum_Ep = sum_Ep + inj*dEp;

		Ep = Ep*Ep_int;

	}

	double perdidas = losses(Ee, p, st, i);

	double total = sum_Ep / perdidas;
				
	return total;

}