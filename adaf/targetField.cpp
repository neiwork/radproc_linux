#include "targetField.h"
//#include "read.h"
#include "globalVariables.h"
#include "thermalLuminosities.h"
#include "adafFunctions.h"

#include <fluminosities/thermalSync.h>
#include <fluminosities/thermalBremss.h>
#include <fluminosities/luminosityHadronic.h>
#include <fluminosities/blackBody.h>

#include <fparameters/SpaceIterator.h>
//#include <fparameters/Dimension.h>
//#include <fparameters/parameters.h>

void targetField(State& st, Matrix lumOut, Vector energies)
{
	size_t jE=0;
	st.photon.ps.iterate([&](const SpaceIterator& itE) {
		size_t jR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itER) {
			double r=itER.val(DIM_R);
			double thetaH = st.thetaH.get(itER);
			double rB2=r*sqrt(paso_r);
			double area=2.0*pi*rB2*rB2*(2.0*cos(thetaH)+sin(thetaH)*sin(thetaH));
			st.photon.distribution.set(itER,lumOut[jE][jR]/
								(area*cLight*planck*energies[jE])); //erg^⁻1 cm^-3
			jR++;
		},{itE.coord[DIM_E],-1,0});
		jE++;
	},{-1,0,0});
}