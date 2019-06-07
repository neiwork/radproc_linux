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

void targetField(State& st)
{

	
	st.photon.ps.iterate([&](const SpaceIterator& itER) {
			
		double E = itER.val(DIM_E);
		double frecuency=E/planck;
	
	
		double r=itER.val(DIM_R);
		double thetaH = st.thetaH.get(itER);
		double rB1=r/sqrt(paso_r);
		double rB2=r*sqrt(paso_r);
		double area=2.0*pi*rB2*rB2*(2.0*cos(thetaH)+sin(thetaH)*sin(thetaH));
		double fluxToLum=area;
		double vol=(rB2*rB2*rB2-rB1*rB1*rB1)*(4.0/3.0)*pi*cos(thetaH);
		double emissToLum=vol*4.0*pi;


		double temp=st.tempElectrons.get(itER);
		double temp_i=st.tempIons.get(itER);
		double magf=st.magf.get(itER);
		double dens_i=st.denf_i.get(itER);
		double dens_e=st.denf_e.get(itER);


		double jSy = jSync(E,temp,magf,dens_e);
		double lumRJ = bb_RJ(frecuency,temp) * fluxToLum;
		double lumSy = jSy*emissToLum;

		double lumBr = jBremss(E,temp,dens_i,dens_e)*emissToLum;
		
		double lumPp = 0.0;
		if (temp_i > 1.0e11 && frecuency > 1.0e20 && frecuency < 1.0e26) {
			double jpp = luminosityHadronic(E,dens_i,temp_i);
			lumPp = jpp*emissToLum;
		}
		
		double lumOut = lumRJ+lumSy+lumBr+lumPp;
		
		st.photon.distribution.set(itER,lumOut/(fluxToLum*cLight*planck*E));
			
	},{-1,-1,0}); //{itE.coord[DIM_E],-1,0});
//		jE++;	
	//},{-1,0,0});
}