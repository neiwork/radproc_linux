#include <iostream>

#include "redshiftFunction.h"
#include "adafFunctions.h"

extern double schwRadius;
extern size_t nR, nRcd;

Vector redshift_to_inf;							// Vector [jR] with the redshift between a shell and infinity.
Vector redshift_CD_to_inf;						// Vector [jRcd] with the redshift between the CD and infinity.
Matrix redshift;								// Matrix [jR][j'R] with the redshift between shells.
Matrix redshift_CD_to_RIAF;						// Matrix [jRcd][jR] with the redshift between the CD and a shell
												// of the RIAF.
Matrix redshift_RIAF_to_CD;						// Matrix [jR][jRcd] with the redshift between a shell of the
												// RIAF and the CD.

void redshiftFactor(State& st)
{
	redshift_to_inf.resize(nR,1.0);
	redshift_CD_to_inf.resize(nRcd,1.0);
	matrixInit(redshift,nR,nR,1.0);
	matrixInit(redshift_CD_to_RIAF,nRcd,nR,1.0);
	matrixInit(redshift_RIAF_to_CD,nR,nRcd,1.0);
	
	size_t jR=0;
	st.photon.ps.iterate([&](const SpaceIterator& itR) {
		double r = itR.val(DIM_R);
		double vr = radialVel(r);
		double beta = vr/cLight;
		if (abs(beta) >= 1.0) beta = -0.9;
		double redshift_factor = sqrt( (1.0-schwRadius/r) * (1.0-beta*beta) );
		redshift_factor = (redshift_factor > 0.0) ? redshift_factor : 1.0;
		redshift_to_inf[jR] = redshift_factor;
		//redshift_to_inf[jR] = 1.0;
		cout << "r [2M] = " << r/schwRadius << "\t (1+z)^-1 = " << redshift_to_inf[jR] << endl;
		size_t jjR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itRR) {
			double rr = itRR.val(DIM_R);
			double vrr = radialVel(rr);
			double relative_velocity = abs(vr-vrr) / (1.0 - vr*vrr/cLight2);
			double relative_beta = relative_velocity/cLight;
			if (abs(relative_beta) >= 1.0) relative_beta = 0.9;
			double doppler_factor = sqrt( (1.0-relative_beta) / (1.0+relative_beta) );
			double grav_factor_shells = sqrt( (1.0-schwRadius/r) / (1.0-schwRadius/rr) );
			redshift[jR][jjR] = grav_factor_shells * doppler_factor;
			//redshift[jR][jjR] = 1.0;
			jjR++;
		},{0,-1,0});
		
		size_t jRcd=0;
		st.photon.ps.iterate([&](const SpaceIterator& itRcd) {
			double rCD = itRcd.val(DIM_Rcd);
			if (rCD < r) beta = -beta;
			double doppler_factor = sqrt( (1.0-beta) / (1.0+beta) );
			double redshift_factor_CD = sqrt( (1.0-schwRadius/r) / (1.0-schwRadius/rCD) );
			redshift_RIAF_to_CD[jR][jRcd] = redshift_factor_CD * doppler_factor;
			jRcd++;
		},{0,0,-1});
		
		jR++;
	},{0,-1,0});
	
	size_t jRcd=0;
	st.photon.ps.iterate([&](const SpaceIterator& itRcd) {
		double rCD = itRcd.val(DIM_Rcd);
		double betaCD = rCD*keplAngVel(rCD)/cLight;
		double redshift_factor_CD_to_inf = sqrt( (1.0-schwRadius/rCD) * (1.0-betaCD*betaCD) );
		redshift_CD_to_inf[jRcd] = redshift_factor_CD_to_inf;
		
		size_t jjR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itRR) {
			double rr = itRR.val(DIM_R);
			double vrr = radialVel(rr);
			double beta = -vrr/cLight;
			if (abs(beta) >= 1.0) beta = 0.9;
			if (rCD < rr) beta = -beta;
			double doppler_factor = sqrt( (1.0-beta) / (1.0+beta) );
			double grav_factor = sqrt( (1.0-schwRadius/rCD) / (1.0-schwRadius/rr) );
			redshift_CD_to_RIAF[jRcd][jjR] = grav_factor * doppler_factor;
			jjR++;
		},{0,-1,itRcd.coord[DIM_Rcd]});
		
		jRcd++;
	},{0,0,-1});
	
}