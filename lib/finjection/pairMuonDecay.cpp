#include "pairMuonDecay.h"


#include <fparameters/parameters.h>
#include <fmath/RungeKutta.h>
//#include <fmath\interpolation.h>
#include <fmath/physics.h>
//#include <algorithm>

double cMuonDec(double Gep, double Ge)
{
	//double Ge = E/(electronmass*cLight2); 

	double Gmu_min  = Ge*Gep-sqrt((P2(Ge)-1.0)*(P2(Gep)-1.0));

	return Gmu_min;
}

double dMuonDec(double Gep, double Ge)     //Gep=gama E prima
{
	///double Ge = E/(electronmass*cLight2); 

	double Gmu_max = Ge*Gep+sqrt((P2(Ge)-1.0)*(P2(Gep)-1.0));

	return Gmu_max;               
}

double fMuonDec(double Gep, double Gmu, Particle& c, const SpaceCoord& distCoord)   //x=Gep; y=Gmu; c = muon
{ 	
	
	double Emu = Gmu*muonMass*cLight2;
	
	double Tdec = muonMeanLife*Emu/(muonMass*cLight2);
	
	double Nmu = 0.0;

	if (Emu > c.emin() && Emu < c.emax()){
		Nmu = c.distribution.interpolate({ { 0, Emu} }, &distCoord);
	}
	
	double Qmu = Nmu*muonMass*cLight2/Tdec;  //N(G) = N(E)*mc^2

	double Gep_max = 104;
		
	double P = 2.0*P2(Gep)*(3.0-2.0*Gep/Gep_max)/P3(Gep_max);
	double f = 0.5*P*Qmu/sqrt((P2(Gep)-1.0)*(P2(Gmu)-1.0));

	return f;

}



double pairMuonDecay(double E, Particle& c, const SpaceCoord& distCoord)
{	
	
	double inf = 2.0;
	double sup = 104.0;

	double Gmu_min = c.emin()/(muonMass*cLight2);
	double Gmu_max = c.emax()/(muonMass*cLight2);
	
	double Ge = E /(electronMass*cLight2); //aparece solo en los limites

		double integral  = RungeKutta(inf,sup, 
		[Ge,Gmu_min](double Gep) {return max(Gmu_min,cMuonDec(Gep, Ge)); },  //limite inferior
		[Ge,Gmu_max](double Gep) {return max(Gmu_max,dMuonDec(Gep, Ge)); },	 //limite superior
		[&c,&distCoord](double Gep, double Gmu) {return fMuonDec(Gep, Gmu, c, distCoord); });
		
	//double integral  = RungeKutta(a,b,bind(cMuonDec,_1,E),bind(dMuonDec,_1,E),bind(fMuonDec,_1,_2,E,creator));


	double inj = integral/(electronMass*cLight2);

	return inj;

}


	


	
	
