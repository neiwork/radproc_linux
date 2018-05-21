#include "coppiBlandford.h"

#include "modelParameters.h"
#include <fparameters/parameters.h>
#include <fmath/RungeKutta.h>
#include <fmath/physics.h>

#include <boost/property_tree/ptree.hpp>

double rate(double w, double g)
{
	double beta = sqrt(1.0-1.0/P2(g));

	double cte = 3.0*cLight*thomson/(32.0*P2(g*w)*beta);

	double min = 2.0*g*w*(1.0-beta);
	double max = 2.0*g*w*(1.0+beta);
	int nX = 10;

	double sum = 0.0;
	double x_int = pow( (max/min), (1.0/nX) );

	double x = min;
	for (int i = 0; i < nX; ++i)   
	{
		double dx = x*(x_int - 1);
		double f = (1.0-4.0/x -8.0/P2(x))*log(1.0+x) +0.5 + 8.0/x - 0.5/P2(1.0+x) ;

		sum = sum + f*dx;
		x = x*x_int;
	}


	double result = sum*cte;
	return result;
}


double w_ave(double wp, double g)  //<w>
{
	double beta = sqrt(1.0-1.0/P2(g));
	

	double min = -1.0;
	double max = 1.0;
	int nMu = 100;

	double sum = 0.0;
	double mu_int = (max-min)/nMu;

	double mu = min;
	for (int i = 0; i < nMu; ++i)   
	{
		//double dmu = mu*(mu_int - 1);
		
		double x = g*wp*(1.0-beta*mu);
		double a = P2(g)*beta*wp*(mu-beta);

		double factor1 = 3.0*cLight*thomson*x/(8.0*g*wp*rate(wp,g));

		double factor2 = 2.0*(a+g*(x-2.0)-(2.0*g+a)/x - a*(6.0/P2(x)+3.0/P3(x)))/(1.0+2.0*x);

		double factor3 = log(1.0+2.0*x)*(g - a + 3.0*a*(1.0/x+1.0/P2(x)))/P2(x);

		double factor4 = (1.0 - 1.0/P2(1.0+2.0*x))*( a*(3.0/P2(x)+1.0/P3(x)) +(a+g)/x +2.0*g )/(2.0*x);

		double factor5 = (1.0-1.0/P3(1.0+2.0*x))/3.0  * (g + a*(1.0/x+1.0/P2(x)))/3.0 -2.0*a/P3(x); //ojo que lo cambie

		double f = factor1*(factor2+factor3+factor4+factor5);

		sum = sum + f*mu_int/2.0;

		mu = mu + mu_int;
	}


	double result = sum;
	return result;
}

double w2_ave(double wp, double g)  //<w^2>
{
	double beta = sqrt(1.0-1.0/P2(g));

	double result = 0.0;

	if(wp*g < 1.0)
	{
		result = 14.0*P2(wp*P2(g))*(1.0-176.0*g*wp/35.0)/5.0; 
	}
	else
	{
		double a1 = 64.0*beta*P2(g*wp);
		double a2 = (6.0*g*wp*(1.0+beta)*(2.0*P2(g)+1.0)+6.0*P2(g)+3.0)*log(2.0*g*wp*(1.0+beta)+1);
		double a3 = (6.0*g*wp*(1.0-beta)*(2.0*P2(g)+1.0)+6.0*P2(g)+3.0)*log(2.0*g*wp*(1.0-beta)+1);
		double a4 = 9.0*P2(1.0-P2(g))*P3(g)*wp/32.0 - (58.0*P2(g)+1.0)/(64.0*g*wp) +7.0*(2.0*P2(g)*(1.0-P2(beta))+1)/32.0;

		result = (a2 - a3 + a4)/rate(wp,g)/a1;
	}

	return result;
}




double scatteredDist(double w, double wp, double g)  //P'(w,w',g)
{
	
	double w_max = pow(10.0,GlobalConfig.get<double>("model.particle.photon.dim.energy.max"))*1.6e-12/(electronMass*cLight2);
	double w_min = pow(10.0,GlobalConfig.get<double>("model.particle.photon.dim.energy.min"))*1.6e-12/(electronMass*cLight2);

	double wAve = w_ave(wp,g);
	double delta_w =  std::abs( w2_ave(wp,g) - P2(wAve));  //<w^2> -<w>^2

	double D = std::min( wAve-w_min , w_max - wAve );
	D = std::min(sqrt(3.0*delta_w),D);


	double heaviside = D-std::abs(w-wAve);

	double result = 0.0;
	if(heaviside >= 0)
	{
		result = 1.0/(2.0*D);
	}

	return result;

}


double f_uno(double g, double w, const Particle& creator, const SpaceCoord& distCoord)   //funcion a integrar  u=Ee
{    
	
	
	double Erest = creator.mass*cLight2;
	double Eval = g*Erest; 
	double distCreator;
	if (Eval < creator.emin() || Eval> creator.emax()){
		distCreator = 0.0;
	}
	else{
		distCreator = creator.distribution.interpolate({ { 0, Eval } }, &distCoord); 
	}
	  
	double function = rate(w, g)* (distCreator*Erest);


	return function;
}


double f_dos(double wp, double g, double E_norm, const Particle& creator, const SpaceCoord& distCoord, const ParamSpaceValues& tpf)   //funcion a integrar  u=Ee
{    
	
	
	double Erest = creator.mass*cLight2;
	double Eval = g*Erest;
	double distCreator;
	if (Eval < creator.emin() || Eval> creator.emax()){
		distCreator = 0.0;
	}
	else{
		distCreator = creator.distribution.interpolate({ { 0, Eval } }, &distCoord); 
	}
	  

	double Nph = tpf.interpolate({ { 0, wp*Erest } }, &distCoord);  

	double w = E_norm;

	double function =  scatteredDist(w, wp, g)*rate(w,g)*(distCreator*Erest)*(Nph*Erest);


	return function;
}


double thComptonCoppi(double E, const Particle& creator, const SpaceCoord& distCoord, const ParamSpaceValues& tpf, double phEmin)
{
	
	double Erest = creator.mass*cLight2;
	double E_norm = E/Erest;
	
	double gMin = creator.emin()/Erest;
	double gMax = creator.emax()/Erest;
	double eMin = phEmin/Erest;
	double eMax = E/Erest;



	double Nph = tpf.interpolate({ { 0, E } }, &distCoord);  
	

	double integral1 = RungeKuttaSimple(gMin,gMax, 
		[E_norm,&creator,&distCoord, tpf](double g){
			return f_uno(g, E_norm, creator, distCoord);  
		});     


	double integral2 = RungeKutta(eMin,eMax, //emax deberia ser std::min(eMax,E_norm)?
		[gMin](double wp){
			return gMin;
		}, 
		[gMax](double wp){
			return gMax;
		}, 
		[E_norm,&creator,&distCoord, tpf](double wp, double g){
			return f_dos(wp,  g,  E_norm,creator, distCoord, tpf); 
		});       

	double luminosity = integral1*Nph + integral2;

	return luminosity;

	}
	
	