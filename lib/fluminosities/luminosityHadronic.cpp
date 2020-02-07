#include "luminosityHadronic.h"

#include <fparameters/parameters.h>
#include <fmath/RungeKutta.h>
#include <fmath/interpolation.h>
#include <flosses/crossSectionInel.h>
#include <fmath/physics.h>
#include <algorithm>

#include <boost/math/special_functions/bessel.hpp>

/*double fpp(double Ep, double E)         //funcion a integrar   x=Eproton; E=Epion
{
	double L = log(Ep/1.6); //el 1.6 son TeV en erg
	double x = E/Ep; 
	double Bg = 1.3+0.14*L+0.011*L*L;
	double beta = 1.0/(1.79+0.11*L+0.008*L*L);
	double kappa = 1.0/(0.801+0.049*L+0.014*L*L);
	
	double equis_b = pow(x,beta);
	double factor = 1-equis_b;
	double f;
	if (factor =! 0)	{
		
		double f1 = Bg*log(x)/x;
		double f2 = factor/(1.0+kappa*equis_b*factor);
		double f3 = 1.0/log(x) - 4.0*beta*equis_b/factor - 4.0*kappa*beta*equis_b*(1.0-2.0*equis_b)/(1.0+kappa*equis_b*factor);
		
		f      =  f1*pow(f2,4)*f3;
	}
	else	{
		f = 0.0;
	}
	return f;		
}*/

double cHadron(double E)  //limite inferior
{	
	double Kpi = 0.17;
	double thr = 0.0016; //1GeV

	return std::max(E+P2(neutralPionMass*cLight2)/(4*E),thr*Kpi); //== Ekin > Ethr
}

double dHadron()         //limite superior   
{
	return 1.6e-12*pow(10.0,13.0);   //esto es un infinito                                
}


/*double fHadron(double Ep, double E, const Particle& p, const SpaceCoord& psc) //funcion a integrar   x=Ecreator; L=L(Ega)
{	
	double Kpi = 0.17;
	double eval = p.mass*cLight2+Ep/Kpi;
	double Ekin = Ep/Kpi;
	double distCreator=0.0;
	if (Ep < p.emax() && Ep >= p.emin()) {
		distCreator = p.distribution.interpolate({ { 0, Ep } }, &psc);
	}
	
	double thr = 0.0016; //1GeV
	//double sigma = 30e-27*(0.95+0.06*log(Ekin/thr));
	
	double l = log10((protonMass*cLight2+Ep/Kpi)/1.6);
	double sigma = 1.e-27 * (34.3+1.88*l+0.25*l*l);
	double pionEmiss = sigma*distCreator*fpp(Ep,E)/Ep; //Kpi;  //sigma = crossSectionHadronicDelta(Ekin)
	
	double result = pionEmiss; ///sqrt(P2(Epi)-P2(neutralPionMass*cLight2));
	return result;
}*/

double heaviside(double x,double a,double b)
{
	return (a <= x && x <= b ? 1.0 : 0.0);
}

double auxf3(double dGeV, double sGeV, double gx, double GammaGeV, double isoGeV)
{
	double pGeV = protonMass*cLight2/1.6e-3;
	double piGeV = neutralPionMass*cLight2/1.6e-3;
	double gd = (sGeV + dGeV*dGeV - pGeV*pGeV)/(2.0*dGeV*sqrt(sGeV));
	double betad = sqrt(1.0-1.0/(gd*gd));
	double epi = (dGeV*dGeV+piGeV*piGeV-pGeV*pGeV)/(2.0*dGeV);
	double ppi = sqrt(epi*epi-piGeV*piGeV);
	double aux1 = 0.5/(betad*gd*ppi);   // REVISAR ESTE p_pi
	double aux2 = P2(dGeV-isoGeV)+GammaGeV*GammaGeV;
	double h = heaviside(gx*piGeV,gd*(epi-betad*ppi),gd*(epi+betad*ppi));
	
	return aux1*h/aux2;    // [GeV^-3]
	
}

double inclusiveSigma(double sGeV)
{
	double pGeV = protonMass*cLight2/1.6e-3;
	double piGeV = neutralPionMass*cLight2/1.6e-3;
	double eta=sqrt(P2(sGeV-P2(piGeV)-P2(2.0*pGeV))-4.0*P2(piGeV*2.0*pGeV))/(2.0*piGeV*sqrt(sGeV));
	double pthr = 0.78; //[GeV]
	double sigma=0.0;;
	
	// s = 2mp(Ek+2mp) = 2mp(E+mp)
	// E = s/2mp - mp
	// p = sqrt(E**2-m**2)
	double p = sqrt(P2(0.5*sGeV/pGeV-pGeV)-P2(pGeV));
	
	if (p >= pthr && p <= 0.96) {
		sigma = 0.032*eta*eta+0.04*pow(eta,6)+0.047*pow(eta,8);
	} else if (p > 0.96 && p <= 1.27) {
		sigma = 32.6*pow(p-0.8,3.21);
	} else if (p > 1.27 && p <= 8.0) {
		sigma = 5.4*pow(p-0.8,0.81);
	} else if (p > 8.0) {
		sigma = 32.0 * log(p) + 48.5 / sqrt(p) - 59.5;
	}
	return 1.0e-27 * sigma;     // [cm^2]
}
/*
double dsigma(double gx, double gr, double sGeV)
{
	double GammaGeV = 0.0575;
	double pGeV = protonMass*cLight2/1.6e-3;
	double piGeV = neutralPionMass * cLight2 / 1.6e-3;
	double isoGeV = 1.236;
	double atan1 = atan((sqrt(sGeV)-pGeV-isoGeV)/GammaGeV);
	double atan2 = atan((pGeV+piGeV-isoGeV)/GammaGeV);
	double aux1 = GammaGeV/(atan1-atan2);
	
	double Min = pGeV+piGeV;
	double Max = sqrt(sGeV)-pGeV;
	double integ= integSimpson(log(Min),log(Max),[&](double log_dGeV)
					{
						double dGeV = exp(log_dGeV);
						return auxf3(dGeV,sGeV,gx,GammaGeV,isoGeV)*dGeV;
					},30);
	return inclusiveSigma(sGeV)*aux1*integ*piGeV;  // [cm^2]
}*/

double dsigma(double gx, double gr, double sGeV)
{
	double GammaGeV = 0.0575;
	double pGeV = protonMass*cLight2/1.6e-3;
	double piGeV = neutralPionMass * cLight2 / 1.6e-3;
	double isoGeV = 1.236;
	double atan1 = atan((sqrt(sGeV)-pGeV-isoGeV)/GammaGeV);
	double atan2 = atan((pGeV+piGeV-isoGeV)/GammaGeV);
	double aux1 = GammaGeV/(atan1-atan2);
	
	double Min = pGeV+piGeV;
	double Max = sqrt(sGeV)-pGeV;
	
	double integ= RungeKuttaSimple(Min,Max,[&](double dGeV)
					{return auxf3(dGeV,sGeV,gx,GammaGeV,isoGeV);});
	
	return inclusiveSigma(sGeV)*aux1*integ*piGeV;  // [cm^2]
}

double auxf2(double gx, double gr, double epi, double normtemp, double sGeV)
{
	double piGeV = neutralPionMass*cLight2/1.6e-3;
	double g = epi/piGeV;
	double beta = sqrt(1.0-1.0/(g*g));
	double betax = sqrt(1.0-1.0/(gx*gx));
	double q = sqrt(2.0*(gr+1.0))/normtemp;
	
	double f1 = exp(-q*g*gx*(1.0-beta*betax))-exp(-q*g*gx*(1.0+beta*betax));
	
	return f1/(betax*gx) * dsigma(gx,gr,sGeV);
}

/*
double auxf(double gr, double epi, double normtemp)
{
	double piGeV = neutralPionMass*cLight2 / 1.6e-3;
	double pGeV = protonMass*cLight2/1.6e-3;
	double sGeV = 2.0*P2(pGeV)*(gr+1.0);
	double ji = (sGeV-4.0*P2(pGeV)+P2(piGeV))/(2.0*sqrt(sGeV));
	double Max = ji/piGeV;
	double integral = integSimpson(0.0,log(Max),[gr,epi,normtemp,sGeV](double log_gx)
						{
							double gx = exp(log_gx);
							return auxf2(gx,gr,epi,normtemp,sGeV)*gx;
						},30);
	double result = (gr*gr-1.0) / sqrt(2.0*(gr+1.0)) * integral;
}*/

double auxf(double gr, double epi, double normtemp)
{
	double piGeV = neutralPionMass*cLight2 / 1.6e-3;
	double pGeV = protonMass*cLight2/1.6e-3;
	double sGeV = 2.0*P2(pGeV)*(gr+1.0);
	double ji = (sGeV-4.0*P2(pGeV)+P2(piGeV))/(2.0*sqrt(sGeV));
	double Max = ji/piGeV;
	double integral = RungeKuttaSimple(1.0,Max,[gr,epi,normtemp,sGeV](double gx)
						{return auxf2(gx,gr,epi,normtemp,sGeV);});
	//double integral = qMidPointLog(1.0,Max,[gr,epi,normtemp,sGeV](double gx)
	//					{return auxf2(gx,gr,epi,normtemp,sGeV);},1.0e-2);
	double result = (gr*gr-1.0) / sqrt(2.0*(gr+1.0)) * integral;
}

/*
double fPion(double epi, double density, double temp)
{
	double normtemp = boltzmann* temp / (protonMass*cLight2);
	double k2 = boost::math::cyl_bessel_k(2, 1.0/normtemp);
	double piGeV = neutralPionMass*cLight2/1.6e-3;
	double constant = cLight*density*density/(4.0*piGeV*normtemp*k2*k2);
	return constant * integSimpson(0.0,log(1.0e3),[epi,normtemp](double log_gr)
						{
							double gr = exp(log_gr);
							return auxf(gr,epi,normtemp)*gr;
						},30);	//  [cm-3 s-1 GeV-1]
}*/

double fPion(double epi, double density, double temp)
{
	double normtemp = boltzmann* temp / (protonMass*cLight2);
	double k2 = boost::math::cyl_bessel_k(2, 1.0/normtemp);
	double piGeV = neutralPionMass*cLight2/1.6e-3;
	double constant = cLight*density*density/(4.0*piGeV*normtemp*k2*k2);
	
	//return constant * RungeKuttaSimple(1.0,1.0e3,[&](double gr)
	//					{return auxf(gr,epi,normtemp);});         //  [cm-3 s-1 GeV-1]
	return constant * qImpropLog(1.0,1.0e3,[epi,normtemp](double gr)
						{return auxf(gr,epi,normtemp);},1.0e-2);
}


double fHadron(double epi, double density, double temp)
{
	double piGeV = neutralPionMass*cLight2/1.6e-3;
	double qpi = fPion(epi,density,temp);
	return qpi / sqrt(P2(epi)-P2(piGeV));   // [cm-3 s-1 GeV-2]
}

/*
double luminosityHadronic(double E,
	const double density, double temp)
{
	double Kpi = 0.17;
	double thr = 0.0016; //1GeV

	//double Max  = dHadron();   //esto es un infinito 
	double Min  = cHadron(E);
	//double Min = E+P2(0.5*neutralPionMass*cLight2)/E;
	Min = Min / 1.6e-3;
	double Max = 1.0e3;
	double integral = 2.0*integSimpson(log(Min),log(Max),[density,temp](double logepi)
						{
							double epi = exp(logepi);
							return fHadron(epi,density,temp)*epi;
						},30); // [cm-3 s-1 GeV-1]
		
	double jpp = integral * E*planck * 0.25/pi; // [erg s^-1 Hz^-1 cm^-3]
	jpp = jpp / (1.6e-3);
	return jpp;
}
*/

double luminosityHadronic(double E,
	const double density, double temp)
{
	double Kpi = 0.17;
	double thr = 0.0016; //1GeV

	//double Max  = dHadron();   //esto es un infinito 
	double Min  = cHadron(E);
	//double Min = E+P2(0.5*neutralPionMass*cLight2)/E;
	Min = Min / 1.6e-3;
	double Max = 1.0e3;  // [en GeV]
	
	//double integral = 2.0*cLight*density*RungeKuttaSimple(Min, Max, 
	//	[&E,&p,&psc](double x) {return fHadron(x, E, p, psc); });    //integra entre Emin y Emax
	
	double integral = 2.0*RungeKuttaSimple(Min,Max,[&](double epi)
							{return fHadron(epi,density,temp);}); // [cm-3 s-1 GeV-1]
	//double integral = 2.0*qMidPointLog(Min,Max,[density,temp](double epi)
	//						{return fHadron(epi,density,temp);},1.0e-2);
	//double integral = 2.0*qromoLog(Min,Max,[density,temp](double epi){return fHadron(epi,density,temp);});
	/*double integral = 2.0*cLight*density*RungeKutta(p.emin(), p.emax(),
		[E](double u){
			return cHadron(u,E);
		}, 
		[E](double u){
			return dHadron(u,E);
		}, 
		[&p,&psc](double u, double t){
			return fHadron(u,t,p, psc);  
		}); */
		
	double jpp = integral * E*planck * 0.25/pi; // [erg s^-1 Hz^-1 cm^-3]
	jpp = jpp / (1.6e-3);
	return jpp;
}



/*
double fPP(double x, double E, Particle& creator )         //funcion a integrar   x=Eproton; E=Epion
{
	//DataInjection* data = (DataInjection*)voiddata;
	//double E = data->E;    //esta E corresponde a la energia del foton emitido; E=Ega
	//double mass   = data->mass;
	//Vector& Ncreator = data->Ncreator;
	//Vector& Ecreator = data->Ecreator;
	//const double mass   = particle.mass;
	//const Vector& Ncreator = creator.distribution.values;
	//const Vector& Ecreator = creator.eDim()->values;
	double L      = log(x/1.6); //el 1.6 son TeV en erg
	double ap     = 3.67+0.83*L+0.075*P2(L);
	double se     = crossSectionHadronic(x);
	double Bp     = ap+0.25;
	double r      = 2.6*pow(ap,-0.5);
	double alpha  = 0.98*pow(ap,-0.5);
	double equis  = E/x;
	double factor = 1-pow(equis,alpha);
	double f;
	if (factor =! 0)	{
		f      = 4*alpha*Bp*pow(equis,(alpha-1))*pow((factor/(1+r*pow(equis,alpha)*factor)),4)
      	         *(1/factor+r*(1-2*pow(equis,alpha))/(1+r*pow(equis,alpha)*factor))*
     	         pow((1-chargedPionMass*cLight2/(equis*x)),0.5);
	}
	else	{
		f = 0.0;
	}*/