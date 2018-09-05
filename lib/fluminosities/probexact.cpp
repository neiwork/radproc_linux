#include "probexact.h"
#include <math.h>
#include <fmath/physics.h>

double rate(double om, double gamma)
{
	double beta=sqrt(1.0-1.0/(gamma*gamma));
    
    if (beta < 1.0e-2 && om < 1.0e-2) {
        return cLight*thomson*(1.0-2.0*gamma*om);
    } else {
        double cte=3.0*cLight*thomson/(32.0*(gamma*om)*(gamma*om)*beta);
        double min=2.0*gamma*om*(1.0-beta);
        double max=2.0*gamma*om*(1.0+beta);
        int nX=100;
        double sum=0.0;
        double x_int=pow(max/min,1.0/nX);
        double x = min;
        for (int i=0;i<nX;++i)   
        {
            double dx = x*(x_int-1.0);
            double f=(1.0-4.0/x-8.0/(x*x))*log(1.0+x)+0.5+8.0/x-0.5/((1.0+x)*(1.0+x));
            sum += f*dx;
            x = x*x_int;
        }
        double result = cte*sum;
        return result;
    }
}

double dPdz(double om,double omprim,double gamma,double z) {
    
    double beta=sqrt(1.0-1.0/(gamma*gamma));
    double rho=om/omprim;
    rho=10.0;
    double epsi=omprim/gamma;
    double ka=gamma/om;
    double aux=rho*(beta*beta+epsi*epsi+2.0*beta*epsi*z);
    double y0=(epsi+beta*z)*(rho+epsi*rho-1.0+beta*z)/aux;
    double aux1=rho*rho*beta*beta;
    double aux2=2.0*rho*epsi*(1.0-rho)*(1.0-beta*z);
    double aux3=P2(rho-1.0+beta*z);
    
    double delta=beta*sqrt(1.0-z*z)*sqrt(aux1+aux2-aux3)/aux;
    double a=1.0-beta*z-(1.0-y0)/ka;
    double b=delta/ka;
    
    double num = 3.0*thomson*cLight*om;
    double r=rate(omprim,gamma);
    double den = pow(2.0*gamma,4)*omprim*omprim*r*sqrt(beta*beta+epsi*epsi+
                                2.0*beta*epsi*z)*(1.0-beta*z);
    double factor = 2.0*y0*ka-a*ka*ka+pow(a*a-b*b,-0.5) * (1.0+y0*y0-2.0*a*y0*ka+a*a*ka*ka) +
                                1.0/(ka*ka*(1.0-beta*z)) * (ka*ka*(1.0 + pow(a*a-b*b,-0.5)*a*(2.0*b-a)/(a-b))+
                                pow(a*a-b*b,-1.5)*(a*(1.0-y0)*(1.0-y0)+2.0*ka*b*b*(1.0-y0)-b*a*a*ka*ka));
    
    double result=num * factor / den;
    if (result > 0.0) return result;
    else return 0.0;
}

double probexact(double om,double omprim,double gamma) 
{
    double rho = om/omprim;
    double beta = sqrt(1.0-1.0/(gamma*gamma));
    double epsi= omprim/gamma;
    double d=1.0+epsi-epsi*rho;
    double aux1=sqrt(d*d-1.0/(gamma*gamma));
    double auxmax= (1.0 - rho*(d-aux1))/beta;
    double auxmin=(1.0-rho*(d+aux1))/beta;
    double zmax= (0.99 > auxmax ? auxmax : 1.0);
    double zmin= (-0.99 < auxmin ? auxmin : -1.0);
    int nZ=100;
    double dz=(zmax-zmin)/nZ;
    double sum=0.0;
    double z=zmin;
    for (int i=0;i<nZ;i++) {
        double dprob=dPdz(om,omprim,gamma,z);
        sum += dprob * dz;
        z += dz;
    }
    return (sum > 0.0 ? sum : 0.1);
}