#include "probabilityMatrix.h"
#include "modelParameters.h"
#include <fparameters/Dimension.h>
#include <fparameters/SpaceIterator.h>
#include <fparameters/parameters.h>
#include <iostream>
#include <boost/property_tree/ptree.hpp>

double optDep2(double r1_sup,double tprim,double pprim,double y0,double z0, ParamSpaceValues denf_e, SpaceIterator i)
{
    double suma = 0.0;
    double r1_inf = 0.0;  //limites para la integral sobre r1
    int nr1 = 10;
    double r1 = r1_inf;

    for (int h=0;h<nr1+1;h++) { //aca comienza la integral
        double dr1 = r1_sup/nr1;
        //coord de P desde O'
        double x1 = r1*sin(tprim)*cos(pprim);
        double y1 = r1*sin(tprim)*sin(pprim);
        double z1 = r1*cos(tprim);
        //coord de P desde O, que son la grilla de los psv
        double x2 = x1; //+x0=0
        double y2 = y1 + y0;
        double z2 = z1 + z0;
        double r2 = sqrt(x2*x2+y2*y2+z2*z2);
        double theta2=0.0;
        double den = 0.0;
        try {
            den = denf_e.interpolate({ {DIM_R, r2} , {DIM_THETA, theta2} }, &(const SpaceCoord&) i);
        } catch (std::runtime_error& e) {
            std::cout << "WARNING: " << e.what() << std::endl;
        }
        suma+=thomson*den*dr1;
        r1+=dr1;
    }
    return suma;
}
                            
void probabilityMatrix2(State& st, Matrix& a)
{
	double rg=GlobalConfig.get<double>("rg");
	Particle& p=st.electron; 
	int nR = p.ps[DIM_R].size();  //ver el -1
	matrixInit(a, nR, nR, 0.0);

    /*size_t z_i=0;
    p.ps.iterate([&](const SpaceIterator& i) {
		for (size_t z_j=1;z_j<nR;z_j++) { 
				const double r_i=p.ps[DIM_R][z_i];
				const double r_j=p.ps[DIM_R][z_j];
                const double r_jant=p.ps[DIM_R][z_j-1];
				double dr=std::abs(r_i-r_j);
                double dr2=std::abs(r_i-r_jant);
                const double theta = 0.0;
                double y0 = r_i*cos(theta);
                double z0 = r_i*sin(theta);
                double tau_1=rg*optDep2(dr,0.0,pi/2,y0,z0,st.denf_e,i);
                double tau_2=rg*optDep2(dr2,0.0,pi/2,y0,z0,st.denf_e,i);
                double P=exp(-tau_2)-exp(-tau_1);
                P=tau_1-tau_2;
                a[z_i][z_j] = P;
        }
        z_i+=1;
    }, {0,-1,0});*/
	
	
	//size_t z_i=0;
	//const double theta = 0.0;
	
    p.ps.iterate([&](const SpaceIterator& i) {
		
		const double r_i = i.val(DIM_R);
		const double theta = i.val(DIM_THETA);
		
		int z_i = i.coord[DIM_R];

		for (size_t z_j=z_i+1 ; z_j < nR; z_j++) { 
			
				const double r_j = p.ps[DIM_R][z_j];
				
				double dr=std::abs(r_i-r_j);
				
				const double r_jant = p.ps[DIM_R][z_j-1];
				
				double dr2=std::abs(r_i-r_jant);
				
				double y0 = r_i*cos(theta);
				double z0 = r_i*sin(theta);
				double tau_1=rg*optDep2(dr,0.0,pi/2,y0,z0,st.denf_e,i);
				double tau_2=rg*optDep2(dr2,0.0,pi/2,y0,z0,st.denf_e,i);
				double P=exp(-tau_2)-exp(-tau_1);
				P=tau_1-tau_2;
				a[z_i][z_j] = P;
				
        }
        //z_i+=1;
    }, {0,-1,0});
}
