#include "probabilityMatrix.h"

#include "modelParameters.h"
#include <fparameters/Dimension.h>
#include <fparameters/SpaceIterator.h>

#include <fparameters/parameters.h>

#include <iostream>
#include <boost/property_tree/ptree.hpp>


double optDep(double r1_sup,double tprim,double pprim,double y0,double z0, ParamSpaceValues denf_e, SpaceIterator i)
{
    
    double suma = 0.0;
    
    double r1_inf = 1.0e-3;  //limites para la integral sobre r1
    int nr1 = 10;
    double r1_int = pow((r1_sup / r1_inf), (1.0 / nr1));
                
    double r1 = r1_inf;

    for (int h = 0; h < nr1+1; h++) { //aca comienza la integral

        double dr1 = r1*(r1_int-1.0);
                        
        //coord de P desde O'
        double x1 = r1*sin(tprim)*cos(pprim);
        double y1 = r1*sin(tprim)*sin(pprim);
        double z1 = r1*cos(tprim);

        //coord de P desde O, que son la grilla de los psv
        double x2 = x1; //+x0=0
        double y2 = y1 + y0;
        double z2 = z1 + z0;

        double r2 = sqrt(x2*x2+y2*y2+z2*z2);
        //double phi2;
        //if(x2 == 0.0) {
        //	phi2 = pi/2.0;}
        //	else{ double phi2 = atan(y2/x2);}
            
        double theta2 = std::abs(asin(z2/r2));  //arcsen


        double den = 0.0;
        try {
            den = denf_e.interpolate({ {DIM_R, r2} , {DIM_THETA, theta2} }, &(const SpaceCoord&) i);
        } catch (std::runtime_error& e) {
            std::cout << "WARNING: " << e.what() << std::endl;
            //return 0.0;
        }		

        suma += 
            thomson*den*dr1;
                                        


        r1=r1*r1_int;

    }
    return suma;
}
                 
                            
void probabilityMatrix(State& st, Matrix& a)
{
    
	double rg = GlobalConfig.get<double>("rg");
	Particle& p = st.electron; 

	Vector theta_prim(3, pi/4.0);
	Vector phi_prim(6, pi/4.0);
	
	for(int t_i = 0; t_i < theta_prim.size(); t_i++){
		theta_prim[t_i] = theta_prim[t_i]*(t_i+1);
	}
    
    for(int t_i = 0; t_i < phi_prim.size(); t_i++){
        if (t_i < phi_prim.size()/2){
            phi_prim[t_i] = phi_prim[t_i]*(t_i+1);
        }else{
            phi_prim[t_i] = phi_prim[t_i]*(t_i+2);
        }
	}


	int nR = p.ps[DIM_R].size();  //ver el -1

	//Matrix a;
	matrixInit(a, nR, nR, 0.0);  //inicializo la matriz con 0, y le doy dimensiones nRxnR


	for (size_t z_i = 0; z_i < nR; z_i++) { 
		
		
		for (size_t z_j = 0; z_j < nR; z_j++) { 

			double aveP = 0.0;
			double Pij = 0.0;

			if(z_i != z_j){

				double weight, P = 0.0;

				const double r_i = p.ps[DIM_R][z_i];
				const double r_j = p.ps[DIM_R][z_j];
				double dr = std::abs(r_i-r_j);  


				p.ps.iterate([&](const SpaceIterator& i) { //iterate sobre theta nada mas

					const double theta = i.val(DIM_THETA);
					//nt t_ix = i.coord[DIM_theta]; //posicion en la coordenada theta
					//if (t_ix % 2 == 0) {}  //quizas pueda usar esta condicion para hacerlo en la mitad de los thetas

					//el sistema primado ya esta definido, con origen en O'=(r,theta,phi=pi/2) o (x0,y0,z0)
					//double x0 = 0.0; //por phi = pi/2
					double y0 = r_i*cos(theta); //cambiar porque el theta del toro se mide desde el ecuador
					double z0 = r_i*sin(theta);

					weight = sin(theta)/(4.0*pi);  //el peso para la probabilidad es el angulo solido, ver cual es el theta q debo usar
                    
                    //desde aca se tiran 3x3(para cada lado) fotones en las direcciones de theta_prim,phi_prim
                    for (size_t t_i = 0; t_i < theta_prim.size(); t_i++) { 
                        
                        
                        double P = 0.0;
                                                
                        if(z_j > z_i){//hacia r crecientes
                        
                            for (size_t p_i = 0; p_i < phi_prim.size()/2; p_i++) { 
                                
                                double r1_sup = dr/(sin(theta_prim[t_i])*sin(phi_prim[p_i]));
                                double opticalDepth = rg*optDep(r1_sup,theta_prim[t_i],phi_prim[p_i],y0,z0,st.denf_e,i);
                                
                                P += 1.0 - exp(-opticalDepth);
                            }
                        } 
                        else{//hacia r menores
                            
                            for (size_t p_i = phi_prim.size()/2; p_i < phi_prim.size(); p_i++) { 
                                
                                double r1_sup = std::abs(dr/(sin(theta_prim[t_i])*sin(phi_prim[p_i])));
                                
                                double opticalDepth = rg*optDep(r1_sup,theta_prim[t_i],phi_prim[p_i],y0,z0,st.denf_e,i);
                                
                                P += 1.0 - exp(-opticalDepth);
							
                            }							
							
						}						                     
                                                   
						aveP += P;//*weight; 
					}


				}, { 0, z_i, -1 });  //fijo cualquier energia; recorro solo los theta
			}//z_i != z_j
	
			Pij += aveP/(theta_prim.size()+phi_prim.size()+p.ps[DIM_THETA].size());

			a[z_i][z_j] = Pij;
		}	
	}
	
}


//	Vector Ee(nE, 0.0);
