void comptonNew2(Matrix& lumOut, double lumIn, Vector energies, size_t nG, size_t nE,
				size_t nOm, double normtemp, Vector probab, size_t jEprim, size_t jR) 
{
	gsl_function gsl_extrinf;
	gsl_extrinf.function = &extrinf;
	
    double omPrim=energies[jEprim]/(electronMass*cLight2);
    double var_int=pow(energies[nE-1]/energies[0],1.0/nE);
    double dnuprim=energies[jEprim]/planck * (var_int-1.0);
    double alf=0.0;
    double *abscissas,*weights;
    abscissas=dvector(1,nG);
    weights=dvector(1,nG);
    double sum2=0.0;
    gaulag(abscissas,weights,nG,alf);
    Vector lumOut1(nE,0.0);
	
    for (size_t jG=1;jG<=nG;jG++) {
        double gamma = abscissas[jG]*normtemp+1.0;
        double beta = sqrt(1.0-1.0/(gamma*gamma));
        double cWeights = weights[jG]*gamma*sqrt(gamma*gamma-1.0);
		double eps = omPrim/gamma;
		
		struct two_d_params extrinf_params = {eps,gamma};
		gsl_extrinf.params = &extrinf_params;
		int status1,status2;
		
		double omMin = omPrim * brent(&gsl_extrinf,0.0,1.0,&status1,&status2);
		double omMaxAbs = omPrim + (gamma - 1.0);
		double omMax = omPrim*(1.0+beta)/(1.0-beta+2.0*omPrim/gamma);
		double aux = beta/(1.0+gamma*(1.0+beta));
		omMax = (eps < aux) ? omMax : omMaxAbs;
		
		double var_int_om=pow(omMax/omMin,1.0/nOm);
		double om=omMin;
		
		size_t jE=1;
		for (size_t jOm=0;jOm<nOm;jOm++) {
			double lumOutAux = lumIn * double(probab[((jR*nG+jG)*nE+jE)*nOm+jOm])*(om/omPrim);
			double lim1=sqrt(energies[1]*energies[0]);
			double lim2;
			double energy = om * (electronMass*cLight2);
			int count=0;
			while (count == 0 && jE < nE-1) {
				lim2 = sqrt(energies[jE]*energies[jE+1]);
				if ( energy < lim2 && energy > lim1) {
					double dnu = (lim2-lim1)/planck;
					lumOut1[jE] += lumOutAux * cWeights /dnu;
					count++;
				} else {
					lim1=lim2;
					jE++;
				}
			}
			om *= var_int_om;
		}
		sum2 += cWeights;
	}
	
    for (size_t jE=0;jE<nE;jE++) {
        lumOut[jE][jR] += lumOut1[jE] * dnuprim /sum2;
    }
}