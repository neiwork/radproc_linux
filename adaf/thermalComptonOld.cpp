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

void thermalCompton2(const State& st, Matrix& lumOut, Matrix scattADAF, Matrix scattCD, 
					Matrix absCD, Matrix& lumOutIC, Vector energies)
{
	show_message(msgStart,Module_thermalCompton);
	
	size_t nOm = 20;
	size_t nG = 6;
	Vector p(nR*nG*nE*nOm,0.0);
	
	size_t jR = 0;
	if (calculateProbs) {
		st.photon.ps.iterate([&](const SpaceIterator& itR) {
			double temp = st.tempElectrons.get(itR);
			double normtemp = boltzmann*temp/(electronMass*cLight2);
			if (normtemp >= 0.1) {
				comptonRedistribution(p,nG,nE,nOm,jR,normtemp,energies,lumOut);
			}
			jR++;
		},{0,-1,0});
		vectorWrite("comptonProb.txt",p,nR*nG*nE*nOm);
	} else {
		vectorRead("comptonProb.txt",p,nR*nG*nE*nOm);
	}
	
	Matrix lumOutLocal,lumCD,lumOutRefl;
	Vector lumInIC(nE,0.0);
	matrixInitCopy(lumOutLocal,nE,nR,lumOut);
	
	double res;
	size_t it=1;
	do {
		cout << "Iteration number = " << it << endl;
		res=0.0;

		Vector lumOld(nE,0.0);
		for (size_t jE=0;jE<nE;jE++) {
			for (size_t jR=0;jR<nR;jR++) {
				lumOld[jE] += lumOut[jE][jR];
			}
			/*if (it>1) {
				for (size_t jRcd=0;jRcd<nRcd;jRcd++) {
					lumOld[jE] += lumCD[jE][jRcd]+lumOutRefl[jE][jRcd];
				}
			}*/
		}

		fill(lumInIC.begin(),lumInIC.end(),0.0);
		matrixInit(lumOutIC,nE,nR,0.0);
		
		//coldDiskLuminosity(st,absCD,lumOut,lumOutRefl,lumCD,energies);
		
		// PARA CADA CELDA
		size_t jR=0;
		st.photon.ps.iterate([&](const SpaceIterator& itR) {
			// CALCULAMOS LinC PARA CADA E        	 function LinC
			for (size_t jE=0;jE<nE;jE++) {
				for (size_t jjR=0;jjR<nR;jjR++)
					lumInIC[jE] += scattADAF[jjR][jR]*lumOut[jE][jjR];
				//for (size_t jjRcd=0;jjRcd<nRcd;jjRcd++)
				//	lumInIC[jE] += scattCD[jjRcd][jR]*lumCD[jE][jjRcd];
			}
			double temp = st.tempElectrons.get(itR);
			double normtemp = boltzmann*temp/(electronMass*cLight2);
			if (normtemp >= 0.1) 
			{
				for (size_t jjE=0;jjE<nE;jjE++) {
					double frecuency=energies[jjE]/planck;
					if (lumInIC[jjE]*frecuency > 1.0e15)
						comptonNew2(lumOutIC,lumInIC[jjE],energies,nG,nE,nOm,normtemp,p,jjE,jR);
				}
			}
			jR++;
		},{0,-1,0});
		
		for (size_t jE=0;jE<nE;jE++) {
			double lumNew = 0.0;
			for (size_t jR=0;jR<nR;jR++) {
				double lumAux = lumOutLocal[jE][jR]+lumOutIC[jE][jR];
				lumNew += lumAux;
				lumOut[jE][jR] = lumAux;
			}
			/*for (size_t jRcd=0;jRcd<nRcd;jRcd++) {
				lumNew += lumCD[jE][jRcd]+lumOutRefl[jE][jRcd];
			}*/
			if (lumNew > 0.0 && lumOld[jE] > 0.0)
				res += abs(log10(lumNew/lumOld[jE]));
		}
		
		res /= nE;
		cout << "Residuo = " << res << endl;
		++it;
	} while (res > 1.0e-1);
	show_message(msgEnd,Module_thermalCompton);
}