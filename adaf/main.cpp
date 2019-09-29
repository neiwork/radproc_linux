#include <stdio.h>

#include "write.h"
#include "messages.h"
#include "modelParameters.h"
#include "State.h"
#include "comptonScattMatrix.h"
#include "thermalProcesses.h"
#include "globalVariables.h"

#include "thermalDistribution.h"
#include "radiativeLosses.h"
#include "injection.h"
#include "injectionNeutrons.h"
#include "distributionNeutrons.h"
#include "distribution.h"
#include "processes.h"
#include "absorption.h"
#include "blob.h"
#include "oneZoneTimeDependent.h"
#include "pairProcesses.h"

#include <fparameters/parameters.h>
#include <fparameters/SpaceIterator.h>
#include <inout/ioutil.h>
#include <iostream>
#include <boost/property_tree/ptree.hpp>

using namespace std;
int main()
{
	string folder{prepareOutputfolder()};
	try {
		GlobalConfig = readConfig();
		prepareGlobalCfg();
		
		show_message(msgStart, Module_state);
		State model(GlobalConfig.get_child("model"));
		show_message(msgEnd, Module_state);
		
		//tAccBlob = 0.0;
		//double tMax = 8000;
		//double dt = tMax/10;
		//for (int i=0;i<=10;i++) {
		
		//blob(model);

		if (calculateComptonScatt) {
			show_message(msgStart,Module_comptonScattMatrix);
			comptonScattMatrix(model);
			show_message(msgEnd,Module_comptonScattMatrix);
		} else {
			comptonScattMatrixRead(model);
		}
		
		//injection(model.ntElectron, model);
		//distributionFast(model.ntElectron, model);
		//writeEandRParamSpace("electronDis",model.ntElectron.distribution,0);
		
		if (calculateThermal)
			thermalProcesses(model,"lum.txt");
		
		if (calculateFlare) {
			timeAfterFlare = 0.0;
			double tMax = 4.0*3600;
			double dt = tMax/100;
			ofstream file;
			file.open("lightCurve.txt");
			for (int i=0;i<=100;i++) {
				oneZoneDist(model.ntPair,model);
				flareEmission(model,model.ntPair,file);
				timeAfterFlare += dt;
				if (i==1) writeEandRParamSpace("electronDis",model.ntPair.distribution,0);
			}
			file.close();
		}
		//blobEmission(model);
		
//***********nonthermal particles**************		
		
		if (calculateNonThermal) {
            
			if (calculateLosses) {
				show_message(msgStart, Module_radLosses);
				radiativeLosses(model.ntElectron, model, "electronLosses.txt");
				//radiativeLosses(model.ntProton, model, "protonLosses.txt");
				show_message(msgEnd, Module_radLosses);
			}
		
			if (calculateNTdistributions) {

			//nt electrons

				if (calculateFlare)
					injectionBurst(model.ntElectron,model);
				else
					injection(model.ntElectron, model);

				writeEandRParamSpace("electronInj",model.ntElectron.injection,0);
				distributionFast(model.ntElectron, model);
				//distributionSimplified(model.ntElectron,model);
				writeEandRParamSpace("electronDis",model.ntElectron.distribution,0);
				
				//nt protons
				/*injection(model.ntProton,model);
				writeEandRParamSpace("protonInj", model.ntProton.injection, 0);
				distributionFast(model.ntProton,model);
				writeEandRParamSpace("protonDis", model.ntProton.distribution, 0);
				*/
				if (calculateNeutronInj) {
					injectionNeutrons(model);
					radiativeLossesNeutron(model.ntNeutron,model,"neutronLosses.txt");
					if (calculateNeutronDis)
						distributionNeutrons(model);
					if (calculateJetDecay)
						jetNeutronDecay(model);
				}
				
				if (calculateNonThermalLum) {
					processes(model,"ntLuminosities.txt");
				}
				
				/*
				injectionChargedPion(model.ntChargedPion,model);
				writeEandRParamSpace("pionInj",model.ntChargedPion.injection,0);
				distributionFast(model.ntChargedPion,model);
				writeEandRParamSpace("pionDist",model.ntChargedPion.distribution,0);
				injectionMuon(model.ntMuon,model);
				writeEandRParamSpace("muonInj",model.ntMuon.injection,0);
				distributionFast(model.ntMuon,model);
				writeEandRParamSpace("muonDist",model.ntMuon.distribution,0);

				injectionPair(model.ntPair,model);
				writeEandRParamSpace("secondaryPairsInj",model.ntPair.injection,0);
				//writeEandRParamSpace("secondaryPairsInj2",model.ntPair.distribution,0);
					
				distributionFast(model.ntPair, model);
				writeEandRParamSpace("secondaryPairsDist",model.ntPair.distribution,0);
				pairProcesses(model, "ntPairtLum.txt");
				*/
			}
		}
	}
	catch (std::runtime_error& e)
	{
		std::cout << "ERROR: " << e.what() << std::endl;
	}
	return 0;
}




/* prueba interpol

		model.ntElectron.ps.iterate([&](const SpaceIterator& i) {
		
			double E = i.val(0);
			double u = model.ntElectron.emin()*10.0;
			
			double Qe = model.ntElectron.injection.interpolate({ { 0, u } }, &i.coord); 
			
			double res = Qe;
			
			
		});*/