#include "write.h"
#include "State.h"
#include "modelParameters.h"
#include "globalVariables.h"
#include "adafFunctions.h"
#include <fparameters/SpaceIterator.h>
#include <fparameters/Dimension.h>
#include <fparameters/parameters.h>
#include <fmath/RungeKutta.h>
#include <boost/property_tree/ptree.hpp>

//namespace {
	double safeLog10( double x ) {
		return x>0.0 ? log10(x) : -100.0;
	}
//}

std::string dataName(std::string id) {
	return id + ".dat";
}

void generateViewScript(std::string path) {
	std::string filename = path.substr(path.find("\\") + 1);
	std::string folder = path.substr(0,path.find("\\"));
	std::ofstream file;
	file.open((folder+"/plots/plot-"+filename+".bat").c_str(), std::ios::out);
	file << "@../../plot-svg-and-view.bat " + filename;
	file.close();
}


void writeAllSpaceParam(const std::string& filename, const ParamSpaceValues& data)
{
	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);
	//const ParamSpace* a = &(data.ps);
	
	// a.iterate;  no me deja hacer esta operacion
	
	data.ps.iterate([&file, &data](const SpaceIterator& i){
		double logE = log10(i.val(DIM_E) / 1.6e-12);
		double logR = (i.val(DIM_R));
		double logT = (i.val(DIM_Rcd));
		double logQ = safeLog10(data.get(i)); //log10(salida.values(i));  // values(i));
//		salida.values(i);


		file << logE << '\t' << logR << '\t' << logT << '\t' << 
			logQ << std::endl;
			//logQ << std::endl;
	});

	file.close();
	generateViewScript(filename);
}

void writeEandRParamSpace(const std::string& filename, const ParamSpaceValues& data, int t,int vol)
{

	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);

	// version acotada
	//double time = log10(data.ps[1][t]);
	
	data.ps.iterate([&](const SpaceIterator& i) {
		
		double logE = log10(i.val(DIM_E) / EV_TO_ERG);
		double r = i.val(DIM_R);
		double voll = (vol == 1) ? volume(r) : 1.0;
		double logQ = safeLog10(data.get(i)*voll);
		file << logE << '\t' << i.coord[DIM_R] << '\t' << safeLog10(r/schwRadius) << '\t' << logQ << std::endl;
			
	}, { -1, -1, t });  

	file.close();
	generateViewScript(filename);
}

void writeRParamSpace(const std::string& filename, const ParamSpaceValues& data, int t, int s)
{

	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);
	double Emin = data.ps[DIM_E].first();
	double Emax = data.ps[DIM_E].last();
	data.ps.iterate([&](const SpaceIterator& iR) {
		double tot = integSimpson(log(Emin),log(Emax),
						[&](double loge)
						{
							double e = exp(loge);
							double nPh = data.interpolate({{DIM_E,e}},&iR.coord);
							return e*nPh;
						},100);
		file << safeLog10(iR.val(DIM_R)/schwRadius) << '\t' << safeLog10(tot) << std::endl;
	},{t,-1,s});  
	file.close();
	generateViewScript(filename);
}

void writeEandTParamSpace(const std::string& filename, const ParamSpaceValues& data, int r)
{

	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);

	// version acotada
	double logR = log10(data.ps[1][r]);
	
	file << "log(r)=" << logR << '\t' ;

	for (size_t t_ix = 0; t_ix < data.ps[2].size(); t_ix++) {
		double time = data.ps[2][t_ix];
		file << "t=" << log10(time) << '\t';
	}


	for (int E_ix = 0; E_ix < data.ps[0].size(); E_ix++) {

		file << std::endl;

		double logE = log10(data.ps[0][E_ix]/1.6e-12);

		file << logE << '\t';

		data.ps.iterate([&file, &data](const SpaceIterator& i){

			//double logR = log10(i.val(DIM_R));
			//double time = i.val(DIM_T);
			double logQ = safeLog10(data.get(i));

			file << logQ << '\t';
			;
		}, { E_ix, r, -1 });  //el -1 indica que las E se recorren, no quedan fijas
		//las otras dos dimensiones quedan fijas en las posiciones r y t (recordar que la primera es 0 )
	}
	file.close();
	generateViewScript(filename);
}


void writeRandTParamSpace(const std::string& filename, const ParamSpaceValues& data, int E)
{

	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);


	data.ps.iterate([&file, &data](const SpaceIterator& i){

		double r = i.val(DIM_R);
		//double theta = i.val(DIM_THETA);
		
		file << r << '\t' << '\t' << data.get(i) << std::endl;
		
	}, { E, -1, -1 });  //el -1 indica que las E se recorren, no quedan fijas
	//las otras dos dimensiones quedan fijas en las posiciones r y t (recordar que la primera es 0 )
	
	file.close();
	generateViewScript(filename);
}


void writeEnergyFunction(const std::string& filename, const ParamSpaceValues& data, int r, int t)
{

	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);

	// version acotada
	double logR = log10(data.ps[1][r]);
	// version larga
	double logT = log10(data.ps.dimensions[2]->values[t]);

	file << "log(r)=" << logR << '\t' << "log(t)=" << logT << std::endl;
	data.ps.iterate([&file, &data](const SpaceIterator& i){

		double logE = log10(i.val(DIM_E) / 1.6e-12);
		double logQ = safeLog10(data.get(i));

		file << logE << '\t' << logQ << std::endl;
		;
	}, { -1, r, t });  //el -1 indica que las E se recorren, no quedan fijas
	//las otras dos dimensiones quedan fijas en las posiciones r y t (recordar que la primera es 0 )

	file.close();
	generateViewScript(filename);
}




void writeMatrix(const std::string& filename, Particle& p, Matrix& a)
{
	std::ofstream file;
	file.open(dataName(filename).c_str(), std::ios::out);

	int nR = p.ps[DIM_R].size();  //ver el -1

	file << '\t';
	
	for (size_t z_j = 0; z_j < nR; z_j++) { 
		const double r_j = p.ps[DIM_R][z_j];
		file << r_j << '\t' ;
	}
	
	file << std::endl; 
			
	for (size_t z_i = 0; z_i < nR; z_i++) { 
		
		const double r_i = p.ps[DIM_R][z_i];
				
		file << r_i << '\t';
		
		for (size_t z_j = 0; z_j < nR; z_j++) { 
			const double r_j = p.ps[DIM_R][z_j];
	
			file << a[z_i][z_j] << '\t' ;
		}
		
		
		
		file << std::endl;  
	}
	
	file.close();
	generateViewScript(filename);
}

void writeFields(State& st) {
	std::ofstream fields;
	fields.open("fields.dat",std::ios::out);
	fields  << "r [Rs]" 	<< "\t"
			<< "MdotRIAF" 	<< "\t"
			<< "MdotCD" 	<< "\t"
			<< "Te"			<< "\t"
			<< "Ti"			<< "\t"
			<< "H/R"		<< "\t"
			<< "ne"			<< "\t"
			<< "v/c"		<< "\t"
			<< "B"			<< endl;

	st.photon.ps.iterate([&](const SpaceIterator& iR) {
		double r = iR.val(DIM_R);
		fields << r/schwRadius << "\t"
			   << gAcc(r) << "\t"
			   << accRateColdDisk(r)/accRateOut << "\t"
			   << st.tempElectrons.get(iR) << "\t"
			   << st.tempIons.get(iR) << "\t"
			   << height_fun(r)/r << "\t"
			   << electronDensity(r) << "\t"
			   << abs(radialVel(r))/cLight <<"\t"
			   << st.magf.get(iR) << endl;
	},{0,-1,0});
	fields.close();
}