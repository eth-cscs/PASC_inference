/** @file test_integration.cpp
 *  @brief test the EntropyIntegration class
 *
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include <vector>

#ifndef USE_PETSC
 #error 'This example is for PETSC'
#endif

#ifndef USE_CUBA
 #error 'This example is for CUBA'
#endif

#ifndef USE_DLIB
 #error 'This example is for DLIB'
#endif

#define DEFAULT_XDIM 1
#define DEFAULT_TYPE 1
#define DEFAULT_PRINTINFO true
#define DEFAULT_FILENAME_IN "data/test_image/usi_text/usi_250_150_02.bin"
#define DEFAULT_EPS 1e-6
#define DEFAULT_KM 2

using namespace pascinference;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("TEST ENTROPYINTEGRATION", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_filename_in", boost::program_options::value< std::string >(), "name of input file with image data (vector in PETSc format) [string]")
		("test_eps", boost::program_options::value<double>(), "integration precision [double]")
		("test_xdim", boost::program_options::value<int>(), "number of values in every pixel [1=greyscale, 3=rgb]")
		("test_type", boost::program_options::value< int >(), "type of integration [0=Dlib, 1=Cuba]")
		("test_Km", boost::program_options::value< int >(), "number of moments [int]")


		("test_printinfo", boost::program_options::value<bool>(), "print informations about created objects [bool]");

	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	}

	int xdim, type, Km;
	double eps;
	bool printinfo;

	std::string filename_in;

	consoleArg.set_option_value("test_xdim", &xdim, DEFAULT_XDIM);
	consoleArg.set_option_value("test_Km", &Km, DEFAULT_KM);
	consoleArg.set_option_value("test_type", &type, DEFAULT_TYPE);
	consoleArg.set_option_value("test_printinfo", &printinfo, DEFAULT_PRINTINFO);
	consoleArg.set_option_value("test_filename_in", &filename_in, DEFAULT_FILENAME_IN);
	consoleArg.set_option_value("test_eps", &eps, DEFAULT_EPS);

	coutMaster << "- PROBLEM INFO ----------------------------" << std::endl;
#ifdef USE_CUDA
	coutMaster << " computing on GPU" << std::endl;
#else
	coutMaster << " computing on CPU" << std::endl;
#endif
	coutMaster << " test_filename_in            = " << std::setw(50) << filename_in << " (name of input file with image data)" << std::endl;

	coutMaster << " test_xdim                   = " << std::setw(50) << xdim << " (dimension of variables)" << std::endl;
	coutMaster << " test_Km                     = " << std::setw(50) << Km << " (number of moments)" << std::endl;
	coutMaster << " test_eps                    = " << std::setw(50) << eps << " (integration precision)" << std::endl;

	coutMaster << " test_type                   = " << std::setw(50) << type << " (type of integration [0=Dlib/1=Cuba])" << std::endl;
	coutMaster << " test_printinfo              = " << std::setw(50) << print_bool(printinfo) << " (print informations about created objects)" << std::endl;

	coutMaster << "-------------------------------------------" << std::endl;
	coutMaster << std::endl;

	/* this test is only for one cluster */
	int K = 1;
	int T = 10;

	/* start logging */
	std::ostringstream oss;
	oss << "log/test_integration.txt";
	logging.begin(oss.str());
	oss.str("");

	/* say hello */
	coutMaster << "- start program" << std::endl;

/* EntropyData */
	coutMaster << "--- PREPARING ENTROPYDATA ---" << std::endl;
	EntropyData<PetscVector> *entropydata;

	entropydata = new EntropyData<PetscVector>(T, xdim, K, Km);

	if(printinfo) entropydata->print(coutMaster);

/* EntropyIntegration */
	coutMaster << "--- PREPARING ENTROPYINTEGRATION ---" << std::endl;
	EntropyIntegration<PetscVector> *entropyintegration;
	if(type == 0){
		entropyintegration = new EntropyIntegrationDlib<PetscVector>(entropydata, eps);
	}
	if(type == 1){
		entropyintegration = new EntropyIntegrationCuba<PetscVector>(entropydata, eps);
	}

	if(printinfo) entropyintegration->print(coutMaster);

	/* say bye */
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize<PetscVector>();

	return 0;
}

