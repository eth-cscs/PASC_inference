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

#define DEFAULT_K 1
#define DEFAULT_DATA_TYPE 2
#define DEFAULT_XDIM 1
#define DEFAULT_TYPE 1
#define DEFAULT_PRINTINFO true
#define DEFAULT_FILENAME_IN "data/entropy_small_data.bin"
#define DEFAULT_EPS 1e-6
#define DEFAULT_KM 2

using namespace pascinference;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("TEST ENTROPYINTEGRATION", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_filename_in", boost::program_options::value< std::string >(), "name of input file with signal data (vector in PETSc format) [string]")
		("test_eps", boost::program_options::value<double>(), "integration precision [double]")
		("test_xdim", boost::program_options::value<int>(), "dimension of data [int]")
		("test_data_type", boost::program_options::value< int >(), "type of input/output vector [0=TRn, 1=TnR, 2=nTR]")
		("test_K", boost::program_options::value<int>(), "number of clusters [int]")
		("test_type", boost::program_options::value< int >(), "type of integration [0=Dlib, 1=Cuba]")
		("test_Km", boost::program_options::value< int >(), "number of moments [int]")


		("test_printinfo", boost::program_options::value<bool>(), "print informations about created objects [bool]");

	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	}

	int xdim, type, Km, K, data_type;
	double eps;
	bool printinfo;

	std::string filename_in;

	consoleArg.set_option_value("test_xdim", &xdim, DEFAULT_XDIM);
	consoleArg.set_option_value("test_data_type", &data_type, DEFAULT_DATA_TYPE);
	consoleArg.set_option_value("test_K", &K, DEFAULT_K);
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
	coutMaster << " test_data_type              = " << std::setw(50) << Decomposition<PetscVector>::get_type_name(data_type) << " (type of output vector [" << Decomposition<PetscVector>::get_type_list() << "])" << std::endl;

	coutMaster << " test_xdim                   = " << std::setw(50) << xdim << " (dimension of variables)" << std::endl;
	coutMaster << " test_K                      = " << std::setw(50) << K << " (number of clusters)" << std::endl;
	coutMaster << " test_Km                     = " << std::setw(50) << Km << " (number of moments)" << std::endl;
	coutMaster << " test_eps                    = " << std::setw(50) << eps << " (integration precision)" << std::endl;

	coutMaster << " test_type                   = " << std::setw(50) << type << " (type of integration [0=Dlib/1=Cuba])" << std::endl;
	coutMaster << " test_printinfo              = " << std::setw(50) << print_bool(printinfo) << " (print informations about created objects)" << std::endl;

	coutMaster << "-------------------------------------------" << std::endl;
	coutMaster << std::endl;

	/* start logging */
	std::ostringstream oss;
	oss << "log/test_integration.txt";
	logging.begin(oss.str());
	oss.str("");

	/* say hello */
	coutMaster << "- start program" << std::endl;

/* prepare preliminary time-series data (to get the size of the problem T) */
	coutMaster << "--- PREPARING PRELIMINARY DATA ---" << std::endl;
	SignalData<PetscVector> mydata(filename_in);

/* Decomposition in time */
	Decomposition<PetscVector> decomposition(mydata.get_Tpreliminary()/(double)xdim, 1, K, xdim, GlobalManager.get_size());

	/* print info about decomposition */
	if(printinfo) decomposition.print(coutMaster);

/* prepare time-series data */
	coutMaster << "--- APPLY DECOMPOSITION TO DATA ---" << std::endl;
	mydata.set_decomposition(decomposition, data_type);

	/* print information about loaded data */
	if(printinfo) mydata.print(coutMaster);

/* EntropyData */
	coutMaster << "--- PREPARING ENTROPYDATA ---" << std::endl;
	EntropyData<PetscVector> *entropydata;

	entropydata = new EntropyData<PetscVector>(&decomposition, Km);

	/* set data to entropy */
	entropydata->set_x(mydata.get_datavector());

	/* for testing purposes create Gamma=1 - this step is typically performed by model */
	Vec gamma_Vec;
	decomposition.createGlobalVec_gamma(&gamma_Vec);
	GeneralVector<PetscVector> gammavector(gamma_Vec);
    TRYCXX( VecSet(gammavector.get_vector(),1.0));
	entropydata->set_gamma(&gammavector);

    /* prepare vector where store computed lambda - this step is typically performed by model */
	Vec lambda_Vec;
	TRYCXX( VecCreate(PETSC_COMM_SELF,&lambda_Vec) );
	#ifdef USE_CUDA
		TRYCXX(VecSetType(lambda_Vec, VECSEQCUDA));
	#else
		TRYCXX(VecSetType(lambda_Vec, VECSEQ));
	#endif
	TRYCXX( VecSetSizes(lambda_Vec, K*entropydata->get_number_of_moments(), PETSC_DECIDE) );
	TRYCXX( VecSetFromOptions(lambda_Vec) );
	GeneralVector<PetscVector> lambda(lambda_Vec);
    entropydata->set_lambda(&lambda);

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


/* ----------- PERFORM INTEGRATION --------- */
    entropyintegration->compute(   );

	/* say bye */
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize<PetscVector>();

	return 0;
}
