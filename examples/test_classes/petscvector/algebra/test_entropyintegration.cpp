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
#define DEFAULT_TYPE 0
#define DEFAULT_XDIM 1
#define DEFAULT_PRINTINFO true
#define DEFAULT_PRINTCONTENT true
#define DEFAULT_EPS 1e-6
#define DEFAULT_KM 2

using namespace pascinference;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("TEST ENTROPYINTEGRATION", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_filename_in", boost::program_options::value< std::string >(), "name of input file with lambda values (vector in PETSc format) [string]")
		("test_eps", boost::program_options::value<double>(), "integration precision [double]")
		("test_xdim", boost::program_options::value<int>(), "dimension of data [int]")
		("test_K", boost::program_options::value<int>(), "number of clusters [int]")
		("test_type", boost::program_options::value< int >(), "type of integration [0=Dlib, 1=Cuba, 2=CudaVegas]")
		("test_Km", boost::program_options::value< int >(), "number of moments [int]")
		("test_printinfo", boost::program_options::value<bool>(), "print informations about created objects [bool]")
		("test_printcontent", boost::program_options::value<bool>(), "print computed integrals [bool]");

	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	}

	int xdim, type, Km, K;
	double eps;
	bool printinfo, printcontent, loadlambda;

	std::string filename_in;

	consoleArg.set_option_value("test_xdim", &xdim, DEFAULT_XDIM);
	consoleArg.set_option_value("test_K", &K, DEFAULT_K);
	consoleArg.set_option_value("test_Km", &Km, DEFAULT_KM);
	consoleArg.set_option_value("test_type", &type, DEFAULT_TYPE);
	consoleArg.set_option_value("test_printinfo", &printinfo, DEFAULT_PRINTINFO);
	consoleArg.set_option_value("test_printcontent", &printcontent, DEFAULT_PRINTCONTENT);
	consoleArg.set_option_value("test_eps", &eps, DEFAULT_EPS);

    if(!consoleArg.set_option_value("test_filename_in", &filename_in)){
        loadlambda = false;
    } else {
        loadlambda = true;
    }


	coutMaster << "- PROBLEM INFO ----------------------------" << std::endl;
#ifdef USE_CUDA
	coutMaster << " computing on GPU" << std::endl;
#else
	coutMaster << " computing on CPU" << std::endl;
#endif
	coutMaster << " test_filename_in            = " << std::setw(50) << filename_in << " (name of input file with lambda values (vector in PETSc format))" << std::endl;

	coutMaster << " test_xdim                   = " << std::setw(50) << xdim << " (dimension of variables)" << std::endl;
	coutMaster << " test_K                      = " << std::setw(50) << K << " (number of clusters)" << std::endl;
	coutMaster << " test_Km                     = " << std::setw(50) << Km << " (number of moments)" << std::endl;
	coutMaster << " test_eps                    = " << std::setw(50) << eps << " (integration precision)" << std::endl;

	coutMaster << " test_type                   = " << std::setw(50) << type << " (type of integration [0=Dlib/1=Cuba/2=CudaVegas])" << std::endl;
	coutMaster << " test_printinfo              = " << std::setw(50) << print_bool(printinfo) << " (print informations about created objects)" << std::endl;
	coutMaster << " test_printcontent           = " << std::setw(50) << print_bool(printcontent) << " (print computed moments and integrals)" << std::endl;

	coutMaster << "-------------------------------------------" << std::endl;
	coutMaster << std::endl;

	/* start logging */
	std::ostringstream oss;
	oss << "log/test_integration.txt";
	logging.begin(oss.str());
	oss.str("");

	/* say hello */
	coutMaster << "- start program" << std::endl;

/* Decomposition in time */
	Decomposition<PetscVector> decomposition(0, 0, K, xdim, GlobalManager.get_size());

	/* print info about decomposition */
	if(printinfo) decomposition.print(coutMaster);

/* EntropyData */
	coutMaster << "--- PREPARING ENTROPYDATA ---" << std::endl;
	EntropyData<PetscVector> *entropydata;

	entropydata = new EntropyData<PetscVector>(&decomposition, Km);

    if(printinfo){
        coutMaster << std::endl;
        coutMaster << "number of moments: " << entropydata->get_number_of_moments() << std::endl;
        coutMaster << std::endl;
    }

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

    if(loadlambda){
        lambda.load_global(filename_in);
    } else {
        TRYCXX( VecSet(lambda.get_vector(), 1.0) );
    }

    entropydata->set_lambda(&lambda);

	if(printinfo) entropydata->print(coutMaster);

    if(printcontent){
        coutMaster << "loaded lambda:" << std::endl;
        coutMaster << lambda << std::endl;
    }

/* EntropyIntegration */
	coutMaster << "--- PREPARING ENTROPYINTEGRATION ---" << std::endl;
	EntropyIntegration<PetscVector> *entropyintegration;
	if(type == 0){
		entropyintegration = new EntropyIntegrationDlib<PetscVector>(entropydata, eps);
	}
	if(type == 1){
		entropyintegration = new EntropyIntegrationCuba<PetscVector>(entropydata, eps);
	}
	if(type == 1){
		entropyintegration = new EntropyIntegrationCudaVegas<PetscVector>(entropydata, eps);
	}

	if(printinfo) entropyintegration->print(coutMaster);

    coutMaster << std::endl;
    coutMaster << "----------------------------------------------" << std::endl;
    coutMaster << std::endl;

/* ----------- prepare vector for integration --------- */
	Vec integrals_Vec;
	TRYCXX( VecCreate(PETSC_COMM_SELF,&integrals_Vec) );
	#ifdef USE_CUDA
		TRYCXX(VecSetType(integrals_Vec, VECSEQCUDA));
	#else
		TRYCXX(VecSetType(integrals_Vec, VECSEQ));
	#endif
	TRYCXX( VecSetSizes(integrals_Vec, K*entropyintegration->get_number_of_integrals(), PETSC_DECIDE) );
	TRYCXX( VecSetFromOptions(integrals_Vec) );
	GeneralVector<PetscVector> integrals(integrals_Vec);

/* ----------- PERFORM INTEGRATION --------- */
    Timer integration_timer;
    integration_timer.restart();
    integration_timer.start();
    entropyintegration->compute(integrals);
    integration_timer.stop();

    coutMaster << "   integrals computed in   : " << integration_timer.get_value_sum() << "s " << std::endl;

    if(printcontent){
        coutMaster << "   computed integrals:" << std::endl;
        coutMaster << integrals << std::endl;
    }

    coutMaster << std::endl;
    coutMaster << "----------------------------------------------" << std::endl;
    coutMaster << std::endl;

	/* say bye */
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize<PetscVector>();

	return 0;
}
