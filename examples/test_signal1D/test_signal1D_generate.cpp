/** @file test_signal1D_generate.cpp
 *  @brief generate testing k-means sample
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"

#include <vector>

#ifndef USE_PETSC
 #error 'This example is for PETSC'
#endif

#define DEFAULT_TPERIOD 10
#define DEFAULT_REPEAT_NMB 5
#define DEFAULT_NOISE 0.1
 
using namespace pascinference;

int myget_cluster_id_period(int tperiod, int Tperiod){
	/* returns id of cluster based on t */
	double step = Tperiod/11.0;
	int cluster_id=0;
	if( ( (tperiod >= 1*step) && (tperiod < 3*step) ) || ( (tperiod >= 4*step) && (tperiod < 5*step) ) || ( (tperiod >= 8*step) && (tperiod < 9*step) ) ) {
		cluster_id = 1;
	}
	if( ( (tperiod >= 3*step) && (tperiod < 4*step) ) || ( (tperiod >= 7*step) && (tperiod < 8*step) ) || (tperiod >= 9*step) ) {
		cluster_id = 2;
	}
	return cluster_id;
}

int myget_cluster_id_period(int t, int Tperiod){
	/* compute t in period */
	int tperiod = t - (int)(t/(double)Tperiod)*Tperiod;
	return myget_cluster_id_period(tperiod,Tperiod);
}

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_filename", boost::program_options::value< std::string >(), "name of output file with signal data (vector in PETSc format) [string]")
		("test_filename_solution", boost::program_options::value< std::string >(), "name of output file with original signal data without noise (vector in PETSc format) [string]")
		("test_filename_gamma0", boost::program_options::value< std::string >(), "name of output file with initial gamma approximation (vector in PETSc format) [string]")
		("test_T", boost::program_options::value< int >(), "length of one period of time-series [int]")
		("test_repeat_nmb", boost::program_options::value< int >(), "number of periods in time-series [int]")
		("test_noise", boost::program_options::value< double >(), "parameter of noise [double]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	}

	/* start to measure time */
    Timer timer_all;
    timer_all.start();

	int Tperiod;
	int repeat_nmb;
	double noise;
	std::string filename;
	std::string filename_solution;
	std::string filename_gamma0;

	consoleArg.set_option_value("test_filename", &filename, "data/samplesignal.bin");
	consoleArg.set_option_value("test_filename_solution", &filename_solution, "data/samplesignal_solution.bin");
	consoleArg.set_option_value("test_filename_gamma0", &filename_solution, "data/samplesignal_gamma0.bin");
	consoleArg.set_option_value("test_Tperiod", &Tperiod, DEFAULT_TPERIOD);
	consoleArg.set_option_value("test_repeat_nmb", &repeat_nmb, DEFAULT_REPEAT_NMB);
	consoleArg.set_option_value("test_noise", &noise, DEFAULT_NOISE);

	/* set decomposition in space */
	int DDT_size = GlobalManager.get_size();

	coutMaster << "- PROBLEM INFO ----------------------------" << std::endl;
	coutMaster << " ranks_per_node              = " << std::setw(30) << ranks_per_node << " (number of MPI processes on one node)" << std::endl;
	coutMaster << " test_Tperiod                = " << std::setw(30) << Tperiod << " (length of one period of time-series)" << std::endl;
	coutMaster << " test_repeat_nmb             = " << std::setw(30) << repeat_nmb << " (number of periods in time-series)" << std::endl;
	coutMaster << " test_noise                  = " << std::setw(30) << noise << " (parameter of noise)" << std::endl;
	coutMaster << " test_filename               = " << std::setw(30) << filename << " (name of input file with signal data)" << std::endl;
	coutMaster << " test_filename_solution      = " << std::setw(30) << filename_solution << " (name of input file with original signal data without noise)" << std::endl;
	coutMaster << " test_filename_gamma0        = " << std::setw(30) << filename_gamma0 << " (name of input file with initial gamma approximation)" << std::endl;
	coutMaster << "-------------------------------------------" << std::endl;


	/* start logging */
	std::ostringstream oss;
	oss << "log/" << filename_out << ".txt";
	logging.begin(oss.str());
	oss.str("");

	/* say hello */
	coutMaster << "- start program" << std::endl;

	int K = 3; /* number of cluster, in this example fixed */
	int T = Tperiod*repeat_nmb;
	double mu[3] = {0.0,1.0,2.0};
	coutMaster << " T                           = " << std::setw(30) << T << " (length of time-series)" << std::endl;
	coutMaster << " mu                          = " << std::setw(30) << print_array(mu,3) << " (mean values in clusters)" << std::endl;


	/* allocate vector of data */
	Vec x;
	TRYCXX( VecCreate(PETSC_COMM_WORLD,&x) );
	TRYCXX( VecSetSizes(x,PETSC_DECIDE,T) );
	TRYCXX( VecSetType(x, VECMPI) );
	TRYCXX( VecSetFromOptions(x) );

	/* ownership range */
	int low, high;
	TRYCXX( VecGetOwnershipRange(x, &low, &high) );

	/* vector for data with noise */
	Vec x_data;
	TRYCXX( VecDuplicate(x,&x_data) );

	/* vector for gamma0 */
	Vec gamma0;
	TRYCXX( VecCreate(PETSC_COMM_WORLD,&gamma0) );
	TRYCXX( VecSetSizes(gamma0,(high-low)*K,T*K) );
	TRYCXX( VecSetType(gamma0, VECMPI) );
	TRYCXX( VecSetFromOptions(gamma0) );
	
	/* add noise */
	PetscRandom rctx;
	TRYCXX( PetscRandomCreate(PETSC_COMM_WORLD,&rctx) );
	TRYCXX( PetscRandomSetFromOptions(rctx) );
	TRYCXX( VecSetRandom(x_data,rctx) );
	TRYCXX( VecSetRandom(gamma0,rctx) );
	TRYCXX( PetscRandomDestroy(&rctx) );

	/* fill local array */
	double x_arr;
	TRYCXX( VecGetArray(x,&x_arr));
	for(int t=low;t<high;t++){
		x_arr[t-low] = mu[myget_cluster_id_period(t, Tperiod)];
		x_arr[t-low] = mu[myget_cluster_id_period(t, Tperiod)];
	}
	TRYCXX( VecRestoreArray(x,&x_arr));

	TRYCXX( VecView(x, PETSC_VIEWER_STDOUT_WORLD) );
	TRYCXX( VecView(x_data, PETSC_VIEWER_STDOUT_WORLD) );
	
	/* add solution data to noise to create x_data */
	TRYCXX( VecAYPX(x_data, noise, x) );
	
	/* save generated vectors */
	PetscViewer mviewer;
	TRYCXX( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
	TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename_solution, FILE_MODE_WRITE, &mviewer) );
	TRYCXX( VecView(x, mviewer) );
	TRYCXX( PetscViewerDestroy(&mviewer) );

	TRYCXX( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
	TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &mviewer) );
	TRYCXX( VecView(x_data, mviewer) );
	TRYCXX( PetscViewerDestroy(&mviewer) );

	TRYCXX( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
	TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename_gamma0, FILE_MODE_WRITE, &mviewer) );
	TRYCXX( VecView(gamma0, mviewer) );
	TRYCXX( PetscViewerDestroy(&mviewer) );


	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize<PetscVector>();

	return 0;
}

