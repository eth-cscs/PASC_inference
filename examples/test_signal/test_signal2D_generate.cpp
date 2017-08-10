/** @file test_signal2D_generate.cpp
 *  @brief generate testing k-means sample
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"

#include <vector>
#include <random>

#ifndef USE_PETSC
 #error 'This example is for PETSC'
#endif

#define DEFAULT_TPERIOD 10
#define DEFAULT_REPEAT_NMB 5
#define DEFAULT_NOISE 0.1
#define DEFAULT_K 3
#define DEFAULT_DATA_TYPE 0

#define DEFAULT_FILENAME_DATA "data/test_signal/signal2D_data.bin"
#define DEFAULT_FILENAME_SOLUTION "data/test_signal/signal2D_solution.bin"
#define DEFAULT_FILENAME_GAMMA0 "data/test_signal/signal2D_gamma0.bin"
#define DEFAULT_GENERATE_DATA true
#define DEFAULT_GENERATE_GAMMA0 true

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

int myget_cluster_id(int t, int Tperiod){
	/* compute t in period */
	int tperiod = t - (int)(t/(double)Tperiod)*Tperiod;
	return myget_cluster_id_period(tperiod,Tperiod);
}

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_filename_data", boost::program_options::value< std::string >(), "name of output file with signal data (vector in PETSc format) [string]")
		("test_filename_solution", boost::program_options::value< std::string >(), "name of output file with original signal data without noise (vector in PETSc format) [string]")
		("test_filename_gamma0", boost::program_options::value< std::string >(), "name of output file with initial gamma approximation (vector in PETSc format) [string]")
		("test_Tperiod", boost::program_options::value< int >(), "length of one period of time-series [int]")
		("test_repeat_nmb", boost::program_options::value< int >(), "number of periods in time-series [int]")
		("test_K", boost::program_options::value< int >(), "number of clusters for gamma0 [int]")
		("test_data_type", boost::program_options::value< int >(), "type of output vector [0=TRn, 1=TnR, 2=nTR]")
		("test_noise", boost::program_options::value< double >(), "parameter of noise [double]")
		("test_generate_data", boost::program_options::value< bool >(), "generate solution and data with noise [bool]")
		("test_generate_gamma0", boost::program_options::value< bool >(), "generate gamma0 [bool]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize<PetscVector>(argc, argv)){
		return 0;
	}

	/* start to measure time */
    Timer timer_all;
    timer_all.start();

	int K;
	int Tperiod;
	int repeat_nmb;
	double noise;
    int data_type;
	std::string filename_data;
	std::string filename_solution;
	std::string filename_gamma0;
	bool generate_data, generate_gamma0;

	consoleArg.set_option_value("test_filename_data", &filename_data, DEFAULT_FILENAME_DATA);
	consoleArg.set_option_value("test_filename_solution", &filename_solution, DEFAULT_FILENAME_SOLUTION);
	consoleArg.set_option_value("test_filename_gamma0", &filename_gamma0, DEFAULT_FILENAME_GAMMA0);
	consoleArg.set_option_value("test_Tperiod", &Tperiod, DEFAULT_TPERIOD);
	consoleArg.set_option_value("test_K", &K, DEFAULT_K);
	consoleArg.set_option_value("test_data_type", &data_type, DEFAULT_DATA_TYPE);
	consoleArg.set_option_value("test_repeat_nmb", &repeat_nmb, DEFAULT_REPEAT_NMB);
	consoleArg.set_option_value("test_noise", &noise, DEFAULT_NOISE);
	consoleArg.set_option_value("test_generate_data", &generate_data, DEFAULT_GENERATE_DATA);
	consoleArg.set_option_value("test_generate_gamma0", &generate_gamma0, DEFAULT_GENERATE_GAMMA0);

	/* set decomposition in space */
	int DDT_size = GlobalManager.get_size();

	coutMaster << "- PROBLEM INFO ----------------------------" << std::endl;
	coutMaster << " test_Tperiod                = " << std::setw(50) << Tperiod << " (length of one period of time-series)" << std::endl;
	coutMaster << " test_repeat_nmb             = " << std::setw(50) << repeat_nmb << " (number of periods in time-series)" << std::endl;
	coutMaster << " test_K                      = " << std::setw(50) << K << " (number of clusters for gamma0)" << std::endl;
	coutMaster << " test_data_type              = " << std::setw(50) << Decomposition<PetscVector>::get_type_name(data_type) << " (type of output vector [" << Decomposition<PetscVector>::get_type_list() << "])" << std::endl;
	coutMaster << " test_noise                  = " << std::setw(50) << noise << " (parameter of noise)" << std::endl;
	coutMaster << " test_filename_data          = " << std::setw(50) << filename_data << " (name of input file with signal data)" << std::endl;
	coutMaster << " test_filename_solution      = " << std::setw(50) << filename_solution << " (name of input file with original signal data without noise)" << std::endl;
	coutMaster << " test_filename_gamma0        = " << std::setw(50) << filename_gamma0 << " (name of input file with initial gamma approximation)" << std::endl;
	coutMaster << " test_generate_data          = " << std::setw(50) << generate_data << " (generate solution and data with noise)" << std::endl;
	coutMaster << " test_generate_gamma0        = " << std::setw(50) << generate_gamma0 << " (generate gamma0)" << std::endl;
	coutMaster << "-------------------------------------------" << std::endl;

	/* say hello */
	coutMaster << "- start program" << std::endl;

	int T = Tperiod*repeat_nmb;
	int xdim = 2;
	int K_true = 3;
	double mu[xdim*K_true] = {1.0,1.0, 2.0,2.0, 3.0,1.5};
	coutMaster << " T                           = " << std::setw(50) << T << " (length of time-series)" << std::endl;
	coutMaster << " mu                          = " << std::setw(50) << print_array(mu,6) << " (mean values in cluster [mu1_x,mu1_y,mu2_x,mu2_y,mu3_x,mu3_y])" << std::endl;

	/* allocate vector of data */
	Vec x;
	TRYCXX( VecCreate(PETSC_COMM_WORLD,&x) );
	TRYCXX( VecSetSizes(x,xdim*T,xdim*T) );
	TRYCXX( VecSetType(x, VECSEQ) );
	TRYCXX( VecSetFromOptions(x) );

	/* vector for data with noise */
	int cluster_id;
	if(generate_data){
		std::default_random_engine generator(time(0));
		std::normal_distribution<double> distribution(0.0,noise);

		Vec x_data;
		TRYCXX( VecDuplicate(x,&x_data) );

		/* fill array */
		double *x_arr;
		double *x_data_arr;
		TRYCXX( VecGetArray(x,&x_arr));
		TRYCXX( VecGetArray(x_data,&x_data_arr));
		for(int t=0;t<T;t++){
			for(int n=0;n<xdim;n++){
                double noise_value = distribution(generator);

                int idx;
                if(data_type == 0 | data_type==1){
                    /* type TRn,TnR with R=1 -> Tn */
                    idx = t*xdim + n;
                }
                if(data_type==2){
                    /* type TRn,TnR with R=1 -> nT */
                    idx = n*T + t;
                }

				cluster_id = myget_cluster_id(t, Tperiod);
				x_arr[idx] = mu[cluster_id*xdim + n];
				x_data_arr[idx] = x_arr[idx] + noise_value;
			}
		}
		TRYCXX( VecRestoreArray(x,&x_arr));
		TRYCXX( VecRestoreArray(x_data,&x_data_arr));

		/* save generated vectors */
		PetscViewer mviewer;
		TRYCXX( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
		TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename_solution.c_str(), FILE_MODE_WRITE, &mviewer) );
		TRYCXX( VecView(x, mviewer) );
		TRYCXX( PetscViewerDestroy(&mviewer) );
		coutMaster << " - new solution vector saved" << std::endl;

		TRYCXX( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer) );
		TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename_data.c_str(), FILE_MODE_WRITE, &mviewer) );
		TRYCXX( VecView(x_data, mviewer) );
		TRYCXX( PetscViewerDestroy(&mviewer) );
		coutMaster << " - new data vector with noise saved" << std::endl;

	}

	if(generate_gamma0){
		/* vector for gamma0 */
		Vec gamma0;
		TRYCXX( VecCreate(PETSC_COMM_WORLD,&gamma0) );
		TRYCXX( VecSetSizes(gamma0,T*K,T*K) );
		TRYCXX( VecSetType(gamma0, VECSEQ) );
		TRYCXX( VecSetFromOptions(gamma0) );

		/* generate random gamma0 */
		PetscRandom rctx2;
		TRYCXX( PetscRandomCreate(PETSC_COMM_WORLD,&rctx2) );
		TRYCXX( PetscRandomSetFromOptions(rctx2) );
		TRYCXX( VecSetRandom(gamma0,rctx2) );
		TRYCXX( PetscRandomDestroy(&rctx2) );

		/* save generated vector */
		PetscViewer mviewer2;
		TRYCXX( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer2) );
		TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename_gamma0.c_str(), FILE_MODE_WRITE, &mviewer2) );
		TRYCXX( VecView(gamma0, mviewer2) );
		TRYCXX( PetscViewerDestroy(&mviewer2) );
		coutMaster << " - new random gamma0 vector saved" << std::endl;

	}

	/* say bye */
	coutMaster << "- end program" << std::endl;

	Finalize<PetscVector>();

	return 0;
}

