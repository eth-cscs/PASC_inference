/** @file test_dot.cu
 *  @brief test the speed of dot product computation
 *
 *  @author Lukas Pospisil
 */

#include <iostream>
#include <list>
#include <algorithm>

#include "pascinference.h"

typedef petscvector::PetscVector PetscVector;

using namespace pascinference;

extern int pascinference::DEBUG_MODE;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_filename1", boost::program_options::value< std::string >(), "name of input file with vector1 data (vector in PETSc format) [string]")
		("test_filename2", boost::program_options::value< std::string >(), "name of input file with vector2 data (vector in PETSc format) [string]")
		("test_n", boost::program_options::value< int >(), "how many times to test operation [int]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	/* load console arguments */
	std::string filename1;
	std::string filename2;
	int n;
	consoleArg.set_option_value("test_filename1", &filename1, "data/small_data.bin");
	consoleArg.set_option_value("test_filename2", &filename2, "data/small_solution.bin");
	consoleArg.set_option_value("test_n", &n, 100);

	/* print settings */
	coutMaster << " test_filename1      = " << std::setw(30) << filename1 << " (name of input file with vector1 data)" << std::endl;
	coutMaster << " test_filename2      = " << std::setw(30) << filename2 << " (name of input file with vector2 data)" << std::endl;
	coutMaster << " test_n              = " << std::setw(30) << n << " (how many times to test operation)" << std::endl;

	coutMaster << std::endl;

	/* load vector */
	Timer mytimer;
	mytimer.restart();
	mytimer.start();

	Vec myvector1_Vec;
	Vec myvector2_Vec;
	GeneralVector<PetscVector> myvector1(myvector1_Vec);
	myvector1.load_global(filename1);
	GeneralVector<PetscVector> myvector2(myvector2_Vec);
	myvector2.load_global(filename2);

	mytimer.stop();
	coutMaster << "- time load      : " << mytimer.get_value_last() << " s" << std::endl;

	coutMaster << "- dimension of v1: " << myvector1.size() << std::endl;
	coutMaster << "- dimension of v2: " << myvector2.size() << std::endl;

	
	/* we will measure the time of operation */
	double dot_result;
	for(int i=0;i<n;i++){
		dot_result = dot(myvector1,myvector2);
		#ifdef USE_CUDA
			gpuErrchk( cudaDeviceSynchronize() );
		#endif
		TRYCXX( PetscBarrier(NULL) );
		MPI_Barrier(MPI_COMM_WORLD);

	}
	mytimer.stop();


	coutMaster << "- result         : " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << dot_result << std::endl;
	coutMaster << "- time total     : " << mytimer.get_value_last() << " s" << std::endl;
	coutMaster << "- time average   : " << mytimer.get_value_last()/(double)n << " s" << std::endl;

	coutMaster << std::endl;

	Finalize();

	return 0;
}


