/** @file test_petscvector_dot.cpp
 *  @brief test the speed of dot product computation
 *
 *  @author Lukas Pospisil
 */

#include <iostream>
#include <list>
#include <algorithm>

#include "pascinference.h"

/* decide which vector to use for computation */
typedef petscvector::PetscVector TestVector;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_filename1", boost::program_options::value< std::string >(), "name of input file with vector1 data (vector in PETSc format) [string]")
		("test_filename2", boost::program_options::value< std::string >(), "name of input file with vector2 data (vector in PETSc format) [string]")
		("test_n", boost::program_options::value< double >(), "how many times to test operation [int]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	/* load console arguments */
	std::string filename1;
	std::string filename2;
	double n;
	consoleArg.set_option_value("test_filename1", &filename1, "data/small_data.bin");
	consoleArg.set_option_value("test_filename2", &filename2, "data/small_solution.bin");
	consoleArg.set_option_value("test_n", &n, 100.0);

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

	TRYCXX( VecCreate(PETSC_COMM_WORLD,&myvector1_Vec) );
	TRYCXX( VecCreate(PETSC_COMM_WORLD,&myvector2_Vec) );

#ifdef USE_CUDA
	TRYCXX( VecSetType(myvector1_Vec, VECMPICUDA) );
	TRYCXX( VecSetType(myvector2_Vec, VECMPICUDA) );

#else
	TRYCXX( VecSetType(myvector1_Vec, VECMPI) );
	TRYCXX( VecSetType(myvector2_Vec, VECMPI) );
#endif

	/* prepare viewer to load from file */
	PetscViewer mviewer1;
	PetscViewer mviewer2;
	TRYCXX( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer1) );
	TRYCXX( PetscViewerCreate(PETSC_COMM_WORLD, &mviewer2) );
	TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_WORLD ,filename1.c_str(), FILE_MODE_READ, &mviewer1) );
	TRYCXX( PetscViewerBinaryOpen(PETSC_COMM_WORLD ,filename2.c_str(), FILE_MODE_READ, &mviewer2) );
	
	/* load vector from viewer */
	TRYCXX( VecLoad(myvector1_Vec, mviewer1) );
	TRYCXX( VecLoad(myvector2_Vec, mviewer2) );

	/* destroy the viewer */
	TRYCXX( PetscViewerDestroy(&mviewer1) );
	TRYCXX( PetscViewerDestroy(&mviewer2) );

	mytimer.stop();
	coutMaster << "- time load      : " << mytimer.get_value_last() << " s" << std::endl;

	int size1, size2;
	TRYCXX( VecGetSize(myvector1_Vec, &size1) );
	TRYCXX( VecGetSize(myvector2_Vec, &size2) );

	coutMaster << "- dimension of v1: " << size1 << std::endl;
	coutMaster << "- dimension of v2: " << size2 << std::endl;


	mytimer.start();
	double dot_result;
	for(int i=0;i<n;i++){
		TRYCXX( VecDot( myvector1_Vec, myvector2_Vec, &dot_result) );
		TRYCXX( VecDot( myvector1_Vec, myvector2_Vec, &dot_result) );
		TRYCXX( VecDot( myvector1_Vec, myvector2_Vec, &dot_result) );
	}
	mytimer.stop();
	coutMaster << "- result        : " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << dot_result << std::endl;
	coutMaster << "- time total     : " << mytimer.get_value_last() << " s" << std::endl;
	coutMaster << "- time average   : " << mytimer.get_value_last()/(double)n << " s" << std::endl;
	coutMaster << std::endl;

	
	PetscScalar Mdots_val[3];
	Vec Mdots_vec[3];

	Mdots_vec[0] = myvector2_Vec;
	Mdots_vec[1] = myvector2_Vec;
	Mdots_vec[2] = myvector2_Vec;

	coutMaster << "--- start to compute ---" << std::endl;
	
	/* we will measure the time of operation */
	mytimer.start();
	for(int i=0;i<n;i++){
		TRYCXX( VecMDot( myvector1_Vec, 3, Mdots_vec, Mdots_val) );
	}
	mytimer.stop();

	coutMaster << "---- computed ----" << std::endl;

	coutMaster << "- result1        : " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << Mdots_val[0] << std::endl;
	coutMaster << "- result2        : " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << Mdots_val[1] << std::endl;
	coutMaster << "- result3        : " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << Mdots_val[2] << std::endl;
	coutMaster << "- time total     : " << mytimer.get_value_last() << " s" << std::endl;
	coutMaster << "- time average   : " << mytimer.get_value_last()/(double)n << " s" << std::endl;
	coutMaster << std::endl;

	coutMaster << std::endl;

	Finalize();

	return 0;
}


