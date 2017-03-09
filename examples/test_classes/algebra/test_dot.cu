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

//extern int pascinference::DEBUG_MODE;

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
	
	/* we will measure the time of operation */
	mytimer.start();
	double dot_result;
	for(int i=0;i<n;i++){

		TRYCXX( VecDot(myvector1_Vec, myvector2_Vec, &dot_result) );

	}
	mytimer.stop();

	coutMaster << "- result         : " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << dot_result << std::endl;
	coutMaster << "- time total     : " << mytimer.get_value_last() << " s" << std::endl;
	coutMaster << "- time average   : " << mytimer.get_value_last()/(double)n << " s" << std::endl;
	coutMaster << std::endl;

	coutMaster << "TEST of the wrapper" << std::endl;
	GeneralVector<PetscVector> myvector1(myvector1_Vec);
	GeneralVector<PetscVector> myvector2(myvector2_Vec);

	mytimer.start();
	double dot_result2;
	for(int i=0;i<n;i++){
		dot_result2 = dot(myvector1,myvector2);
	}
	mytimer.stop();

	coutMaster << "- result         : " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << dot_result2 << std::endl;
	coutMaster << "- time total     : " << mytimer.get_value_last() << " s" << std::endl;
	coutMaster << "- time average   : " << mytimer.get_value_last()/(double)n << " s" << std::endl;
	coutMaster << std::endl;


	coutMaster << "TEST of the wrapper2" << std::endl;

	typedef GeneralVector<PetscVector> (&pVector);
	pVector b1 = *(&myvector1);
	pVector b2 = *(&myvector2);

	mytimer.start();
	double dot_result3;
	for(int i=0;i<n;i++){
		dot_result3 = dot(b1,b2);
	}
	mytimer.stop();

	coutMaster << "- result         : " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << dot_result3 << std::endl;
	coutMaster << "- time total     : " << mytimer.get_value_last() << " s" << std::endl;
	coutMaster << "- time average   : " << mytimer.get_value_last()/(double)n << " s" << std::endl;
	coutMaster << std::endl;


	coutMaster << "TEST of the Vec = Vec" << std::endl;

	Vec myvector3_Vec;
	Vec myvector4_Vec;
	
	myvector3_Vec = b1.get_vector();
	myvector4_Vec = b2.get_vector();

	Vec myvector5_Vec;
	Vec myvector6_Vec;
	TRYCXX( VecDuplicate(myvector3_Vec, &myvector5_Vec) );
	TRYCXX( VecDuplicate(myvector4_Vec, &myvector6_Vec) );
	TRYCXX( VecCopy(myvector3_Vec, myvector5_Vec) );
	TRYCXX( VecCopy(myvector4_Vec, myvector6_Vec) );

	double dot_result4;
	mytimer.start();
	for(int i=0;i<n;i++){
		TRYCXX( VecDot(myvector5_Vec, myvector6_Vec, &dot_result4) );
	}
	mytimer.stop();

	coutMaster << "- result         : " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << dot_result4 << std::endl;
	coutMaster << "- time total     : " << mytimer.get_value_last() << " s" << std::endl;
	coutMaster << "- time average   : " << mytimer.get_value_last()/(double)n << " s" << std::endl;
	coutMaster << std::endl;


	coutMaster << "TEST of the IS" << std::endl;
	
	int localsize1;
	int localsize2;
	
	TRYCXX( VecGetLocalSize(b1.get_vector(), &localsize1) );
	TRYCXX( VecGetLocalSize(b2.get_vector(), &localsize2) );
	
	int start1, end1;
	int start2, end2;
	TRYCXX( VecGetOwnershipRange(b1.get_vector(), &start1, &end1) );
	TRYCXX( VecGetOwnershipRange(b2.get_vector(), &start2, &end2) );
	
	coutMaster << "- create stride" << std::endl;
	IS myis1;
	IS myis2;
	TRYCXX( ISCreateStride(PETSC_COMM_WORLD, localsize1, end2-1, -1, &myis1) );
	TRYCXX( ISCreateStride(PETSC_COMM_WORLD, localsize2, end2-1, -1, &myis2) );

	Vec new1_Vec;
	Vec new2_Vec;
	Vec myvector73;
	TRYCXX( VecDuplicate(b2.get_vector(),&myvector73));

	coutMaster << "- get subvector" << std::endl;
	TRYCXX( VecGetSubVector(b1.get_vector(), myis1, &new1_Vec) );
	TRYCXX( VecGetSubVector(myvector73, myis2, &new2_Vec) );

	coutMaster << "- copy" << std::endl;
	TRYCXX( VecCopy(new1_Vec, new2_Vec) );
	
	coutMaster << "- restore" << std::endl;
	TRYCXX( VecRestoreSubVector(b1.get_vector(), myis1, &new1_Vec) );
	TRYCXX( VecRestoreSubVector(myvector73, myis2, &new2_Vec) );
	
	const char *mytype1;
	const char *mytype2;
	coutMaster << "- get type" << std::endl;
	TRYCXX( VecGetType(b2.get_vector(), &mytype1) );
	TRYCXX( VecGetType(myvector73, &mytype2) );
	
	coutMaster << " - type1: " << mytype1 << std::endl;
	coutMaster << " - type2: " << mytype2 << std::endl;
	
	coutMaster << "- compute" << std::endl;
	double dot_result6;
	mytimer.start();
	for(int i=0;i<n;i++){
		TRYCXX( VecDot(b2.get_vector(), myvector73, &dot_result6) );
	}
	mytimer.stop();


	coutMaster << "- result         : " << std::setprecision(std::numeric_limits<long double>::digits10 + 1) << dot_result6 << std::endl;
	coutMaster << "- time total     : " << mytimer.get_value_last() << " s" << std::endl;
	coutMaster << "- time average   : " << mytimer.get_value_last()/(double)n << " s" << std::endl;

	coutMaster << std::endl;

	Finalize();

	return 0;
}


