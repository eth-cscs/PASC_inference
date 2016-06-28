/** @file test_mat_blocklaplace_free_to_dense.cu
 *  @brief print block laplace free-matrix as dense
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"

#include "matrix/blocklaplace.h"
#include "data/kmeansdata.h"

#ifndef USE_PETSCVECTOR
 #error 'This example is for PETSCVECTOR'
#endif

#ifdef USE_CUDA
    #include <../src/vec/vec/impls/seq/seqcuda/cudavecimpl.h>

#endif

typedef petscvector::PetscVector PetscVector;
 
using namespace pascinference;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_K", boost::program_options::value<int>(), "number of clusters")
		("test_T", boost::program_options::value<int>(), "dimension of the problem")
		("test_n", boost::program_options::value<int>(), "number of tests");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	/* load console arguments */
	int T, K, n;
	
	consoleArg.set_option_value("test_T", &T, 10);
	consoleArg.set_option_value("test_K", &K, 2);
	consoleArg.set_option_value("test_n", &n, 10);

	coutMaster << " T        = " << std::setw(7) << T << " (length of time-series)" << std::endl;
	coutMaster << " K        = " << std::setw(7) << K << " (number of clusters)" << std::endl;
	coutMaster << " n        = " << std::setw(7) << n << " (number of tests)" << std::endl;

	/* start logging */
	std::ostringstream oss_name_of_file_log;
	oss_name_of_file_log << "results/test_mat_blocklaplace_free_to_dense_log_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_name_of_file_log.str());

	/* say hello */
	coutMaster << "- start program" << std::endl;

	/* create base vector */
	/* try to make a global vector of length T and then get size of local portion */
	Vec x_mpi_Vec, x_seq_Vec, y_mpi_Vec, y_seq_Vec;
	int Tlocal, Tbegin, Tend;
	
	Vec layout;
	TRY( VecCreate(PETSC_COMM_WORLD,&layout) );
	TRY( VecSetSizes(layout,PETSC_DECIDE,T) );
	TRY( VecSetFromOptions(layout) );
	/* get the ownership range - now I know how much I will calculate from the time-series */
	TRY( VecGetLocalSize(layout,&Tlocal) );
	TRY( VecGetOwnershipRange(layout,&Tbegin,&Tend) );
	TRY( VecDestroy(&layout) ); /* destroy testing vector - it is useless now */

	#ifndef USE_GPU
		/* kmeans data will help us with distribution of SEQ to MPI */
		KmeansData<PetscVector> mydata(T,K);
		x_mpi_Vec = mydata.get_datavector()->get_vector();

		TRY( VecDuplicate(x_mpi_Vec, &y_mpi_Vec) );

		TRY( VecCreateSeq(PETSC_COMM_SELF,K*T, &x_seq_Vec) );
		TRY( VecDuplicate(x_seq_Vec, &y_seq_Vec) );

	#else
		TRY( VecCreate(PETSC_COMM_WORLD,&x_Vec) );
		TRY( VecSetType(x_Vec, VECMPICUDA) );
		TRY( VecSetSizes(x_Vec,K*Tlocal,K*T) );
		TRY( VecSetFromOptions(x_Vec) );

		TRY( VecCreateMPI(PETSC_COMM_WORLD,K*Tlocal,K*T, &y_VecCPU) ); /* for printing I need vector on CPU */
 
//		TRY( VecCreateMPICUDA(PETSC_COMM_WORLD,K*Tlocal,K*T, &x_Vec) );
	#endif

	/* create general vectors */
	GeneralVector<PetscVector> x_mpi(x_mpi_Vec);
	GeneralVector<PetscVector> y_mpi(y_mpi_Vec);
	BlockLaplaceMatrix<PetscVector> A_mpi(x_mpi, K);

	Timer timer1;
	timer1.restart();

	if(GlobalManager.get_rank == 0){
		GeneralVector<PetscVector> x_seq(x_seq_Vec);
		GeneralVector<PetscVector> y_seq(y_seq_Vec);
		BlockLaplaceMatrix<PetscVector> A_seq(x_seq, K);

		Timer timer2;
		timer2.restart();
	}
	
	/* prepare random generator */
	PetscRandom rnd;
	TRY( PetscRandomCreate(PETSC_COMM_WORLD,&rnd) );
	TRY( PetscRandomSetType(rnd,PETSCRAND) );
	TRY( PetscRandomSetFromOptions(rnd) );

	TRY( PetscRandomSetSeed(rnd,13) );

	coutMaster << "running main fun" << std::endl;
	int ni;
	double mynorm_seq;
	double mynorm_mpi;
	
	for(ni = 0; ni < n; ni++){
		if(GlobalManager.get_rank()==0){
			TRY( VecSetRandom(x_seq_Vec, rnd) );
		}

//		TRY( VecAssemblyBegin(x_Vec) );
//		TRY( VecAssemblyEnd(x_Vec) );

		timer1.start();
			y_seq = A_seq*x_seq;
		timer1.stop();

		timer2.start();
//			y_mpi = A_mpi*x_mpi;
		timer2.stop();

		TRY( VecDot(y_seq_Vec,y_seq_Vec,&mynorm_seq) );
//		TRY( VecDot(y_mpi_Vec,y_mpi_Vec,&mynorm_mpi) );

		coutMaster << "n=" << ni << ": ";
		coutMaster << "mpi=" << mynorm_mpi << ", ";
		coutMaster << "seq=" << mynorm_seq;
		coutMaster << std::endl;
		
	}

	/* destroy the random generator */
	TRY( PetscRandomDestroy(&rnd) );


	coutMaster << "T=" << std::setw(9) << T << ", K="<< std::setw(4) << K << ", time_mpi=" << std::setw(10) << timer2.get_value_sum() << ", time_seq=";
	coutAll << std::setw(10) << timer1.get_value_sum();
	coutAll.synchronize();
	coutMaster << std::endl;

	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize();

	return 0;
}

