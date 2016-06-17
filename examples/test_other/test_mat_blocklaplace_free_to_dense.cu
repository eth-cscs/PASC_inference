/** @file test_mat_blocklaplace_free_to_dense.cu
 *  @brief print block laplace free-matrix as dense
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include "matrix/blocklaplace.h"

#ifndef USE_PETSCVECTOR
 #error 'This example is for PETSCVECTOR'
#endif

#ifdef USE_CUDA
    #include <../src/vec/vec/impls/seq/seqcuda/cudavecimpl.h>

#endif

typedef petscvector::PetscVector PetscVector;
 
using namespace pascinference;

#ifdef USE_CUDA
	__global__ void write_to_array(double* arr, int idx, double value){
		arr[idx] = value;
	}
#endif

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_K", boost::program_options::value<int>(), "number of clusters")
		("test_T", boost::program_options::value<int>(), "dimension of the problem")
		("test_view_matrix", boost::program_options::value<bool>(), "print dense matrix or not");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	/* load console arguments */
	int T, K;
	bool test_view_matrix;
	
	consoleArg.set_option_value("test_T", &T, 10);
	consoleArg.set_option_value("test_K", &K, 2);
	consoleArg.set_option_value("test_view_matrix", &test_view_matrix, true);

	coutMaster << " T        = " << std::setw(7) << T << " (length of time-series)" << std::endl;
	coutMaster << " K        = " << std::setw(7) << K << " (number of clusters)" << std::endl;

	/* start logging */
	std::ostringstream oss_name_of_file_log;
	oss_name_of_file_log << "results/test_mat_blocklaplace_free_to_dense_log_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_name_of_file_log.str());

	/* say hello */
	coutMaster << "- start program" << std::endl;

	/* create base vector */
	/* try to make a global vector of length T and then get size of local portion */
	Vec x_Vec, y_VecCPU;
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
		TRY( VecCreateMPI(PETSC_COMM_WORLD,K*Tlocal,K*T, &x_Vec) );
	#else
		TRY( VecCreate(PETSC_COMM_WORLD,&x_Vec) );
		TRY( VecSetType(x_Vec, VECMPICUDA) );
		TRY( VecSetSizes(x_Vec,K*Tlocal,K*T) );
		TRY( VecSetFromOptions(x_Vec) );

		TRY( VecCreateMPI(PETSC_COMM_WORLD,K*Tlocal,K*T, &y_VecCPU) ); /* for printing I need vector on CPU */

//		TRY( VecCreateMPICUDA(PETSC_COMM_WORLD,K*Tlocal,K*T, &x_Vec) );
	#endif

	GeneralVector<PetscVector> x(x_Vec);
	GeneralVector<PetscVector> y(x); /* result, i.e. y = A*x */

	/* create matrix */
	BlockLaplaceMatrix A(x, K);
	A.print(coutMaster,coutAll);
	coutAll.synchronize();

	Timer timer1;
	timer1.restart();
	
	/* print full matrix as a multiplication A*e_n */
	double *values;
	double *gamma_arr;
	int k,t,tlocal,idx;
	
	for(int row_idx = 0; row_idx < T*K; row_idx++){
		k = floor(row_idx/(double)T);
		t = row_idx - k*T;
		
		/* prepare e_k in layout */
		TRY( VecSet(x_Vec,0) );
		TRY( VecAssemblyBegin(x_Vec) );
		TRY( VecAssemblyEnd(x_Vec) );

		#ifndef USE_GPU
			TRY( VecGetArray(x_Vec,&gamma_arr) );

			if(t >= Tbegin && t < Tend){
				tlocal = t - Tbegin;
				gamma_arr[tlocal + k*Tlocal] = 1;
			}

			TRY( VecRestoreArray(x_Vec,&gamma_arr) );
		#else
			TRY( VecCUDAGetArrayReadWrite(x_Vec,&gamma_arr) );

			if(t >= Tbegin && t < Tend){
				write_to_array<<<1, 1>>>(gamma_arr, tlocal + k*Tlocal, 1.0);
				gpuErrchk( cudaDeviceSynchronize() );
			}

			TRY( VecCUDARestoreArrayReadWrite(x_Vec,&gamma_arr) );
		#endif

		TRY( VecAssemblyBegin(x_Vec) );
		TRY( VecAssemblyEnd(x_Vec) );

	timer1.start();
		y = A*x;
	timer1.stop();
		
		if(test_view_matrix){
			/* get array of values */
			#ifndef USE_GPU
				TRY( VecGetArray(y.get_vector(), &values) );
			#else
				/* from GPU to CPU and then get array */
				TRY( VecCopy(y.get_vector(),y_VecCPU) );
				TRY( VecGetArray(y_Vec, &values) );
			#endif
		
			/* print row */
			TRY( PetscPrintf(PETSC_COMM_WORLD, "%*d: ", 3, row_idx) );

			/* interpret results from new layout to standart one */
			for(k=0;k<K;k++){
				for(t=0;t<Tlocal;t++){
					idx = t+k*Tlocal;
					if(abs(values[idx]) > 0.0001){
						TRY( PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%*.*f,",6, 2, values[idx]) );
					} else {
						TRY( PetscSynchronizedPrintf(PETSC_COMM_WORLD, "      ,") );
					}
				}
				TRY( PetscSynchronizedFlush(PETSC_COMM_WORLD, NULL) );
			}
			TRY( PetscPrintf(PETSC_COMM_WORLD, "\n") );

			/* restore array with values */
			#ifndef USE_GPU
				TRY( VecRestoreArray(y.get_vector(), &values) );
			#else
				TRY( VecRestoreArray(y_Vec, &values) );
			#endif

		}
	}


	coutMaster << "T = " << std::setw(9) << T << ", K = "<< std::setw(4) << K << ", time = " << std::setw(10) << timer1.get_value_sum()/(double)(K) << std::endl;

	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize();

	return 0;
}

