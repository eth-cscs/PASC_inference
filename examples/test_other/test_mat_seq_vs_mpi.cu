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
	consoleArg.set_option_value("test_n", &n, 1);

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
	Vec x_mpi_Vec, x_seq_Vec, y_mpi_Vec, y_seq_Vec, values_Vec;

	#ifndef USE_GPU
		/* kmeans data will help us with distribution of SEQ to MPI */
		KmeansData<PetscVector> data_values(T,K);
		values_Vec = data_values.get_datavector()->get_vector();
		TRY( VecDuplicate(values_Vec, &x_mpi_Vec) );
		TRY( VecDuplicate(values_Vec, &y_mpi_Vec) );

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

	GeneralVector<PetscVector> x_seq(x_seq_Vec);
	GeneralVector<PetscVector> y_seq(y_seq_Vec);
	BlockLaplaceMatrix<PetscVector> A_seq(x_seq, K);

	Timer timer2;
	timer2.restart();
	
	/* prepare random generator */
	PetscRandom rnd;
	TRY( PetscRandomCreate(PETSC_COMM_WORLD,&rnd) );
	TRY( PetscRandomSetType(rnd,PETSCRAND) );
	TRY( PetscRandomSetFromOptions(rnd) );

	TRY( PetscRandomSetSeed(rnd,13) );

	double *x_seq_arr, *x_mpi_arr;
	double *y_seq_arr, *y_mpi_arr;

	int ni;
	double mynorm_seq;
	double mynorm_mpi;

	KmeansData<PetscVector> data_mpi(NULL, &x_mpi, NULL, T);
	KmeansData<PetscVector> data_seq(NULL, &x_seq, NULL, T);

	coutMaster << "running main fun" << std::endl;

	GeneralVector<PetscVector> values(values_Vec);

	int Tlocal = data_values.get_Tlocal();
	int Tbegin = data_values.get_Tbegin();
	int Tend = data_values.get_Tend();

	coutMaster << " Tlocal   = " << std::setw(7) << Tlocal << std::endl;
	coutMaster << " Tbegin   = " << std::setw(7) << Tbegin << std::endl;
	coutMaster << " Tend     = " << std::setw(7) << Tend << std::endl;

	for(ni = 0; ni < n; ni++){
		TRY( VecSetRandom(values_Vec, rnd) );
		TRY( VecAssemblyBegin(values_Vec) );
		TRY( VecAssemblyEnd(values_Vec) );

		/* scatter global values to local */
		data_mpi.load_gammavector(values);
		data_seq.load_gammavector(values);
		
		timer1.start();
			y_seq = A_seq*x_seq;
		timer1.stop();

		timer2.start();
			y_mpi = A_mpi*x_mpi;
		timer2.stop();

		/* get arrays for print */
		TRY( VecGetArray(x_mpi_Vec, &x_mpi_arr) );
		TRY( VecGetArray(x_seq_Vec, &x_seq_arr) );
		TRY( VecGetArray(y_mpi_Vec, &y_mpi_arr) );
		TRY( VecGetArray(y_seq_Vec, &y_seq_arr) );

		coutMaster << "y_seq (master): " << std::endl;
		for(int t=0; t<T; t++){
			coutMaster << " t = " << t << ": x = [";
			for(int k=0;k<K;k++){
				coutMaster << x_seq_arr[k*T+t];
				if(k < K-1){
					coutMaster << ",";
				}
			}
			coutMaster << "], y = [";
			for(int k=0;k<K;k++){
				coutMaster << y_seq_arr[k*T+t];
				if(k < K-1){
					coutMaster << ",";
				}
			}
			coutMaster << "]" << std::endl;
			coutMaster.synchronize();
		}
		
		coutMaster << "y_mpi: " << std::endl;
		for(int t=0; t<Tlocal; t++){
			coutAll << " t = " << Tbegin+t << ": x = [";
			for(int k=0;k<K;k++){
				coutAll << x_mpi_arr[k*Tlocal+t];
				if(k < K-1){
					coutAll << ",";
				}
			}
			coutAll << "], y = [";
			for(int k=0;k<K;k++){
				coutAll << y_mpi_arr[k*Tlocal+t];
				if(k < K-1){
					coutAll << ",";
				}
			}
			coutAll << "]" << std::endl;
		}
		coutAll.synchronize();


		/* restore arrays for print */
		TRY( VecRestoreArray(y_mpi_Vec, &y_mpi_arr) );
		TRY( VecRestoreArray(y_seq_Vec, &y_seq_arr) );
		TRY( VecRestoreArray(x_mpi_Vec, &x_mpi_arr) );
		TRY( VecRestoreArray(x_seq_Vec, &x_seq_arr) );
	
		TRY( VecDot(y_seq_Vec,y_seq_Vec,&mynorm_seq) );
		TRY( VecDot(y_mpi_Vec,y_mpi_Vec,&mynorm_mpi) );

		coutMaster << "n=" << ni << ": ";
		coutMaster << "mpi=" << mynorm_mpi << ", ";
		coutMaster << "seq=" << mynorm_seq;
		coutMaster << std::endl;

	}

	coutMaster << "T=" << std::setw(9) << T << ", K="<< std::setw(4) << K << ", time_mpi=" << std::setw(10) << timer2.get_value_sum() << ", time_seq=" << std::setw(10) << timer1.get_value_sum() << std::endl;

	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize();

	return 0;
}

