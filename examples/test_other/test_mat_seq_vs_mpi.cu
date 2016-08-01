/** @file test_mat_blocklaplace_free_to_dense.cu
 *  @brief print block laplace free-matrix as dense
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"


#include "matrix/blocklaplacefree.h"
#include "matrix/blocklaplacesparse.h"

#include "matrix/blockgraphfree.h"
#include "matrix/blockgraphsparse.h"

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
		("test_matrix_type", boost::program_options::value<int>(), "which matrix type to test [0=BLOCKGRAPHFREE/1=BLOCKGRAPHSPARSE/2=BLOCKLAPLACEFREE/3=BLOCKLAPLACESPARSE]")
		("test_K", boost::program_options::value<int>(), "number of clusters")
		("test_T", boost::program_options::value<int>(), "dimension of the problem")
		("test_n", boost::program_options::value<int>(), "number of tests")
		("test_alpha", boost::program_options::value<double>(), "coefficient of the matrix [double]")
		("test_graph_filename", boost::program_options::value< std::string >(), "name of file with coordinates [string]")
		("test_graph_coeff", boost::program_options::value<double>(), "threshold of the graph [double]")
		("test_print_vectors", boost::program_options::value<bool>(), "print content of vectors or not [bool]");

	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	/* load console arguments */
	int T, K, R, n, matrix_type;
	std::string graph_filename;
	double alpha, graph_coeff;
	bool print_vectors;
	
	consoleArg.set_option_value("test_matrix_type", &matrix_type, 0);
	consoleArg.set_option_value("test_T", &T, 10);
	consoleArg.set_option_value("test_K", &K, 2);
	consoleArg.set_option_value("test_n", &n, 1);
	consoleArg.set_option_value("test_alpha", &alpha, 1.0);
	consoleArg.set_option_value("test_graph_filename", &graph_filename, "data/graph_twonodes.bin");
	consoleArg.set_option_value("test_graph_coeff", &graph_coeff, 1.1);
	consoleArg.set_option_value("test_print_vectors", &print_vectors, true);

	coutMaster << " test_matrix_type          = " << std::setw(20) << matrix_type << " (which matrix type to test)" << std::endl;
	coutMaster << " test_T                    = " << std::setw(20) << T << " (length of time-series)" << std::endl;
	coutMaster << " test_K                    = " << std::setw(20) << K << " (number of clusters)" << std::endl;
	coutMaster << " test_n                    = " << std::setw(20) << n << " (number of tests)" << std::endl;
	coutMaster << " test_alpha                = " << std::setw(20) << alpha << " (coeficient of the matrix)" << std::endl;
	coutMaster << " test_graph_filename       = " << std::setw(20) << graph_filename << " (name of file with coordinates)" << std::endl;
	coutMaster << " test_graph_coeff          = " << std::setw(20) << graph_coeff << " (threshold of the graph)" << std::endl;
	coutMaster << " test_print_vectors        = " << std::setw(20) << print_vectors << " (print content of vectors or not)" << std::endl;

	/* start logging */
	std::ostringstream oss_name_of_file_log;
	oss_name_of_file_log << "results/test_mat_seq_vs_mpi_log_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_name_of_file_log.str());

	/* say hello */
	coutMaster << "- start program" << std::endl;

	/* create graph */
	BGMGraph graph(graph_filename);
	graph.process(graph_coeff);

	if(matrix_type==0 || matrix_type==1){
		/* graph matrix */
		R = graph.get_n();
	}
	if(matrix_type==2 || matrix_type==3){
		/* laplace matrix */
		R = 1;
	}

	coutAll << "T=" << T << std::endl;
	coutAll.synchronize();

	/* create base vector */
	Vec x_mpi_Vec, x_seq_Vec, y_mpi_Vec, y_seq_Vec;
	int Tlocal, Tbegin, Tend;
	Vec layout;
	TRY( VecCreate(PETSC_COMM_WORLD,&layout) );
	TRY( VecSetSizes(layout, PETSC_DECIDE, T ));
	TRY( VecSetFromOptions(layout) );

	/* get the ownership range - now I know how much I will calculate from the time-series */
	TRY( VecGetLocalSize(layout,&Tlocal) );
	TRY( VecGetOwnershipRange(layout,&Tbegin,&Tend) );
	TRY( VecDestroy(&layout) ); /* destroy testing vector - it is useless now */
		
	/* MPI */
	TRY( VecCreate(PETSC_COMM_WORLD,&x_mpi_Vec) );
	TRY( VecSetSizes(x_mpi_Vec,R*K*Tlocal,R*K*T) );
	TRY( VecSetFromOptions(x_mpi_Vec) );
	TRY( VecDuplicate(x_mpi_Vec, &y_mpi_Vec) );

	GeneralVector<PetscVector> x_mpi(x_mpi_Vec);
	GeneralVector<PetscVector> y_mpi(y_mpi_Vec);

	/* SEQ */
	TRY( VecCreateSeq(PETSC_COMM_SELF,K*T*R, &x_seq_Vec) );
	TRY( VecDuplicate(x_seq_Vec, &y_seq_Vec) );
	GeneralVector<PetscVector> x_seq(x_seq_Vec);
	GeneralVector<PetscVector> y_seq(y_seq_Vec);

	/* difference SEQ - MPI */
	Vec diff_seq_mpi_Vec;
	TRY( VecDuplicate(x_seq_Vec, &diff_seq_mpi_Vec) );
	
	/* create matrix */
	GeneralMatrix<PetscVector> *A_mpi, *A_seq;
	if(matrix_type==0){
		A_mpi = new BlockGraphFreeMatrix<PetscVector>(x_mpi, graph, K, alpha);
		A_seq = new BlockGraphFreeMatrix<PetscVector>(x_seq, graph, K, alpha);
	}
	if(matrix_type==1){
		A_mpi = new BlockGraphSparseMatrix<PetscVector>(x_mpi, graph, K, alpha);
		A_seq = new BlockGraphSparseMatrix<PetscVector>(x_seq, graph, K, alpha);
	}
	if(matrix_type==2){
		A_mpi = new BlockLaplaceFreeMatrix<PetscVector>(x_mpi, K, alpha);
		A_seq = new BlockLaplaceFreeMatrix<PetscVector>(x_seq, K, alpha);
	}
	if(matrix_type==3){
		A_mpi = new BlockLaplaceSparseMatrix<PetscVector>(x_mpi, K, alpha);
		A_seq = new BlockLaplaceSparseMatrix<PetscVector>(x_seq, K, alpha);
	}

	/* timers */
	Timer timer_mpi;
	timer_mpi.restart();

	Timer timer_seq;
	timer_seq.restart();
	
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
	double mynorm_err;

	coutMaster << "running main fun" << std::endl;
	coutMaster << " Tlocal   = " << std::setw(7) << Tlocal << std::endl;
	coutMaster << " Tbegin   = " << std::setw(7) << Tbegin << std::endl;
	coutMaster << " Tend     = " << std::setw(7) << Tend << std::endl;

	for(ni = 0; ni < n; ni++){
		/* set random values to values_Vec */
		TRY( VecSetRandom(x_mpi_Vec, rnd) );
		TRY( VecAssemblyBegin(x_mpi_Vec) );
		TRY( VecAssemblyEnd(x_mpi_Vec) );

		/* scatter global values to local */
		/* prepare IS with my indexes in provided vector */
		IS x_sub_IS;
		TRY( ISCreateStride(PETSC_COMM_SELF, T*K*R, 0, 1, &(x_sub_IS)) );
		/* distribte data to all processors */
		Vec x_sub_Vec;
		TRY( VecGetSubVector(x_mpi_Vec, x_sub_IS, &x_sub_Vec) );
		/* copy values */
		TRY( VecCopy(x_sub_Vec, x_seq_Vec) );
		/* restore vector */
		TRY( VecRestoreSubVector(x_mpi_Vec, x_sub_IS, &x_sub_Vec) );

		/* perform muliplication */
		timer_seq.start();
			y_seq = (*A_seq)*x_seq;
		timer_seq.stop();

		timer_mpi.start();
			y_mpi = (*A_mpi)*x_mpi;
		timer_mpi.stop();

		/* scatter y_mpi to one proc to compute diference */
		TRY( VecGetSubVector(y_mpi_Vec, x_sub_IS, &x_sub_Vec) );
		TRY( VecCopy(x_sub_Vec, diff_seq_mpi_Vec) );
		TRY( VecRestoreSubVector(y_mpi_Vec, x_sub_IS, &x_sub_Vec) );
		/* compute difference */
		TRY( VecAXPY(diff_seq_mpi_Vec,-1.0,y_seq_Vec) );
		TRY( ISDestroy(&x_sub_IS) );

		/* get arrays for print */
		TRY( VecGetArray(x_mpi_Vec, &x_mpi_arr) );
		TRY( VecGetArray(x_seq_Vec, &x_seq_arr) );
		TRY( VecGetArray(y_mpi_Vec, &y_mpi_arr) );
		TRY( VecGetArray(y_seq_Vec, &y_seq_arr) );

		if(print_vectors){
			coutMaster << "y_seq (master): " << std::endl;
			for(int k=0; k<K; k++){
				for(int r=0;r<R;r++){
					coutMaster << " k=" << k << ",r=" << r << ": x = [";

					for(int t=0;t<T;t++){
						coutMaster << x_seq_arr[t*K*R+r*K+k];
						if(t < T-1){
							coutMaster << ",";
						}
					}
				
					coutMaster << "], y = [";
					for(int t=0;t<T;t++){
						coutMaster << y_seq_arr[t*K*R+r*K+k];
						if(t < T-1){
							coutMaster << ",";
						}
					}
					coutMaster << "]" << std::endl;
				}
			}
		
			coutMaster << "y_mpi: " << std::endl;
			for(int k=0; k<K; k++){
				for(int r=0;r<R;r++){
					TRY( PetscPrintf(PETSC_COMM_WORLD, " k=%d,r=%d: x = [", k, r) );
					for(int t=0;t<Tlocal;t++){
						TRY( PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%f",x_mpi_arr[t*K*R+r*K+k]) );
						if(t < Tlocal-1){
							TRY( PetscSynchronizedPrintf(PETSC_COMM_WORLD, ",") );
						}
					}
					TRY( PetscSynchronizedFlush(PETSC_COMM_WORLD, NULL) );

					TRY( PetscPrintf(PETSC_COMM_WORLD, "], y = [") );	
					for(int t=0;t<Tlocal;t++){
						TRY( PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%f",y_mpi_arr[t*K*R+r*K+k]) );
						if(t < Tlocal-1){
							TRY( PetscSynchronizedPrintf(PETSC_COMM_WORLD, ",") );
						}
					}
					TRY( PetscSynchronizedFlush(PETSC_COMM_WORLD, NULL) );
				
					TRY( PetscPrintf(PETSC_COMM_WORLD, "]\n") );	

				}
			}
		}

		/* restore arrays for print */
		TRY( VecRestoreArray(y_mpi_Vec, &y_mpi_arr) );
		TRY( VecRestoreArray(y_seq_Vec, &y_seq_arr) );
		TRY( VecRestoreArray(x_mpi_Vec, &x_mpi_arr) );
		TRY( VecRestoreArray(x_seq_Vec, &x_seq_arr) );
	
		TRY( VecDot(y_seq_Vec,y_seq_Vec,&mynorm_seq) );
		TRY( VecDot(y_mpi_Vec,y_mpi_Vec,&mynorm_mpi) );
		TRY( VecDot(diff_seq_mpi_Vec,diff_seq_mpi_Vec,&mynorm_err) );

		std::streamsize ss = std::cout.precision();
		coutMaster << std::setprecision(17);
		coutMaster << "n=" << ni << ": ";
		coutMaster << "err=" << std::setw(10) << mynorm_err << ": ";
		coutMaster << "mpi=" << mynorm_mpi << ", ";
		coutMaster << "seq=" << mynorm_seq;
		coutMaster << std::endl;
		coutMaster << std::setprecision(ss);

	}

	coutMaster << "T=" << std::setw(9) << T << ", K="<< std::setw(4) << K << ", R="<< std::setw(4) << R << ", time_mpi=" << std::setw(10) << timer_mpi.get_value_sum() << ", time_seq=" << std::setw(10) << timer_seq.get_value_sum() << std::endl;

	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize();

	return 0;
}

