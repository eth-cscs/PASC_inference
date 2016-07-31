/** @file test_mat_blockgraph_free_to_dense.cu
 *  @brief print free matrix as dense
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

typedef petscvector::PetscVector PetscVector;
 
using namespace pascinference;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_matrix_type", boost::program_options::value<int>(), "which matrix type to test [0=BLOCKGRAPHFREE/1=BLOCKGRAPHSPARSE/2=BLOCKLAPLACEFREE/3=BLOCKLAPLACESPARSE]")
		("test_T", boost::program_options::value<int>(), "dimension of the problem [int]")
		("test_K", boost::program_options::value<int>(), "number of clusters [int]")
		("test_alpha", boost::program_options::value<double>(), "coefficient of the matrix [double]")
		("test_graph_filename", boost::program_options::value< std::string >(), "name of file with coordinates [string]")
		("test_graph_coeff", boost::program_options::value<double>(), "threshold of the graph [double]")
		("test_view_graph", boost::program_options::value<bool>(), "print content of graph or not [bool]")
		("test_view_matrix", boost::program_options::value<bool>(), "print dense matrix or not [bool]")
		("test_view_matrix_unsorted", boost::program_options::value<bool>(), "print dense matrix in original unsorted form or not [bool]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	/* load console arguments */
	int T, K, R, matrix_type;
	std::string graph_filename;
	bool view_matrix, view_graph, view_matrix_unsorted;
	double alpha, graph_coeff;
	
	consoleArg.set_option_value("test_matrix_type", &matrix_type, 0);
	consoleArg.set_option_value("test_T", &T, 1);
	consoleArg.set_option_value("test_K", &K, 1);
	consoleArg.set_option_value("test_alpha", &alpha, 1.0);
	consoleArg.set_option_value("test_graph_filename", &graph_filename, "data/graph_twonodes.bin");
	consoleArg.set_option_value("test_graph_coeff", &graph_coeff, 1.1);
	consoleArg.set_option_value("test_view_matrix", &view_matrix, true);
	consoleArg.set_option_value("test_view_matrix_unsorted", &view_matrix_unsorted, false);
	consoleArg.set_option_value("test_view_graph", &view_graph, false);

	/* print settings */
	coutMaster << " test_matrix_type          = " << std::setw(20) << matrix_type << " (which matrix type to test)" << std::endl;
	coutMaster << " test_T                    = " << std::setw(20) << T << " (length of time-series)" << std::endl;
	coutMaster << " test_K                    = " << std::setw(20) << K << " (number of clusters)" << std::endl;
	coutMaster << " test_alpha                = " << std::setw(20) << alpha << " (coeficient of the matrix)" << std::endl;
	coutMaster << " test_graph_filename       = " << std::setw(20) << graph_filename << " (name of file with coordinates)" << std::endl;
	coutMaster << " test_graph_coeff          = " << std::setw(20) << graph_coeff << " (threshold of the graph)" << std::endl;
	coutMaster << " test_view_graph           = " << std::setw(20) << view_graph << " (print content of graph or not)" << std::endl;
	coutMaster << " test_view_matrix          = " << std::setw(20) << view_matrix << " (print matrix or not)" << std::endl;
	coutMaster << " test_view_matrix_unsorted = " << std::setw(20) << view_matrix_unsorted << " (print matrix in unsorted form or not)" << std::endl;

	/* start logging */
	std::ostringstream oss_name_of_file_log;
	oss_name_of_file_log << "results/test_mat_to_dense_log_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_name_of_file_log.str());

	/* say hello */
	coutMaster << "- start program" << std::endl;

	/* create graph */
	BGMGraph graph(graph_filename);
	graph.process(graph_coeff);

//	const double coordinates[2] = {1,2};
//	BGMGraph graph(coordinates,2,1);
//	graph.process(100);

//	const double coordinates[3] = {1,2,5};
//	BGMGraph graph(coordinates,3,1);
//	graph.process(2.5);

//	const double coordinates[1] = {1};
//	BGMGraph graph(coordinates,1,1);
//	graph.process_cpu(100);

	if(matrix_type==0 || matrix_type==1){
		/* graph matrix */

		/* print info about graph */
		if(!view_graph){
			graph.print(coutMaster);
		} else {
			graph.print_content(coutMaster);
		}

		R = graph.get_n();
	}
	if(matrix_type==2 || matrix_type==3){
		/* laplace matrix */
		R = 1;
	}


	/* create base vector */
	/* try to make a global vector of length T and then get size of local portion */
	Vec x_Vec;
	int Tlocal, Tbegin, Tend;
	Vec layout;
	TRY( VecCreate(PETSC_COMM_WORLD,&layout) );
	TRY( VecSetSizes(layout,PETSC_DECIDE,T) );
	TRY( VecSetFromOptions(layout) );
	/* get the ownership range - now I know how much I will calculate from the time-series */
	TRY( VecGetLocalSize(layout,&Tlocal) );
	TRY( VecGetOwnershipRange(layout,&Tbegin,&Tend) );
	TRY( VecDestroy(&layout) ); /* destroy testing vector - it is useless now */
		
	TRY( VecCreateMPI(PETSC_COMM_WORLD,R*K*Tlocal,R*K*T, &x_Vec) );

	GeneralVector<PetscVector> x(x_Vec);
	GeneralVector<PetscVector> y(x); /* result, i.e. y = A*x */

	/* create matrix */
	GeneralMatrix<PetscVector> *A;
	if(matrix_type==0){
		A = new BlockGraphFreeMatrix<PetscVector>(x, graph, K, alpha);
	}
	if(matrix_type==1){
		A = new BlockGraphSparseMatrix<PetscVector>(x, graph, K, alpha);
	}
	if(matrix_type==2){
		A = new BlockLaplaceFreeMatrix<PetscVector>(x, K, alpha);
	}
	if(matrix_type==3){
		A = new BlockLaplaceSparseMatrix<PetscVector>(x, K, alpha);
	}

//	A->print(coutMaster);
	A->printcontent(coutMaster);

	/* prepare timer */
	Timer timer1;
	timer1.restart();
	
	/* print full matrix as a multiplication A*e_n */
	double *values;
	double *gamma_arr;
	int k,r,t,col;

	/* print header of matrix print */
	if(view_matrix){
		coutMaster << "------------------------- A -------------------------" << std::endl;
	}

	/* go through rows and multiply with vector of standart basis */
	if(!view_matrix_unsorted){
		for(k=0; k < K; k++){
			for(r=0; r < R; r++){
				for(t=0;t < T; t++){
					/* prepare e_k in layout */
					TRY( VecSet(x_Vec,0) );

					TRY( VecGetArray(x_Vec,&gamma_arr) );
					if(t >= Tbegin && t < Tend){
						gamma_arr[(t*R + r)*K + k - Tbegin*R*K] = 1.0;
					}
					TRY( VecRestoreArray(x_Vec,&gamma_arr) );
					TRY( VecAssemblyBegin(x_Vec) );
					TRY( VecAssemblyEnd(x_Vec) );

					/* perform multiplication */
					timer1.start();
						y = (*A)*x;
					timer1.stop();

					/* print row (in fact, it is column, but matrix is symmetric) */
					if(view_matrix){
						/* get array of values */
						TRY( VecGetArray(y.get_vector(), &values) );
		
						/* print row */
						TRY( PetscPrintf(PETSC_COMM_WORLD, "%*d: ", 3, k*R*T + r*T + t) );

						/* interpret results from new layout to standart one */
						for(int k_col = 0; k_col < K; k_col++){
							for(int r_col = 0; r_col < R; r_col++){
								for(int t_col = 0; t_col < Tlocal; t_col++){
									col = t_col*R*K + r_col*K + k_col;
									if(abs(values[col]) > 0.000001){
										TRY( PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%*.*f,",6, 2, values[col]) );
									} else {
										TRY( PetscSynchronizedPrintf(PETSC_COMM_WORLD, "      ,") );
									}
								}
								TRY( PetscSynchronizedFlush(PETSC_COMM_WORLD, NULL) );
							}
						}
						TRY( PetscPrintf(PETSC_COMM_WORLD, "\n") );
					
						/* restore array with values */
						TRY( VecRestoreArray(y.get_vector(), &values) );
					}
				}
			}
		}
	} else {
		for(col=0; col < K*T*R; col++){
			/* prepare e_k in layout */
			TRY( VecSet(x_Vec,0) );

			TRY( VecGetArray(x_Vec,&gamma_arr) );
			if(col >= Tbegin*R*K && col < Tend*R*K){
				gamma_arr[col - Tbegin*R*K] = 1.0;
			}
			TRY( VecRestoreArray(x_Vec,&gamma_arr) );
			TRY( VecAssemblyBegin(x_Vec) );
			TRY( VecAssemblyEnd(x_Vec) );

			/* perform multiplication */
			timer1.start();
				y = (*A)*x;
			timer1.stop();

			/* print row (in fact, it is column, but matrix is symmetric) */
			if(view_matrix){
				/* get array of values */
				TRY( VecGetArray(y.get_vector(), &values) );
		
				/* print row */
				TRY( PetscPrintf(PETSC_COMM_WORLD, "%*d: ", 3, col) );

				for(int loccol = 0; loccol < K*R*Tlocal; loccol++){
					if(abs(values[loccol]) > 0.000001){
						TRY( PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%*.*f,",6, 2, values[loccol]) );
					} else {
						TRY( PetscSynchronizedPrintf(PETSC_COMM_WORLD, "      ,") );
					}
				}
				TRY( PetscSynchronizedFlush(PETSC_COMM_WORLD, NULL) );
				TRY( PetscPrintf(PETSC_COMM_WORLD, "\n") );
					
				/* restore array with values */
				TRY( VecRestoreArray(y.get_vector(), &values) );
			}
		}
	}

	/* print footer of matrix print */
	if(view_matrix){
		coutMaster << "-----------------------------------------------------" << std::endl;
	}

	/* print final info */
	coutMaster << "T = " << std::setw(9) << T << ", K = "<< std::setw(4) << K << ", R = "<< std::setw(4) << R << ", time = " << std::setw(10) << timer1.get_value_sum()/(double)(T*R*K) << std::endl;

	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize();

	return 0;
}

