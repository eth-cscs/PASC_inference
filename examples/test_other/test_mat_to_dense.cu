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
		("test_DDT", boost::program_options::value<int>(), "decomposition in time [int]")
		("test_DDR", boost::program_options::value<int>(), "decomposition in space [int]")
		("test_alpha", boost::program_options::value<double>(), "coefficient of the matrix [double]")
		("test_graph_filename", boost::program_options::value< std::string >(), "name of file with coordinates [string]")
		("test_graph_coeff", boost::program_options::value<double>(), "threshold of the graph [double]")
		("test_view_graph", boost::program_options::value<bool>(), "print content of graph or not [bool]")
		("test_view_decomposition", boost::program_options::value<bool>(), "print content of decomposition or not [bool]")
		("test_view_matrix", boost::program_options::value<bool>(), "print dense matrix or not [bool]")
		("test_view_matrix_unsorted", boost::program_options::value<bool>(), "print dense matrix in original unsorted form or not [bool]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	/* load console arguments */
	int T, K, matrix_type, DDT_size, DDR_size;
	std::string graph_filename;
	bool view_matrix, view_graph, view_matrix_unsorted, view_decomposition;
	double alpha, graph_coeff;
	int nproc = GlobalManager.get_size();
	
	consoleArg.set_option_value("test_matrix_type", &matrix_type, 1);
	consoleArg.set_option_value("test_T", &T, 1);
	consoleArg.set_option_value("test_K", &K, 1);
	consoleArg.set_option_value("test_DDT", &DDT_size, nproc);
	consoleArg.set_option_value("test_DDR", &DDR_size, 1);
	consoleArg.set_option_value("test_alpha", &alpha, 1.0);
	consoleArg.set_option_value("test_graph_filename", &graph_filename, "data/graph_twonodes.bin");
	consoleArg.set_option_value("test_graph_coeff", &graph_coeff, 1.1);
	consoleArg.set_option_value("test_view_matrix", &view_matrix, true);
	consoleArg.set_option_value("test_view_matrix_unsorted", &view_matrix_unsorted, false);
	consoleArg.set_option_value("test_view_graph", &view_graph, false);
	consoleArg.set_option_value("test_view_decomposition", &view_decomposition, false);

	/* print settings */
	coutMaster << " test_matrix_type          = " << std::setw(30) << matrix_type << " (which matrix type to test)\n";
	coutMaster << " test_T                    = " << std::setw(30) << T << " (length of time-series)\n";
	coutMaster << " test_K                    = " << std::setw(30) << K << " (number of clusters)\n";
	coutMaster << " test_DDT                  = " << std::setw(30) << DDT_size << " (decomposition in time)\n";
	coutMaster << " test_DDR                  = " << std::setw(30) << DDR_size << " (decomposition in space)\n";
	coutMaster << " test_alpha                = " << std::setw(30) << alpha << " (coeficient of the matrix)\n";
	coutMaster << " test_graph_filename       = " << std::setw(30) << graph_filename << " (name of file with coordinates)\n";
	coutMaster << " test_graph_coeff          = " << std::setw(30) << graph_coeff << " (threshold of the graph)\n";
	coutMaster << " test_view_graph           = " << std::setw(30) << view_graph << " (print content of graph or not)\n";
	coutMaster << " test_view_decomposition   = " << std::setw(30) << view_decomposition << " (print content of decomposition or not)\n";
	coutMaster << " test_view_matrix          = " << std::setw(30) << view_matrix << " (print matrix or not)\n";
	coutMaster << " test_view_matrix_unsorted = " << std::setw(30) << view_matrix_unsorted << " (print matrix in unsorted form or not)\n";

	if(DDT_size*DDR_size != nproc){
		coutMaster << "Sorry, DDT*DDR != nproc\n";
		coutMaster << " DDT   = " << DDT_size << "\n";
		coutMaster << " DDR   = " << DDR_size << "\n";
		coutMaster << " nproc = " << nproc << "\n";
		
		return 0;
	}

	/* start logging */
	std::ostringstream oss_name_of_file_log;
	oss_name_of_file_log << "results/test_mat_to_dense_log_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_name_of_file_log.str());

	/* say hello */
	coutMaster << "- start program\n";

	/* create graph */
	BGMGraph graph(graph_filename);
	graph.process(graph_coeff);

//	BGMGraphGrid2D graph(4,4);
//	graph.process_grid();

	/* create decomposition */
	Decomposition *decomposition;

	if(matrix_type==0 || matrix_type==1){ /* graph matrix */
		decomposition = new Decomposition(T, graph, K, 1, DDT_size, DDR_size);

		/* print info about graph */
		if(!view_graph){
			graph.print(coutMaster);
		} else {
			graph.print_content(coutMaster);
		}
	}
	
	if(matrix_type==2 || matrix_type==3){ /* laplace matrix */
		decomposition = new Decomposition(T, 1, K, 1, DDT_size);

		/* print info about graph */
		if(!view_graph){
			graph.print(coutMaster);
		} else {
			graph.print_content(coutMaster);
		}
	}

	/* print info about decomposition */
	if(!view_decomposition){
		decomposition->print(coutMaster);
	} else {
		decomposition->print_content(coutMaster, coutAll);
	}

	/* create vector x */
	Vec x_Vec;
	decomposition->createGlobalVec_gamma(&x_Vec);
	GeneralVector<PetscVector> x(x_Vec);
	
	GeneralVector<PetscVector> y(x); /* result, i.e. y = A*x */

	/* create matrix */
	GeneralMatrix<PetscVector> *A;
	if(matrix_type==0){
//		A = new BlockGraphFreeMatrix<PetscVector>(*decomposition, alpha);
	}
	if(matrix_type==1){
		A = new BlockGraphSparseMatrix<PetscVector>(*decomposition, alpha);
	}
	if(matrix_type==2){
		A = new BlockLaplaceFreeMatrix<PetscVector>(*decomposition, alpha);
	}
	if(matrix_type==3){
		A = new BlockLaplaceSparseMatrix<PetscVector>(*decomposition, alpha);
	}

//	A->print(coutMaster);
	A->printcontent(coutMaster);

	/* prepare timer */
	Timer timer1;
	timer1.restart();
	
	/* print full matrix as a multiplication A*e_n */
	double *values;
	double *gamma_arr;
	int R = decomposition->get_R();
	int k,r,t,col;

	/* print header of matrix print */
	if(view_matrix){
		coutMaster << "------------------------- A -------------------------\n";
	}

	/* go through rows and multiply with vector of standart basis */
	if(!view_matrix_unsorted){
		for(k=0; k < K; k++){
			for(r=0; r < R; r++){
				for(t=0;t < T; t++){
					/* prepare e_k in layout */
					TRY( VecSet(x_Vec,0) );
					TRY( VecSetValue(x_Vec, decomposition->get_idxglobal(t,r,k), 1.0, INSERT_VALUES) );
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
								for(int t_col = 0; t_col < T; t_col++){
									/* it is my element? */
									int Pr_col = decomposition->get_Pr(r_col);
									if( t_col >= decomposition->get_Tbegin() && t_col < decomposition->get_Tend()
										&& Pr_col >= decomposition->get_Rbegin() && Pr_col < decomposition->get_Rend()){
										int t_local = t_col - decomposition->get_Tbegin();
										int r_local = Pr_col - decomposition->get_Rbegin();
										double value = values[t_local*decomposition->get_Rlocal()*K + r_local*K + k_col];

										if(abs(value) > 0.000001){
											TRY( PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%*.*f,",6, 2, value) );
										} else {
											TRY( PetscSynchronizedPrintf(PETSC_COMM_WORLD, "      ,") );
										}
									}
									TRY( PetscSynchronizedFlush(PETSC_COMM_WORLD, NULL) );
								}
							}
						}
						TRY( PetscPrintf(PETSC_COMM_WORLD, "\n") );
					
						/* restore array with values */
						TRY( VecRestoreArray(y.get_vector(), &values) );
					}
				} /* end T */
			} /* end R */
		} /* end K */
	} else {
		for(col=0; col < K*T*R; col++){
			/* prepare e_k in layout */
			TRY( VecSet(x_Vec,0) );
			TRY( VecSetValue(x_Vec, col, 1.0, INSERT_VALUES) );
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

				for(int loccol = 0; loccol < K*decomposition->get_Rlocal()*decomposition->get_Tlocal(); loccol++){
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
		coutMaster << "-----------------------------------------------------\n";
	}

	/* print final info */
	coutMaster << " time = " << std::setw(10) << timer1.get_value_sum()/(double)(decomposition->get_T()*decomposition->get_R()*K) << "\n";

	/* say bye */	
	coutMaster << "- end program\n";

	logging.end();
	Finalize();

	return 0;
}

