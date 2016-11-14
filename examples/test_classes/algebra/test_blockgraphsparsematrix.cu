/** @file test_blockgraphsparsematrix.cpp
 *  @brief test class and methods: BlockGraphSparseMatrix
 *
 *  This file tests the class for manipulation with block diagonal matrix.
 * 
 *  @author Lukas Pospisil
 */

#include <iostream>
#include <list>
#include <algorithm>

#include "pascinference.h"
#include "algebra/matrix/blockgraphsparse.h"

typedef petscvector::PetscVector PetscVector;

using namespace pascinference;

extern int pascinference::DEBUG_MODE;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_T", boost::program_options::value<int>(), "length of time-series [int]")
		("test_K", boost::program_options::value<int>(), "number of clusters [int]")
		("test_DDT", boost::program_options::value<int>(), "decomposition in time [int]")
		("test_DDR", boost::program_options::value<int>(), "decomposition in space [int]")
		("test_alpha", boost::program_options::value<double>(), "coefficient of the matrix [double]")
		("test_graph_filename", boost::program_options::value< std::string >(), "name of file with coordinates [string]")
		("test_graph_dim", boost::program_options::value<int>(), "dimension of graph 1/2/3 [int]")
		("test_graph_coeff", boost::program_options::value<double>(), "threshold of the graph [double]")
		("test_print_graph", boost::program_options::value<bool>(), "print content of graph or not [bool]")
		("test_print_decomposition", boost::program_options::value<bool>(), "print content of decomposition or not [bool]")
		("test_print_matrix", boost::program_options::value<bool>(), "print dense matrix or not [bool]")
		("test_print_matrix_unsorted", boost::program_options::value<bool>(), "print dense matrix in original unsorted form or not [bool]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	/* load console arguments */
	int T, K, DDT_size, DDR_size, graph_dim;
	std::string graph_filename;
	bool print_matrix, print_graph, print_matrix_unsorted, print_decomposition;
	double alpha, graph_coeff;
	int nproc = GlobalManager.get_size();
	
	consoleArg.set_option_value("test_T", &T, 3);
	consoleArg.set_option_value("test_K", &K, 1);
	consoleArg.set_option_value("test_DDT", &DDT_size, nproc);
	consoleArg.set_option_value("test_DDR", &DDR_size, 1);
	consoleArg.set_option_value("test_alpha", &alpha, 1.0);
	consoleArg.set_option_value("test_graph_filename", &graph_filename, "data/graph_twonodes.bin");
	consoleArg.set_option_value("test_graph_dim", &graph_dim, 2);
	consoleArg.set_option_value("test_graph_coeff", &graph_coeff, 1.1);
	consoleArg.set_option_value("test_print_matrix", &print_matrix, true);
	consoleArg.set_option_value("test_print_matrix_unsorted", &print_matrix_unsorted, false);
	consoleArg.set_option_value("test_print_graph", &print_graph, false);
	consoleArg.set_option_value("test_print_decomposition", &print_decomposition, false);

	/* print settings */
	coutMaster << " test_T                     = " << std::setw(30) << T << " (length of time-series)" << std::endl;
	coutMaster << " test_K                     = " << std::setw(30) << K << " (number of clusters)" << std::endl;
	coutMaster << " test_DDT                   = " << std::setw(30) << DDT_size << " (decomposition in time)" << std::endl;
	coutMaster << " test_DDR                   = " << std::setw(30) << DDR_size << " (decomposition in space)" << std::endl;
	coutMaster << " test_alpha                 = " << std::setw(30) << alpha << " (coeficient of the matrix)" << std::endl;
	coutMaster << " test_graph_filename        = " << std::setw(30) << graph_filename << " (name of file with coordinates)" << std::endl;
	coutMaster << " test_graph_dim             = " << std::setw(30) << graph_dim << " (dimension of the graph 1/2/3)" << std::endl;
	coutMaster << " test_graph_coeff           = " << std::setw(30) << graph_coeff << " (threshold of the graph)" << std::endl;
	coutMaster << " test_print_graph           = " << std::setw(30) << print_graph << " (print content of graph or not)" << std::endl;
	coutMaster << " test_print_decomposition   = " << std::setw(30) << print_decomposition << " (print content of decomposition or not)" << std::endl;
	coutMaster << " test_print_matrix          = " << std::setw(30) << print_matrix << " (print matrix or not)" << std::endl;
	coutMaster << " test_print_matrix_unsorted = " << std::setw(30) << print_matrix_unsorted << " (print matrix in unsorted form or not)" << std::endl;
	coutMaster << std::endl;

	if(DDT_size*DDR_size != nproc){
		coutMaster << "Sorry, DDT*DDR != nproc" << std::endl;
		coutMaster << " DDT   = " << std::setw(15) << DDT_size << std::endl;
		coutMaster << " DDR   = " << std::setw(15) << DDR_size << std::endl;
		coutMaster << " nproc = " << std::setw(15) << nproc << std::endl;
		
		return 0;
	}

	/* we will measure the time of operations */
	Timer mytimer;
	mytimer.restart();
	
	/* create graph (i.e. load from filename) */
	mytimer.start();
	 BGMGraph graph(graph_filename,graph_dim);
	mytimer.stop();
	coutMaster << "- time load graph           : " << std::setw(15) << mytimer.get_value_last() << " s" << std::endl;
	
	/* process graph with given coefficient */
	mytimer.start();
	 graph.process(graph_coeff);
	mytimer.stop();
	coutMaster << "- time process graph        : " << std::setw(15) << mytimer.get_value_last() << " s" << std::endl;
	if(!print_graph){
		/* print info about graph */
		coutMaster.push();
		 graph.print(coutMaster);
		coutMaster.pop();
	} else {
		/* print detailed info (content) */
		coutMaster.push();
		 graph.print_content(coutMaster);
		coutMaster.pop();
	}

	/* create decomposition */
	Decomposition *decomposition;
	mytimer.start();
	 decomposition = new Decomposition(T, graph, K, 1, DDT_size, DDR_size);
	mytimer.stop();
	coutMaster << "- time create decomposition : " << std::setw(15) << mytimer.get_value_last() << " s" << std::endl;
	if(!print_decomposition){
		/* print info about decomposition */
		coutMaster.push();
		 decomposition->print(coutMaster);
		coutMaster.pop();
	} else {
		/* print content of decomposition */
		coutMaster.push();
		 decomposition->print_content(coutMaster, coutAll);
		coutMaster.pop();
	}

	/* create vectors x,y in tested multiplication y=A*x */
	Vec x_Vec; /* this is standart PETSc vector */
	mytimer.start();
	 decomposition->createGlobalVec_gamma(&x_Vec); /* this creates the layout of provided vector based on decomposition */
	 GeneralVector<PetscVector> x(x_Vec); /* initialize our PetscVector from provided standart PETSc vector */
	 GeneralVector<PetscVector> y(x); /* initialize result of multiplication duplicating provided vector */
	mytimer.stop();
	coutMaster << "- time create vectors       : " << std::setw(15) << mytimer.get_value_last() << " s" << std::endl;

	/* create matrix A */
	GeneralMatrix<PetscVector> *A;
	mytimer.start();
	 A = new BlockGraphSparseMatrix<PetscVector>(*decomposition, alpha);
	mytimer.stop();
	coutMaster << "- time create matrix        : " << std::setw(15) << mytimer.get_value_last() << " s" << std::endl;
	/* print informations about matrix */
	coutMaster.push();
	 A->print(coutMaster);
//	 A->printcontent(coutMaster);
	coutMaster.pop();

	/* we will print full matrix as a results of multiplication A*e_k, k = 1,...n */
	double *values_arr;
	double *gamma_arr;
	int R = decomposition->get_R(); /* we have T,K, why not to define also R - we obtain it from decomposed graph */

	/* print header of matrix print */
	if(print_matrix){
		coutMaster << "------------------------- A -------------------------" << std::endl << std::endl;;
	}

	/* go through rows and multiply with vector of standart basis */
	if(!print_matrix_unsorted){
		/* use index rearranging based on decomposition (see get_idxglobal() in following code) */
		for(int k=0; k < K; k++){ /* through clusters */
			for(int r=0; r < R; r++){ /* through graph nodes */
				for(int t=0;t < T; t++){ /* through time steps */
					/* prepare e_k in layout */
					TRYCXX( VecSet(x_Vec,0) );
					TRYCXX( VecSetValue(x_Vec, decomposition->get_idxglobal(t,r,k), 1.0, INSERT_VALUES) );
					TRYCXX( VecAssemblyBegin(x_Vec) );
					TRYCXX( VecAssemblyEnd(x_Vec) );

					/* perform multiplication */
					mytimer.start();
						y = (*A)*x; /* A is in this case a pointer, normally there will be y=A*x; */
					mytimer.stop();

					/* print row (in fact, it is column, but matrix is symmetric) */
					if(print_matrix){
						/* get array of values */
						TRYCXX( VecGetArray(y.get_vector(), &values_arr) );
		
						/* print row */
						TRYCXX( PetscPrintf(PETSC_COMM_WORLD, "%*d: ", 3, k*R*T + r*T + t) );

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
										double value = values_arr[t_local*decomposition->get_Rlocal()*K + r_local*K + k_col];

										if(abs(value) > 0.000001){
											TRYCXX( PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%*.*f,",6, 2, value) );
										} else {
											TRYCXX( PetscSynchronizedPrintf(PETSC_COMM_WORLD, "      ,") );
										}
									}
									TRYCXX( PetscSynchronizedFlush(PETSC_COMM_WORLD, NULL) );
								}
							}
						}
						TRYCXX( PetscPrintf(PETSC_COMM_WORLD, "\n") );
					
						/* restore array with values */
						TRYCXX( VecRestoreArray(y.get_vector(), &values_arr) );
					}
				} /* end T */
			} /* end R */
		} /* end K */
	} else {
		/* naiive printing - simply print a content without rearranging coeficient subject to decomposition */
		for(int col=0; col < K*T*R; col++){ /* simply through all columns */
			/* prepare e_k in layout */
			TRYCXX( VecSet(x_Vec,0) );
			TRYCXX( VecSetValue(x_Vec, col, 1.0, INSERT_VALUES) );
			TRYCXX( VecAssemblyBegin(x_Vec) );
			TRYCXX( VecAssemblyEnd(x_Vec) );

			/* perform multiplication */
			mytimer.start();
				y = (*A)*x;
			mytimer.stop();

			/* print row (in fact, it is column, but matrix is symmetric) */
			if(print_matrix){
				/* get array of values */
				TRYCXX( VecGetArray(y.get_vector(), &values_arr) );
		
				/* print row */
				TRYCXX( PetscPrintf(PETSC_COMM_WORLD, "%*d: ", 3, col) );

				for(int loccol = 0; loccol < K*decomposition->get_Rlocal()*decomposition->get_Tlocal(); loccol++){
					if(abs(values_arr[loccol]) > 0.000001){
						TRYCXX( PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%*.*f,",6, 2, values_arr[loccol]) );
					} else {
						TRYCXX( PetscSynchronizedPrintf(PETSC_COMM_WORLD, "      ,") );
					}
				}
				TRYCXX( PetscSynchronizedFlush(PETSC_COMM_WORLD, NULL) );
				TRYCXX( PetscPrintf(PETSC_COMM_WORLD, "\n") );
					
				/* restore array with values */
				TRYCXX( VecRestoreArray(y.get_vector(), &values_arr) );
			}
		}
	}

	/* print footer of matrix print */
	if(print_matrix){
		coutMaster << "-----------------------------------------------------" << std::endl << std::endl;;
	}

	/* print final info */
	coutMaster << "- time multiplication all   : " << std::setw(15) << mytimer.get_value_sum() << " s" << std::endl;
	coutMaster << "- time multiplication avg   : " << std::setw(15) << mytimer.get_value_sum()/(double)(decomposition->get_T()*decomposition->get_R()*K) << " s" << std::endl;

	
	/* say bye */	
	coutMaster << "- end of program" << std::endl;
	coutMaster << std::endl;

	Finalize();

	return 0;
}


