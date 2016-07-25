/** @file test_mat_blockgraph_free_to_dense.cu
 *  @brief print free matrix as dense
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include "matrix/blockgraphfree.h"

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
		("test_T", boost::program_options::value<int>(), "dimension of the problem [int]")
		("test_K", boost::program_options::value<int>(), "number of clusters [int]")
		("test_alpha", boost::program_options::value<double>(), "coefficient of the matrix [double]")
		("test_graph_filename", boost::program_options::value< std::string >(), "name of file with coordinates [string]")
		("test_graph_coeff", boost::program_options::value<double>(), "threshold of the graph [double]")
		("test_view_graph", boost::program_options::value<bool>(), "print content of graph or not [bool]")
		("test_view_matrix", boost::program_options::value<bool>(), "print dense matrix or not [bool]");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	/* load console arguments */
	int T, K, R;
	std::string graph_filename;
	bool view_matrix, view_graph;
	double alpha, graph_coeff;
	
	consoleArg.set_option_value("test_T", &T, 1);
	consoleArg.set_option_value("test_K", &K, 1);
	consoleArg.set_option_value("test_alpha", &alpha, 1.0);
	consoleArg.set_option_value("test_graph_filename", &graph_filename, "data/graph_twonodes.bin");
	consoleArg.set_option_value("test_graph_coeff", &graph_coeff, 1.1);
	consoleArg.set_option_value("test_view_matrix", &view_matrix, true);
	consoleArg.set_option_value("test_view_graph", &view_graph, false);

	/* print settings */
	coutMaster << " test_T              = " << std::setw(20) << T << " (length of time-series)" << std::endl;
	coutMaster << " test_K              = " << std::setw(20) << K << " (number of clusters)" << std::endl;
	coutMaster << " test_alpha          = " << std::setw(20) << alpha << " (coeficient of the matrix)" << std::endl;
	coutMaster << " test_graph_filename = " << std::setw(20) << graph_filename << " (name of file with coordinates)" << std::endl;
	coutMaster << " test_graph_coeff    = " << std::setw(20) << graph_coeff << " (threshold of the graph)" << std::endl;
	coutMaster << " test_view_graph     = " << std::setw(20) << view_graph << " (print content of graph or not)" << std::endl;
	coutMaster << " test_view_matrix    = " << std::setw(20) << view_matrix << " (print matrix or not)" << std::endl;

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

	R = graph.get_n();

	/* print info about graph */
	if(!view_graph){
		graph.print(coutMaster);
	} else {
		graph.print_content(coutMaster);
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
	BlockGraphFreeMatrix<PetscVector> A(x, graph, K, alpha);
	A.print(coutMaster,coutAll);
	coutAll.synchronize();

	/* prepare timer */
	Timer timer1;
	timer1.restart();
	
	/* print full matrix as a multiplication A*e_n */
	double *values;
	double *gamma_arr;
	int k,r,t,tlocal,idx;

	/* print header of matrix print */
	if(view_matrix){
		coutMaster << "------------------------- A -------------------------" << std::endl;
	}

	/* go through rows and multiply with vector of standart basis */
	for(int row_idx = 0; row_idx < T*R*K; row_idx++){
		/* ompute position in matrix */
		k = floor(row_idx/(double)(T*R));
		r = floor((row_idx-k*T*R)/(double)(T));
		t = row_idx - (k*R+r)*T;
		
		/* prepare e_k in layout */
		TRY( VecSet(x_Vec,0) );
		TRY( VecAssemblyBegin(x_Vec) );
		TRY( VecAssemblyEnd(x_Vec) );

		TRY( VecGetArray(x_Vec,&gamma_arr) );
		if(t >= Tbegin && t < Tend){
			tlocal = t - Tbegin;
			gamma_arr[tlocal + (k*R+r)*Tlocal] = 1;
		}
		TRY( VecRestoreArray(x_Vec,&gamma_arr) );
		TRY( VecAssemblyBegin(x_Vec) );
		TRY( VecAssemblyEnd(x_Vec) );

		/* perform multiplication */
		timer1.start();
			y = A*x;
		timer1.stop();

		/* print row (in fact, it is column, but matrix is symmetric) */
		if(view_matrix){
			/* get array of values */
			TRY( VecGetArray(y.get_vector(), &values) );
		
			/* print row */
			TRY( PetscPrintf(PETSC_COMM_WORLD, "%*d: ", 3, row_idx) );

			/* interpret results from new layout to standart one */
			for(k=0;k<K;k++){
				for(r=0;r<R;r++){
					for(t=0;t<Tlocal;t++){
						idx = t+(k*R+r)*Tlocal;
						if(abs(values[idx]) > 0.000001){
							TRY( PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%*.*f,",6, 2, values[idx]) );
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

