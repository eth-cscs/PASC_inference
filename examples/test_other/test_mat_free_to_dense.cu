/** @file test_mat_free_to_dense.cu
 *  @brief print free matrix as dense
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include "matrix/blockgraph.h"

#ifndef USE_PETSCVECTOR
 #error 'This example is for PETSCVECTOR'
#endif

typedef petscvector::PetscVector PetscVector;
 
using namespace pascinference;

void set_stride(Vec x_Vec){
	int low, high;
	TRY( VecGetOwnershipRange(x_Vec, &low, &high) );
	double *values;
	TRY( VecGetArray(x_Vec, &values) );
	for(int i = 0; i < high-low; i++){
		values[i] = low + i;
	}	
	TRY( VecRestoreArray(x_Vec, &values) );
}

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_graph_filename", boost::program_options::value< std::string >(), "name of file with coordinates")
		("test_K", boost::program_options::value<int>(), "number of clusters")
		("test_T", boost::program_options::value<int>(), "dimension of the problem")
		("test_view_graph", boost::program_options::value<bool>(), "print content of graph or not");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	/* load console arguments */
	int T, K, R;
	std::string filename;
	bool test_view_x, test_view_y, test_view_graph;
	
	consoleArg.set_option_value("test_T", &T, 10);
	consoleArg.set_option_value("test_K", &K, 1);
	consoleArg.set_option_value("test_graph_filename", &filename, "../data/Koordinaten_EEG_P.bin");
	consoleArg.set_option_value("test_view_x", &test_view_x, false);
	consoleArg.set_option_value("test_view_y", &test_view_y, false);
	consoleArg.set_option_value("test_view_graph", &test_view_graph, false);

	coutMaster << " T        = " << std::setw(7) << T << " (length of time-series)" << std::endl;
	coutMaster << " K        = " << std::setw(7) << K << " (number of clusters)" << std::endl;
	coutMaster << " filename = " << std::setw(7) << filename << " (name of file with coordinates)" << std::endl;

	/* start logging */
	std::ostringstream oss_name_of_file_log;
	oss_name_of_file_log << "results/projection_log_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_name_of_file_log.str());

	/* say hello */
	coutMaster << "- start program" << std::endl;

	/* create graph */
//	BGM_Graph graph(filename);
//	graph.process_cpu(2.5);

	const double coordinates[2] = {1,2};
	BGM_Graph graph(coordinates,2,1);
	graph.process_cpu(100);

//	const double coordinates[1] = {1};
//	BGM_Graph graph(coordinates,1,1);
//	graph.process_cpu(100);

	R = graph.get_n();

	/* print info about graph */
	if(!test_view_graph){
		graph.print(coutMaster);
	} else {
		graph.print_content(coutMaster);
	}

	/* create base vector */
	/* try to make a global vector of length T and then get size of local portion */
	Vec x_Vec;
	int Tlocal, Tbegin;
	Vec layout;
	TRY( VecCreate(PETSC_COMM_WORLD,&layout) );
	TRY( VecSetSizes(layout,PETSC_DECIDE,T) );
	TRY( VecSetFromOptions(layout) );
	/* get the ownership range - now I know how much I will calculate from the time-series */
	TRY( VecGetLocalSize(layout,&Tlocal) );
	TRY( VecGetOwnershipRange(layout,&Tbegin,NULL) );
	TRY( VecDestroy(&layout) ); /* destroy testing vector - it is useless now */
		
	TRY( VecCreateMPI(PETSC_COMM_WORLD,R*K*Tlocal,R*K*T, &x_Vec) );

	GeneralVector<PetscVector> x(x_Vec);
	GeneralVector<PetscVector> y(x); /* result, i.e. y = A*x */

	/* create matrix */
	BlockGraphMatrix<PetscVector> A(x, graph, K);
	A.print(coutMaster,coutAll);
	coutAll.synchronize();

	/* print full matrix as a multiplication A*e_n */
	double *values;
	int k,r,tlocal,idx;
	for(int row_idx = 0; row_idx < Tlocal*R*K; row_idx++){
		int k = (int)(row_idx/(double)(Tlocal*R));
		int r = (int)((row_idx-k*Tlocal*R)/(double)(Tlocal));
		int tlocal = row_idx - k*Tlocal*R - r*Tlocal;
		idx = Tbegin+tlocal;

		/* prepare e_k */
		TRY( VecSet(x_Vec,0) );
		TRY( VecSetValue(x_Vec, idx, 1.0, INSERT_VALUES) );
		TRY( VecAssemblyBegin(x_Vec) );
		TRY( VecAssemblyEnd(x_Vec) );

		y = A*x;
		
		/* get array of values */
		TRY( VecGetArray(y.get_vector(), &values) );
		
		/* print row */
		TRY( PetscPrintf(PETSC_COMM_WORLD, "%*d: ", 3, i) );

		for(int t=0;t<Tlocal*R*K;t++){
			if(abs(values[t]) > 0.0001){
				TRY( PetscSynchronizedPrintf(PETSC_COMM_WORLD, "%*.*f,", 7, 3, values[t]) );
			} else {
				TRY( PetscSynchronizedPrintf(PETSC_COMM_WORLD, "       ,") );
			}
		}
		TRY( PetscSynchronizedFlush(PETSC_COMM_WORLD, NULL) );

		TRY( PetscPrintf(PETSC_COMM_WORLD, "\n", i) );


		/* restore array with values */
		TRY( VecRestoreArray(y.get_vector(), &values) );

	}

	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize();

	return 0;
}

