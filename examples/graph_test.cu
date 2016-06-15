/** @file graph_test.cu
 *  @brief test the graph matrix multiplication on several architectures
 *
 *  Generate n random vectors of length KN and compute matrix-multiplication.
 *  Measure the computing time
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

void create_base_vector(Vec *x_Vec,int T, int R, int K){
	/* create base vector */
	/* try to make a global vector of length T and then get size of local portion */
	int Tlocal;
	Vec layout;
	TRY( VecCreate(PETSC_COMM_WORLD,&layout) );
	TRY( VecSetSizes(layout,PETSC_DECIDE,T) );
	TRY( VecSetFromOptions(layout) );
	/* get the ownership range - now I know how much I will calculate from the time-series */
	TRY( VecGetLocalSize(layout,&Tlocal) );
	TRY( VecDestroy(&layout) ); /* destroy testing vector - it is useless now */
	
	TRY( VecCreateMPI(PETSC_COMM_WORLD,R*K*Tlocal,R*K*T, x_Vec) );
}

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
//		("test_filename", boost::program_options::value<char*>(), "name of file with coordinates")
		("test_K", boost::program_options::value<int>(), "number of clusters")
		("test_T", boost::program_options::value<int>(), "dimension of the problem");
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	/* load console arguments */
	int T, K, R;
//	char filename[50];
	
	consoleArg.set_option_value("test_T", &T, 10);
	consoleArg.set_option_value("test_K", &K, 1);
//	consoleArg.set_option_value("filename", &filename, "Koordinaten_EEG2.txt");

	coutMaster << " T        = " << std::setw(7) << T << " (length of time-series)" << std::endl;
	coutMaster << " K        = " << std::setw(7) << K << " (number of clusters)" << std::endl;
//	coutMaster << " filename = " << std::setw(7) << filename << " (name of file with coordinates)" << std::endl;

	/* start logging */
	std::ostringstream oss_name_of_file_log;
	oss_name_of_file_log << "results/projection_log_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_name_of_file_log.str());

	/* say hello */
	coutMaster << "- start program" << std::endl;

	/* create graph */
	BGM_Graph graph("../data/Koordinaten_EEG_P.bin");
	graph.process_cpu(2.5);
	R = 1; //graph.get_n();

	/* create base vector */
	Vec x_Vec;
	create_base_vector(&x_Vec, T, R, K);
	GeneralVector<PetscVector> x(x_Vec);
	GeneralVector<PetscVector> y(x); /* result, i.e. y = A*x */

	/* set some values to x */
	set_stride(x_Vec);

	/* create matrix */
	BlockGraphMatrix<PetscVector> A(x, graph, K);

//	TRY( VecDestroy(&x_Vec) );

	A.print(coutMaster,coutAll);
	coutAll.synchronize();

	y = A*x;

	coutMaster << "----------------- x -------------" << std::endl;
	TRY( VecView(x.get_vector(), PETSC_VIEWER_STDOUT_WORLD) );

	coutMaster << "----------------- y -------------" << std::endl;
	TRY( VecView(y.get_vector(), PETSC_VIEWER_STDOUT_WORLD) );

	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize();

	return 0;
}

