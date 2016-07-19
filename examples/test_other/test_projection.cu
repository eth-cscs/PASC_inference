/** @file projection_test.cu
 *  @brief test the projection on several architectures
 *
 *  Generate n random vectors of length KN and compute projection.
 *  Measure the computing time
 *
 *  @author Lukas Pospisil
 */

#include "pascinference.h"
#include "feasibleset/simplex_local.h"


#ifndef USE_PETSCVECTOR
 #error 'This example is for PETSCVECTOR'
#endif

#include "solver/tssolver.h"

typedef petscvector::PetscVector PetscVector;

 
using namespace pascinference;

int main( int argc, char *argv[] )
{
	/* add local program options */
	boost::program_options::options_description opt_problem("PROBLEM EXAMPLE", consoleArg.get_console_nmb_cols());
	opt_problem.add_options()
		("test_T", boost::program_options::value<std::vector<int> >()->multitoken(), "dimensions of the problem")
		("test_K", boost::program_options::value<std::vector<int> >()->multitoken(), "number of clusters")
		("test_n", boost::program_options::value<int>(), "number of tests");	
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	/* --- LOAD CONSOLE ARGUMENTS --- */
	/* which times to try with projection */
	int T;
	std::vector<int> T_list;
	if(!consoleArg.set_option_value("test_T", &T_list)){
		std::cout << "test_T has to be set! Call application with parameter -h to see all parameters" << std::endl;
		return 0;
	}

	/* which number of clusters to try with projection */
	int K;
	std::vector<int> K_list;
	if(!consoleArg.set_option_value("test_K", &K_list)){
		std::cout << "test_K has to be set! Call application with parameter -h to see all parameters" << std::endl;
		return 0;
	}

	/* number of test for each combination T,K */
	int n;
	consoleArg.set_option_value("test_n", &n, 1);

	/* print loaded values */
	int i;
	int T_size = T_list.size();
	int K_size = K_list.size();
	std::ostringstream oss;

	/* print info about what we will compute */
	coutMaster << "- PROBLEM INFO --------------------------------------------------" << std::endl;
	coutMaster << " T      = ";
	for(i=0;i<T_size;i++){
		coutMaster << T_list[i];
		if(i < T_size-1){ 
				coutMaster << ", ";
		}
	}
	coutMaster << " (length of time-series)" << std::endl;
	coutMaster << " K      = ";
	for(i=0;i<K_size;i++){
		coutMaster << K_list[i];
		if(i < K_size-1){ 
				coutMaster << ", ";
		}
	}
	coutMaster << " (number of clusters)" << std::endl;
	coutMaster << " n      = " << n << " (number of tests)" << std::endl;
	coutMaster << "------------------------------------------------------------------" << std::endl;	


	/* start logging */
	std::ostringstream oss_name_of_file_log;
	oss_name_of_file_log << "results/projection_log_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_name_of_file_log.str());

	/* say hello */
	coutMaster << "- start program" << std::endl;

	Timer timer1; /* total projection time for petscseq */
	timer1.restart();	

	Vec x_global;

	SimplexFeasibleSet_Local *feasibleset;  	

	/* prepare random generator */
	PetscRandom rnd;
	TRY( PetscRandomCreate(PETSC_COMM_WORLD,&rnd) );
	TRY( PetscRandomSetType(rnd,PETSCRAND) );
	TRY( PetscRandomSetFromOptions(rnd) );
	TRY( PetscRandomSetSeed(rnd,13) );

	Vec layout;
	int Tlocal;
	GeneralVector<PetscVector> *x;

	int Ti, Ki;
	for(Ki = 0; Ki < K_size; Ki++){
		K = K_list[Ki];
		for(Ti = 0; Ti < T_size; Ti++){
			T = T_list[Ti];

			/* try to make a global vector of length T and then get size of local portion */
			TRY( VecCreate(PETSC_COMM_WORLD,&layout) );
			TRY( VecSetSizes(layout,PETSC_DECIDE,T) );
			TRY( VecSetFromOptions(layout) );
			/* get the ownership range - now I know how much I will calculate from the time-series */
			TRY( VecGetLocalSize(layout,&Tlocal) );
			TRY( VecDestroy(&layout) ); /* destroy testing vector - it is useless now */

			/* create global vector */
			#ifdef USE_CUDA
				TRY( VecCreate(PETSC_COMM_WORLD, &x_global) );
				TRY( VecSetSizes(x_global, K*Tlocal, PETSC_DETERMINE) );
				TRY( VecSetType(x_global, VECMPICUDA) );
				TRY( VecSetFromOptions(x_global) );
				TRY( VecAssemblyBegin(x_global) );
				TRY( VecAssemblyEnd(x_global) );
			#else
				TRY( VecCreate(PETSC_COMM_WORLD, &x_global) );
				TRY( VecSetSizes(x_global, K*Tlocal, PETSC_DETERMINE) );
				TRY( VecSetType(x_global, VECMPI) );
				TRY( VecSetFromOptions(x_global) );
				TRY( VecAssemblyBegin(x_global) );
				TRY( VecAssemblyEnd(x_global) );
			#endif

			/* create general vector from PETSc vector */
			x = new GeneralVector<PetscVector>(x_global);
		
			timer1.restart();

			/* prepare feasible set */
			feasibleset =  new SimplexFeasibleSet_Local(Tlocal,K);

			/* log */
			std::ostringstream oss_print_to_log;

			/* throught all sample cases */
			int j;
			for(j=0;j<n;j++){
				/* --- SET RANDOM VALUES TO GLOBAL VECTOR --- */
				TRY( VecSetRandom(x_global, rnd) );

				/* --- COMPUTE PROJECTION --- */
				timer1.start();
					feasibleset->project(*x);
				timer1.stop();

			}

			/* --- PRINT INFO ABOUT TIMERS --- */
			coutMaster << "K = "<< std::setw(4) << K << ", T = " << std::setw(9) << T << ", time = " << std::setw(10) << timer1.get_value_sum() << std::endl;

			#ifdef USE_CUDA
				oss_print_to_log << "TIME_PETSCVECMPI|";
			#else
				oss_print_to_log << "TIME_PETSCVECMPICUDA|";
			#endif
			oss_print_to_log  << K << "|" << T << "|" << timer1.get_value_sum();

			LOG_DIRECT(oss_print_to_log.str());
			oss_print_to_log.str("");

			/* destroy feasible set */
			free(feasibleset);
	
			/* destroy used vectors */
			TRY( VecDestroy(&x_global) );

		}
	}
	
	/* --- DESTROY --- */
	/* destroy the random generator */
	TRY( PetscRandomDestroy(&rnd) );

	/* say bye */	
	coutMaster << "- end program" << std::endl;

	logging.end();
	Finalize();

	return 0;
}

