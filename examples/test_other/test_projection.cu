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
		("test_Tbegin", boost::program_options::value<int>(), "dimension of the problem")
		("test_Tstep", boost::program_options::value<int>(), "dimension of the problem")
		("test_Tend", boost::program_options::value<int>(), "dimension of the problem")
		("test_Kbegin", boost::program_options::value<int>(), "number of clusters")
		("test_Kstep", boost::program_options::value<int>(), "number of clusters")
		("test_Kend", boost::program_options::value<int>(), "number of clusters")
		("test_n", boost::program_options::value<int>(), "number of tests");	
	consoleArg.get_description()->add(opt_problem);

	/* call initialize */
	if(!Initialize(argc, argv)){
		return 0;
	} 

	/* load console arguments */
	int T_begin,T_step, T_end;
	int K_begin, K_step, K_end;
	int n;
	
	consoleArg.set_option_value("test_Tbegin", &T_begin, 2);
	consoleArg.set_option_value("test_Tend", &T_end, T_begin);
	consoleArg.set_option_value("test_Tstep", &T_step, 1);

	consoleArg.set_option_value("test_Kbegin", &K_begin, 3);
	consoleArg.set_option_value("test_Kend", &K_end, K_begin);
	consoleArg.set_option_value("test_Kstep", &K_step, 1);

	consoleArg.set_option_value("test_n", &n, 1);

	coutMaster << "T_begin:T_step:T_end      = " << std::setw(7) << T_begin << std::setw(7) << T_step << std::setw(7) << T_end << " (length of time-series)" << std::endl;
	coutMaster << "K_begin:K_step:K_end      = " << std::setw(7) << K_begin << std::setw(7) << K_step << std::setw(7) << K_end << " (number of clusters)" << std::endl;

	coutMaster << "n = " << std::setw(7) << n << " (number of tests)" << std::endl;

	/* start logging */
	std::ostringstream oss_name_of_file_log;
	oss_name_of_file_log << "results/projection_log_p" << GlobalManager.get_rank() << ".txt";
	logging.begin(oss_name_of_file_log.str());

	/* say hello */
	coutMaster << "- start program" << std::endl;

	Timer timer1; /* total projection time for petscseq */
	timer1.restart();	

	int K,T;
	Vec x_global;
	Vec x_local;

	SimplexFeasibleSet_Local *feasibleset;  	

	/* prepare random generator */
	PetscRandom rnd;
	TRY( PetscRandomCreate(PETSC_COMM_WORLD,&rnd) );
	TRY( PetscRandomSetType(rnd,PETSCRAND) );
	TRY( PetscRandomSetFromOptions(rnd) );
	TRY( PetscRandomSetSeed(rnd,13) );

	GeneralVector<PetscVector> *x;

	int Ti, Ki;
	for(Ki=K_begin;Ki<=K_end;Ki+=K_step){
	for(Ti=T_begin;Ti<=T_end;Ti+=T_step){

		K = Ki;
		T = pow(10,Ti);

		timer1.restart();

		/* create global vector */
		TRY( VecCreate(PETSC_COMM_WORLD, &x_global) );
		TRY( VecSetSizes(x_global, K*T, PETSC_DETERMINE) );
		TRY( VecSetType(x_global, VECMPI) );
		TRY( VecSetFromOptions(x_global) );
		TRY( VecAssemblyBegin(x_global) );
		TRY( VecAssemblyEnd(x_global) );
		
		/* --- INITIALIZATION --- */

		/* prepare feasible set */
		feasibleset =  new SimplexFeasibleSet_Local(T,K);  	

		#ifdef USE_CUDA
			TRY( VecCreateSeqCUDA(PETSC_COMM_SELF, K*T, &x_local) );
			gpuErrchk(cudaDeviceSynchronize());
		#else
			TRY( VecCreateSeq(PETSC_COMM_SELF, K*T, &x_local)  );
		#endif

		/* log */
		std::ostringstream oss_print_to_log;

		/* throught all sample cases */
		int j;
		for(j=0;j<n;j++){
			/* --- SET RANDOM VALUES TO GLOBAL VECTOR --- */
			TRY( VecSetRandom(x_global, rnd) );
		
			/* --- GET LOCAL VECTOR --- */
			TRY( VecGetLocalVector(x_global,x_local) );

			x = new GeneralVector<PetscVector>(x_local);

			/* --- COMPUTE PROJECTION --- */
			timer1.start();
				feasibleset->project(*x);
			timer1.stop();
	
			/* --- RESTORE GLOBAL VECTOR --- */
			TRY( VecRestoreLocalVector(x_global,x_local) );
		}

		/* --- PRINT INFO ABOUT TIMERS --- */
		coutMaster << "K = "<< std::setw(4) << K << ", T = " << std::setw(9) << T << ", time = " << std::setw(10) << timer1.get_value_sum() << std::endl;

		#ifdef USE_CUDA
			oss_print_to_log << "TIME_PETSCVECSEQ|";
		#else
			oss_print_to_log << "TIME_PETSCVECSEQCUDA|";
		#endif
		oss_print_to_log  << K << "|" << T << "|" << timer1.get_value_sum();

		LOG_DIRECT(oss_print_to_log.str());
		oss_print_to_log.str("");

		/* destroy feasible set */
		free(feasibleset);
	
		TRY( VecDestroy(&x_local) );
	
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

