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


typedef petscvector::PetscVector PetscVector;

 
using namespace pascinference;

int main( int argc, char *argv[] )
{
	Initialize(argc, argv); 

	/* read command line arguments */
	if(argc < 7){
		coutMaster << "1. 2. 3. argument - T_begin:T_step:T_end      - the dimension of the problem" << std::endl;
		coutMaster << "4. 5. 6. argument - K_begin:K_step:K_end      - the number of clusters" << std::endl;
		coutMaster << "7. argument       - n     					 - number of tests" << std::endl;
		coutMaster << std::endl << argv[0] << " T_begin T_step T_end K_begin K_step K_end n" << std::endl;
		return 1;
	}

	int T_begin = atoi(argv[1]);
	int T_step = atoi(argv[2]);
	int T_end = atoi(argv[3]);

	int K_begin = atoi(argv[4]);
	int K_step = atoi(argv[5]);
	int K_end = atoi(argv[6]);

	int n = atoi(argv[7]);
	coutMaster << "T_begin:T_step:T_end      = " << std::setw(7) << T_begin << std::setw(7) << T_step << std::setw(7) << T_end << " (length of time-series)" << std::endl;
	coutMaster << "K_begin:K_step:K_end      = " << std::setw(7) << K_begin << std::setw(7) << K_step << std::setw(7) << K_end << " (number of clusters)" << std::endl;

	coutMaster << "n       = " << std::setw(7) << n << " (number of tests)" << std::endl;
	

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
	GeneralVector<PetscVector> *x;

	/* prepare random generator */
	PetscRandom rnd;
	TRY( PetscRandomCreate(PETSC_COMM_WORLD,&rnd) );
	TRY( PetscRandomSetType(rnd,PETSCRAND) );
	TRY( PetscRandomSetFromOptions(rnd) );
	TRY( PetscRandomSetSeed(rnd,13) );


	for(K=K_begin;K<=K_end;K+=K_step){
	for(T=T_begin;T<=T_end;T+=T_step){

		TRY( VecCreateMPI(PETSC_COMM_WORLD,K*T,PETSC_DETERMINE, &x_global) );
		/* --- INITIALIZATION --- */

		/* prepare feasible set */
		feasibleset =  new SimplexFeasibleSet_Local(T,K);  	

		#ifdef USE_CUDA
			TRY( VecCreateSeqCUDA(PETSC_COMM_SELF, K*T, &x_local) );
		#else
			TRY( VecCreateSeq(PETSC_COMM_SELF, K*T, &x_local) );
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

