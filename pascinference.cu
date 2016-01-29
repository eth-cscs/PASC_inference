/*******************************************************************************
PASC INFERENCE library
Lukas Pospisil, Illia Horenko, Patrick Gagliardini, Will Sawyer
USI Lugano, 2016
lukas.pospisil@usi.ch

*******************************************************************************/

#include "common.h"
#include "problem.h"
#include "gamma.h"
#include "theta.h"
#include "savevtk.h"
#include "qpsolver.h"

/* PROBLEM SETTINGS */
#define DEFAULT_T 5000 /* default length of generated time serie */
#define DEFAULT_K 3 /* default number of clusters */

#define DEBUG_PRINTDATA false /* print values of all data */

#define ALGORITHM_deltaL_eps 0.0001 /*stopping criteria of outer main loop */
#define ALGORITHM_max_s_steps 100 /* max number of outer steps */
#define ALGORITHM_EPSSQUARE 10.0 /* FEM regularization parameter */
#define DEBUG_ALGORITHM_PRINTDATA false /* print values of Theta, Gamma, QPSolver during main cycle */
#define DEBUG_ALGORITHM_PRINTDATA_L true /* print descent of object function in main outer loop */
#define DEBUG_ALGORITHM_PRINTDATA_QPIT false /* print number of QPSolver iteration in every outer step */
#define DEBUG_ALGORITHM_PRINTDATA_GAMMA false /* print values of Gamma during main cycle */
#define DEBUG_ALGORITHM_PRINTDATA_THETA false /* print values of Theta during main cycle */

int main( int argc, char *argv[] )
{
	Timer timer_all; /* global timer for whole application */
	Timer timer_problem; /* for generating the problem */
	Timer timer_gamma; /* for gamma manipulation */
	Timer timer_theta; /* for theta manipulation */
	Timer timer_saveVTK; /* for final export to VTK */

	timer_all.restart();
	timer_problem.restart();
	timer_gamma.restart();
	timer_theta.restart();
	timer_saveVTK.restart();


	timer_all.start(); /* here starts the timer for whole application */
	
	/* parameters of application */
	int dataT = DEFAULT_T; // TODO: do it in a different way
	int gammaK = DEFAULT_K;
	int max_s_steps = ALGORITHM_max_s_steps;

	Initialize(argc,argv); // TODO: load parameters of problem from console input

	/* variables */
	Data data;
	Gamma gamma;
	Theta theta;
	QPSolver qpsolver(&gamma, ALGORITHM_EPSSQUARE);

	int s; /* index of main iterations */
	Scalar L, L_old, deltaL; /* object function value */


	/* say hello */	
	Message("- start program");
	
	/* generate problem */
	timer_problem.start(); /* start timer for generating problem */
 	 generate_problem(&data,dataT);
	timer_problem.stop();
	Message_info_time(" - problem generated in: ",timer_problem.get_value_last());
	
	/* print problem */
	if(DEBUG_PRINTDATA){
		data.print();
	}	

	/* initialize gamma */
	timer_gamma.start(); /* start timer for initializing gamma */
	 gamma.init(data, gammaK);
	 gamma.prepare_random();	/* prepare gammas */
	timer_gamma.stop();
	Message_info_time(" - gamma generated in: ",timer_gamma.get_value_last());

	if(DEBUG_PRINTDATA){ /* print gamma */
		gamma.print();
	}


	/* initialize theta */
	timer_theta.start();
 	 theta.init(data,gamma);
	timer_theta.stop();
	Message_info_time(" - theta prepared in: ",timer_theta.get_value_last());
	
	if(DEBUG_PRINTDATA){ /* print theta */
		theta.print();
	}

	/* initialize QP solvers */
	qpsolver.init();
	if(DEBUG_PRINTDATA){ /* print state of qpsolver */
		qpsolver.print();
	}

	
	/* initialize value of object function */
	L = std::numeric_limits<Scalar>::max(); // TODO: the computation of L should be done in the different way
	
	/* main cycle */
	Message("- run main cycle:");
	for(s=0;s < max_s_steps;s++){
		Message_info_value(" - s = ",s);

		/* --- COMPUTE Theta --- */
		timer_theta.start(); /* start timer for solving Theta-problem */
		 theta.compute(data,gamma);
		timer_theta.stop();
		Message_info_time("  - theta problem solved in: ",timer_theta.get_value_last());

		if(DEBUG_ALGORITHM_PRINTDATA_THETA || DEBUG_ALGORITHM_PRINTDATA){ /* print theta */
			theta.print(2);
		}

		
		/* --- COMPUTE gamma --- */
		timer_gamma.start(); /* start timer for solving gamma-problem */
		 gamma.compute(&qpsolver,data,theta);
		timer_gamma.stop(); 
		Message_info_time("  - gamma problem solved in: ",timer_gamma.get_value_last());
		 
		if(DEBUG_ALGORITHM_PRINTDATA_GAMMA || DEBUG_ALGORITHM_PRINTDATA){
			qpsolver.print(2);
		}

		/* compute stopping criteria */
		L_old = L;
		L = qpsolver.get_function_value();
		deltaL = abs(L - L_old);

		/* print info about cost function */
		if(DEBUG_ALGORITHM_PRINTDATA_L || DEBUG_ALGORITHM_PRINTDATA){
//			Message_info_value("  - L_old       = ",L_old);
			Message_info_value("  - L           = ",L);
//			Message_info_value("  - |L - L_old| = ",deltaL);
		}	

		/* end the main cycle if the change of function value is sufficient */
		if (deltaL < ALGORITHM_deltaL_eps){
			break;
		}
		
	}
	Message("- main cycle finished");

	/* save the solution to VTK */
	if(EXPORT_SAVEVTK){
		Message("- save solution to VTK");
		timer_saveVTK.start();
		 save_VTK(data,gamma);
		timer_saveVTK.stop();
	}

	theta.finalize();
	gamma.finalize();
	data.finalize();
	qpsolver.finalize();

	/* here ends the application */
	timer_all.stop();
	
	/* print info about elapsed time and solution */
	Message_info(      "- final info:");
	Message_info_time( " - time all: ",timer_all.get_value_sum());
	Message_info_time( "  - time problem: ",timer_problem.get_value_sum());
	Message_info_time( "  - time gamma:   ",timer_gamma.get_value_sum());
	Message_info_time( "  - time theta:   ",timer_theta.get_value_sum());
	Message_info_time( "  - time saveVTK: ",timer_saveVTK.get_value_sum());
	Message_info_time( "  - time other:   ",timer_all.get_value_sum() - (timer_problem.get_value_sum() + timer_gamma.get_value_sum() + timer_theta.get_value_sum() +  timer_saveVTK.get_value_sum()));

	Message_info_value(" - number of outer iterations: ",s);
	Message_info_value(" - |L - L_old| = ",deltaL);
	Message_info(" - QPSolver:");
	Message_info_value("  - it =        ", qpsolver.get_it_all());
	Message_info_value("  - hessmult =  ", qpsolver.get_hessmult_all());
	Message_info_time( "  - time =      ", qpsolver.get_time_total());
	Message_info_time( "   - t_project =  ", qpsolver.get_time_projection());
	Message_info_time( "   - t_matmult =  ", qpsolver.get_time_matmult());
	Message_info_time( "   - t_dot =      ", qpsolver.get_time_dot());
	Message_info_time( "   - t_update =   ", qpsolver.get_time_update());
	Message_info_time( "   - t_stepsize = ", qpsolver.get_time_stepsize());
	Message_info_time( "   - t_fs =       ", qpsolver.get_time_fs());
	Message_info_time( "   - t_other =    ", qpsolver.get_time_total() - (qpsolver.get_time_projection() + qpsolver.get_time_matmult() + qpsolver.get_time_dot() + qpsolver.get_time_update() + qpsolver.get_time_stepsize() + qpsolver.get_time_fs()));

	/* say bye */	
	Message("- end program");

	
	Finalize();
	return 0;
}

