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
#define DEFAULT_T 500 /* default length of generated time serie */
#define DEFAULT_K 3 /* default number of clusters */

#define DEBUG_PRINTDATA false /* print values of all data */

#define ALGORITHM_deltaL_eps 0.0001 /*stopping criteria of outer main loop */
#define ALGORITHM_max_s_steps 100 /* max number of outer steps */
#define ALGORITHM_EPSSQUARE 10.0 /* FEM regularization parameter */
#define DEBUG_ALGORITHM_PRINTDATA false /* print values of Theta, Gamma, QPSolver during main cycle */
#define DEBUG_ALGORITHM_PRINTDATA_L true /* print descent of object function in main outer loop */
#define DEBUG_ALGORITHM_PRINTDATA_GAMMA false /* print values of Gamma during main cycle */
#define DEBUG_ALGORITHM_PRINTDATA_THETA false /* print values of Theta during main cycle */

int main( int argc, char *argv[] )
{
	/* parameters of application */
	int dataT = DEFAULT_T; // TODO: do it in a different way
	int gammaK = DEFAULT_K;
	int max_s_steps = ALGORITHM_max_s_steps;

	Initialize(argc,argv); // TODO: load parameters of problem from console input

	/* variables */
	Data data;
	Gamma gamma;
	Theta theta;
	QPSolver qpsolver(&data,&gamma,&theta, ALGORITHM_EPSSQUARE);

	int s; /* index of main iterations */
	Scalar L, L_old, Lgamma, Ltheta, deltaL; /* object function value */

	/* say hello */	
	Message("- start program");
	
	/* generate problem */
	generate_problem(&data,dataT);
	/* print problem */
	if(DEBUG_PRINTDATA){
		data.print();
	}	

	/* initialize gamma */
	gamma.init(data, gammaK);

	/* prepare gammas */
	gamma.prepare_random();	
	if(DEBUG_PRINTDATA){ /* print gamma */
		gamma.print();
	}

	/* initialize theta */
	theta.init(data,gamma);
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
		theta.compute(data,gamma);
		if(DEBUG_ALGORITHM_PRINTDATA_THETA || DEBUG_ALGORITHM_PRINTDATA){ /* print theta */
			theta.print(2);
		}
		qpsolver.compute_b();
		Ltheta = qpsolver.get_function_value();
		
		/* --- COMPUTE gamma --- */
		gamma.compute(&qpsolver,data,theta);
		if(DEBUG_ALGORITHM_PRINTDATA_GAMMA || DEBUG_ALGORITHM_PRINTDATA){
			qpsolver.print(2,false);
		}
		qpsolver.compute_b();
		Lgamma = qpsolver.get_function_value();

		/* compute stopping criteria */
		L_old = L;
		L = qpsolver.get_function_value();
		deltaL = abs(L - L_old);

		/* print info about cost function */
		if(DEBUG_ALGORITHM_PRINTDATA_L || DEBUG_ALGORITHM_PRINTDATA){
			Message_info_value("  - L_old       = ",L_old);
			Message_info_value("  - Ltheta      = ",Ltheta);
			Message_info_value("  - Lgamma      = ",Lgamma);
			Message_info_value("  - L           = ",L);
			Message_info_value("  - |L - L_old| = ",deltaL);
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
		save_VTK(data,gamma);
	}

	theta.finalize();
	gamma.finalize();
	data.finalize();
	qpsolver.finalize();
	
	/* print info about elapsed time and solution */
	Message_info("- final info:");
	Message_info_time(" - time for computation: ",timer.stop());
	Message_info_value(" - number of iterations: ",s);
	Message_info_value(" - |L - L_old| = ",deltaL);

	/* say bye */	
	Message("- end program");
	
	Finalize();
	return 0;
}

