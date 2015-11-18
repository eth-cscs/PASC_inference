#include "common.h"
#include "problem.h"
#include "gamma.h"
#include "theta.h"
#include "savevtk.h"

/* include QPproblem */
#include "qpproblemprojectionstep.h"


PetscMPIInt proc_n, proc_id; /* for MPI_Comm functions */	

int main( int argc, char *argv[] )
{
	Initialize(argc,argv);

	/* variables */
	PetscErrorCode ierr;
	Data data;
	Gamma gamma;
	Theta theta;
	PetscInt s; /* index of main iterations */
	PetscReal L, L_old, deltaL; /* object function value */

	PetscLogDouble time_begin, time_end, time_elapsed; /* elapsed time of computation */

	
	/* viewer */
	PetscViewer my_viewer = PETSC_VIEWER_STDOUT_WORLD;

	/* say hello */	
	ierr = PetscViewerASCIIPrintf(my_viewer,"- start program:\n"); CHKERRQ(ierr);
	
	/* parameters of application */
	PetscInt dataN = 1000;
	PetscScalar eps_sqr = 10;
	/* load parameters from input */
	// TODO: move to settings
	ierr = PetscOptionsInt("-N","set the size of the problem","",dataN,&dataN,NULL); CHKERRQ(ierr);
	ierr = PetscOptionsReal("-eps_sqr","regularization parameter","",eps_sqr,&eps_sqr,NULL); CHKERRQ(ierr);
	
	
	
	/* generate problem */
	ierr = generate_problem(&data,dataN); CHKERRQ(ierr);
	/* print problem */
	if(PRINT_DATA){
		ierr = data.print(my_viewer); CHKERRQ(ierr);
	}	

	/* initialize gamma */
	ierr = gamma.init(data, gammaK); CHKERRQ(ierr);

	/* prepare gammas */
	ierr = gamma.prepare_random(); CHKERRQ(ierr);	
	if(PRINT_DATA){ /* print gamma */
		ierr = gamma.print(my_viewer); CHKERRQ(ierr);
	}	

	/* initialize theta */
	ierr = theta.init(data,gamma); CHKERRQ(ierr);

	/* initialize QP solver */
	QPproblemProjectionstep qpproblem(&data,&gamma,&theta,10);
	gamma.qpproblem = &qpproblem;
	gamma.qpproblem->init();
	
	/* initialize value of object function */
	L = PETSC_MAX_REAL; // TODO: the computation of L should be done in the different way
	
	/* here start to measure time for computation */
	ierr = PetscTime(&time_begin);CHKERRQ(ierr);	
	
	/* main cycle */
	ierr = PetscViewerASCIIPrintf(my_viewer,"- run main cycle:\n"); CHKERRQ(ierr);
	ierr = PetscViewerASCIIPushTab(my_viewer); CHKERRQ(ierr);
	for(s=0;s<max_s_steps;s++){
		ierr = PetscViewerASCIIPrintf(my_viewer,"- s = %d:\n",s); CHKERRQ(ierr);
		ierr = PetscViewerASCIIPushTab(my_viewer); CHKERRQ(ierr);

		/* compute Theta */
		ierr = theta.compute(data,gamma); CHKERRQ(ierr);

		/* compute gamma */
		ierr = gamma.compute(data,theta); CHKERRQ(ierr);
		
		/* update value of object function */
		L_old = L;
		ierr = gamma.get_objectfunc_value(&L); CHKERRQ(ierr);
		deltaL = PetscAbsScalar(L - L_old);
		
		/* if L_new > L_old then make something with solver */
		if(L > L_old){
			ierr = gamma.correctsolver(L - L_old); CHKERRQ(ierr);
		}
		
		/* print info about cost function */
		if(PETSC_TRUE){ 
			ierr = PetscViewerASCIIPrintf(my_viewer,"- L_old       = %f:\n",L_old); CHKERRQ(ierr);
			ierr = PetscViewerASCIIPrintf(my_viewer,"- L           = %f:\n",L); CHKERRQ(ierr);
			ierr = PetscViewerASCIIPrintf(my_viewer,"- |L - L_old| = %f:\n",deltaL); CHKERRQ(ierr);
		}	

		/* end the main cycle if the change of function value is sufficient */
		if (deltaL < deltaL_eps){
			break;
		}
		
		ierr = PetscViewerASCIIPopTab(my_viewer); CHKERRQ(ierr);
	}
	ierr = PetscViewerASCIIPopTab(my_viewer); CHKERRQ(ierr);

	/* here stop to measure time for computation */
	ierr = PetscTime(&time_end);CHKERRQ(ierr);	

	/* save the solution to VTK */
	if(PETSC_TRUE){
		ierr = PetscViewerASCIIPrintf(my_viewer,"- save solution to VTK:\n"); CHKERRQ(ierr);
		ierr = save_VTK(data,gamma); CHKERRQ(ierr);
	}
	
	theta.finalize();
	gamma.finalize();
	data.finalize();
	
	/* print info about elapsed time */
	time_elapsed = time_end - time_begin;
	ierr = PetscViewerASCIIPrintf(my_viewer,"- time for computation: %f s\n", time_elapsed); CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(my_viewer,"- number of iterations: %d\n", s); CHKERRQ(ierr);
	ierr = PetscViewerASCIIPrintf(my_viewer,"- |L - L_old|         = %f:\n",deltaL); CHKERRQ(ierr);
	
	Finalize();
	return 0;
}

