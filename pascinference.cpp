#include "common.h"
#include "problem.h"
#include "gamma.h"
#include "theta.h"
#include "savevtk.h"

PetscMPIInt proc_n, proc_id; /* for MPI_Comm functions */	

int main( int argc, char *argv[] )
{
	Initialize(argc,argv);

	/* variables */
	PetscErrorCode ierr;
	Data data;
	Gamma gamma;
	Theta theta;
	PetscInt s;

	PetscLogDouble time_begin, time_end, time_elapsed; /* elapsed time of computation */
	
	/* viewer */
	PetscViewer my_viewer = PETSC_VIEWER_STDOUT_WORLD;

	/* say hello */	
	ierr = PetscViewerASCIIPrintf(my_viewer,"- start program:\n"); CHKERRQ(ierr);
	
	/* generate problem */
	ierr = get_problem(&data); CHKERRQ(ierr);
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
		/* print Theta */
		if(PRINT_DATA){ /* print gamma */
			ierr = theta.print(my_viewer); CHKERRQ(ierr);
		}	

		/* compute gamma */
		ierr = gamma.compute(data,theta); CHKERRQ(ierr);
		
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
	
	Finalize();
	return 0;
}

