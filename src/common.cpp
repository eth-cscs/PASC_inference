#include "common.h"

/*!
 * initialize the application
 */ 
void Initialize(int argc, char *argv[]){
	/* init the Petsc */
	PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);

}

/*!
 * final call of the application
 */ 
void Finalize(){
	/* finalize the Petsc */
	PetscFinalize();
}

