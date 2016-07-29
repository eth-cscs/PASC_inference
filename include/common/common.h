/** @file common.h
 *  @brief commonly used stuff
 *
 *  This file includes commonly used classes and functions.
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_COMMON_H
#define	PASC_COMMON_H

/* general default values */
#define DEFAULT_DEBUG_MODE 0

#define RANDOM_BY_TIME false  /* if false, then random generator is initialized subject to time, else generated random data are always the same */
#define EXPORT_SAVEVTK true /* export solution to VTK */ //TODO: old and unnecessary?
#define EXPORT_SAVEVTK_filename "output/data.vtk" /* name of file to export VTK */ //TODO: old and unnecessary?

#include "sys/types.h"
#include "sys/sysinfo.h"
#include <iostream>
#include <iomanip> /* setw - formated cout output */
#include <typeinfo> 

#include "common/timer.h"
#include "common/memorycheck.h"
#include "common/globalmanager.h"
#include "common/consoleoutput.h"
#include "common/consoleinput.h"
#include "common/arrayoperation.h"
#include "common/logging.h"
#include "common/mvnrnd.h"

namespace pascinference {

/* global variables */
int DEBUG_MODE; /**< the debug mode of the library */ //TODO: old and unnecessary?

#ifdef USE_PETSCVECTOR
 extern bool PETSC_INITIALIZED;
#endif

/** @brief initialize the library
 * 
 *  Initialize random number generator. 
 *  If we are using Petsc, then Petsc is initialized.
 *
 * @todo process input arguments (use boost library?)
 */
bool Initialize(int, char**); 

/** @brief finalize the library
 * 
 *  If we are using Petsc, then Petsc is finalized.
 *
 */
void Finalize();

#ifdef USE_CUDA
	/* cuda error check */ 
	#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
	inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
	{
		if (code != cudaSuccess) 
		{
			fprintf(stderr,"\n\x1B[31mCUDA error:\x1B[0m %s %s \x1B[33m%d\x1B[0m\n\n", cudaGetErrorString(code), file, line);
			if (abort) exit(code);
		}
	}
#endif

}

/* ------------ IMPLEMENTATION -------------- */

namespace pascinference {

bool Initialize(int argc, char *argv[]){
	/* console arguments */
	if(!consoleArg.init(argc,argv)){
		return false;
	}

	/* initialize random seed: */
	if(RANDOM_BY_TIME){
		srand(time(NULL));
	} else {
		srand(0);
	}

  	/* init Petsc */
  	#ifdef USE_PETSC
//		PetscInitialize(&argc,&argv,PETSC_NULL,PETSC_NULL);
		PetscInitialize(PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);
		petscvector::PETSC_INITIALIZED = true;
	#endif
	
	return true;
}

void Finalize(){
  	/* finalize Petsc */
  	#ifdef USE_PETSC
		PetscFinalize();
		petscvector::PETSC_INITIALIZED = false;
	
	#endif

}

void myround(double in, double *out){
	union myUnion {
		double dValue;
		uint64_t iValue;
	} myValue;
	myValue.dValue=in;

//	myValue.iValue = myValue.iValue*0.001;
//	myValue.iValue = myValue.iValue*1000;

	*out = myValue.dValue;
}

} /* end of namespace */

#endif
