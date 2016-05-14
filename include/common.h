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

#define RANDOM_BY_TIME false /* if false, then random generator is initialized subject to time, else generated random data are always the same */
#define EXPORT_SAVEVTK true /* export solution to VTK */
#define EXPORT_SAVEVTK_filename "output/data.vtk" /* name of file to export VTK */

#include "sys/types.h"
#include "sys/sysinfo.h"
#include <iostream>
#include <iomanip> /* setw - formated cout output */
#include <typeinfo> 

#include "common/timer.h"
#include "common/memorycheck.h"
#include "common/globalmanager.h"
#include "common/consoleoutput.h"
#include "common/arrayoperation.h"
#include "common/logging.h"


namespace pascinference {

/* global variables */
int DEBUG_MODE; /**< the debug mode of the library */

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
void Initialize(int, char**); 

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

#ifdef USE_PETSC
void MyVecDot(Vec x, Vec y, double *out){
	double *x_arr;
	double *y_arr;
	
	int size;
	TRY( VecGetSize(x, &size) );
	
	TRY( VecGetArray(x,&x_arr) );
	TRY( VecGetArray(y,&y_arr) );
	
	*out = dot_arrays(size, x_arr, y_arr);

	TRY( VecRestoreArray(x,&x_arr) );
	TRY( VecRestoreArray(y,&y_arr) );
}

/* w = x.*y */
void MyVecPointwiseMult(Vec w, Vec x, Vec y){
	double *w_arr;
	double *x_arr;
	double *y_arr;
	
	int size;
	TRY( VecGetSize(w, &size) );
	
	TRY( VecGetArray(w,&w_arr) );
	TRY( VecGetArray(x,&x_arr) );
	TRY( VecGetArray(y,&y_arr) );
	
	coutAll << "w: ";
	print_array(coutAll, size, w_arr);
	coutAll << std::endl;

	coutAll << "x: ";
	print_array(coutAll, size, x_arr);
	coutAll << std::endl;

	coutAll << "y: ";
	print_array(coutAll, size, x_arr);

	
	for(int i=0;i<=size;i++){
		w_arr[i] = x_arr[i]*y_arr[i];
	}	

	TRY( VecRestoreArray(w,&w_arr) );
	TRY( VecRestoreArray(x,&x_arr) );
	TRY( VecRestoreArray(y,&y_arr) );
}


#endif

} /* end of namespace */

#endif
