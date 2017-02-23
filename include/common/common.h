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
#include "common/logging.h"
#include "common/mvnrnd.h"
#include "common/shortinfo.h"

namespace pascinference {
namespace common {

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

void myround(double in, double *out);
std::string printbool(bool input);
void arg_parse(const char *args, int *argc, char ***argv);

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
	
__global__ void kernel_warmup(){
	
}

/* include additional petsc-cuda stuff */
#ifdef USE_PETSC
	#include "petsccuda.h"
#endif

#endif

std::vector<std::string> split(const std::string &s, char delim);
template<typename Out> void split(const std::string &s, char delim, Out result);


}
} /* end of namespace */

/* ------------ IMPLEMENTATION -------------- */

namespace pascinference {
namespace common {

/* for loading PETSc options */
char **argv_petsc;
int argc_petsc;

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
		/* Petsc options are provided using "-petsc_options=string" */
		std::string petsc_options_string;
		consoleArg.set_option_value("petsc_options", &petsc_options_string, "");

		/* split string and create argc_petsc and argv_petsc */
		petsc_options_string.insert(0," "); /* add blank space in the beginning (?) */
		std::vector<std::string> petsc_options_vector = split(petsc_options_string, ' ');
		argc_petsc = petsc_options_vector.size();
		
		argv_petsc = new char*[argc_petsc]; /* the first parameter is ignored (?) */
		for(size_t i = 0; i < argc_petsc; ++i){
			argv_petsc[i] = new char[petsc_options_vector[i].size() + 1];
			std::strcpy(argv_petsc[i], petsc_options_vector[i].c_str());
		}

		#ifdef USE_PERMON
			FllopInitialize(&argc_petsc,&argv_petsc,PETSC_NULL);
//			FllopInitialize(PETSC_NULL,PETSC_NULL,PETSC_NULL);
		#else
			PetscInitialize(&argc_petsc,&argv_petsc,PETSC_NULL,PETSC_NULL);
//			PetscInitialize(PETSC_NULL,PETSC_NULL,PETSC_NULL,PETSC_NULL);
		#endif

		petscvector::PETSC_INITIALIZED = true;
	#endif
	
	/* cuda warm up */
	#ifdef USE_CUDA
		kernel_warmup<<<1,1>>>();
		gpuErrchk( cudaDeviceSynchronize() );
	#endif
	
	return true;
}

void Finalize(){
  	/* finalize Petsc */
  	#ifdef USE_PETSC
		/* clean memory of arguments */
		for(size_t i = 0; i < argc_petsc; ++i)
			delete[] argv_petsc[i];

		#ifdef USE_PERMON
			FllopFinalize();
		#else
			PetscFinalize();
		#endif

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

std::string printbool(bool input){
	std::string out;
	if(input){
		out = "true";
	} else {
		out = "false";
	}
	return out;
}

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss;
    ss.str(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}



}
} /* end of namespace */

#endif
