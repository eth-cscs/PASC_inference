/** @file common.h
 *  @brief commonly used stuff
 *
 *  This file includes commonly used classes and functions.
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_COMMON_H
#define	PASC_COMMON_H

#define RANDOM_BY_TIME false  /* if false, then random generator is initialized subject to time, else generated random data are always the same */
#define EXPORT_SAVEVTK true /* export solution to VTK */ //TODO: old and unnecessary?
#define EXPORT_SAVEVTK_filename "output/data.vtk" /* name of file to export VTK */ //TODO: old and unnecessary?

#include "sys/types.h"
#include "sys/sysinfo.h"
#include <iostream>
#include <iomanip>
#include <typeinfo> 
#include <cmath>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>
#include <stack>
#include <limits>

#ifdef USE_PETSC
	#include "petsc.h"
	#include "petscsys.h"

	#include "external/petsc/algebra/vector/petscvector.h"
#endif

#ifdef USE_CUDA
 #include "external/cuda/common/common.cuh"

 /* include additional petsc-cuda stuff */
 #ifdef USE_PETSC
  #include "petsccuda.h"
 #endif
#endif

#include "common/timer.h"
#include "common/memorycheck.h"
//#include "common/powercheck.h"
#include "common/globalmanager.h"
#include "common/consoleoutput.h"
#include "common/consoleinput.h"
#include "common/logging.h"
#include "common/mvnrnd.h"
#include "common/shortinfo.h"



namespace pascinference {
namespace common {

#ifdef USE_PETSC
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

void allbarrier();

void myround(double in, double *out);
std::string printbool(bool input);
void arg_parse(const char *args, int *argc, char ***argv);

std::vector<std::string> split(const std::string &s, char delim);
template<typename Out> void split(const std::string &s, char delim, Out result);


}
} /* end of namespace */

#endif
