/** @file pascinference.h
 *  @brief main header file
 *
 *  Include this header file to work with PASC_INFERENCE library.
 *
 *  @author Lukas Pospisil
 */

#ifndef PASCINFERENCE_H
#define	PASCINFERENCE_H


/* include common c++ header files */
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
#include <stack>
#include <limits>


/* include MINLIN */ 
#ifdef USE_MINLIN
	#include <minlin/minlin.h>
	#include <minlin/modules/threx/threx.h>

	namespace minlin {
		namespace threx {
			MINLIN_INIT
		}
	}
#endif

/* include Petsc */
#ifdef USE_PETSC
	#include "petsc.h"
#endif

/* PetscVector */
#ifdef USE_PETSCVECTOR
	#include "petscvector.h"
#endif

/* cuda stuff */
#ifdef USE_CUDA
	#include <stdio.h> /* printf in cuda */

	#include <cuda.h>
	#include <cuda_runtime_api.h>
	#include <device_launch_parameters.h>
	#include <device_functions.h>
#endif


/* include pascinference stuff */
#include "common.h"
#include "algebra.h"

#include "generaldata.h"
#include "generalsolver.h"
#include "generalfeasibleset.h"


#endif


