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
#include <cmath>
#include <math.h>
#include <iostream>
#include <string>
#include <fstream>
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
	#include "petscsys.h"
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

/* boost stuff */
#ifdef USE_BOOST
	#ifdef USE_CUDA
		#undef _GLIBCXX_ATOMIC_BUILTINS
	#endif

	/* load console parameters with boost */
	#include <boost/program_options.hpp>

#endif

/* include pascinference stuff */
#include "common/common.h"
#include "algebra/algebra.h"

#include "data/generaldata.h"
#include "solver/list.h"
#include "solver/generalsolver.h"
#include "feasibleset/generalfeasibleset.h"
#include "model/generalmodel.h"


#endif


