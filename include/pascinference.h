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

	#ifdef USE_CUDA
	/* anselm hotfix */
		typedef struct _p_PetscCUDAIndices* PetscCUDAIndices;
		typedef struct _p_VecScatterCUDAIndices_StoS* VecScatterCUDAIndices_StoS;
		typedef struct _p_VecScatterCUDAIndices_PtoP* VecScatterCUDAIndices_PtoP;
		PETSC_EXTERN PetscErrorCode VecCUDACopyToGPUSome_Public(Vec,PetscCUDAIndices);
		PETSC_EXTERN PetscErrorCode VecCUDACopyFromGPUSome_Public(Vec,PetscCUDAIndices);
		PETSC_EXTERN PetscErrorCode VecScatterInitializeForGPU(VecScatter,Vec,ScatterMode);
		PETSC_EXTERN PetscErrorCode VecScatterFinalizeForGPU(VecScatter);
		PETSC_EXTERN PetscErrorCode VecCreateSeqCUDA(MPI_Comm,PetscInt,Vec*);
		PETSC_EXTERN PetscErrorCode VecCreateMPICUDA(MPI_Comm,PetscInt,PetscInt,Vec*);
	#endif

	/* PetscVector */
	#include "external/petsc/algebra/vector/petscvector.h"
#endif

/* cuda stuff */
#ifdef USE_CUDA
	#include <stdio.h> /* printf in cuda */

	#include <cuda.h>
	#include <cuda_runtime_api.h>
	#include <device_launch_parameters.h>
	#include <device_functions.h>
#endif

/* PERMON */
#ifdef USE_PERMON
	#include "fllopqp.h"
#endif

/* boost stuff */
#ifdef USE_BOOST
	#ifdef USE_CUDA
		#undef _GLIBCXX_ATOMIC_BUILTINS
	#endif

	/* load console parameters with boost */
	#include <boost/program_options.hpp>
	#include <boost/filesystem.hpp>

#endif

/* include metis */
#ifdef USE_METIS
	#include "metis.h"
#endif

/** 
*  \namespace pascinference
*  \brief main namespace of the library
*
*/
namespace pascinference {
	/** 
	*  \namespace pascinference::common
	*  \brief general commonly-used stuff
	*
	*/
	namespace common {}

	/** 
	*  \namespace pascinference::algebra
	*  \brief for manipulation with algebraic structures
	*
	*/
	namespace algebra {}

	/** 
	*  \namespace pascinference::data
	*  \brief for manipulation with problem data
	*
	*/
	namespace data {}

	/** 
	*  \namespace pascinference::solver
	*  \brief solvers for different types of problems
	*
	*/
	namespace solver {}

	/** 
	*  \namespace pascinference::model
	*  \brief models for solving time-series problems
	*
	*/
	namespace model {}


	// TODO: is this good idea?
	using namespace common;
	using namespace algebra;
	using namespace data;
	using namespace solver;
	using namespace model;
	
}

#include "common/common.h"

#include "algebra/arrayoperation.h"
#include "algebra/matrix/generalmatrix.h"
#include "algebra/vector/generalvector.h"
#include "algebra/feasibleset/generalfeasibleset.h"
#include "common/decomposition.h"

#include "data/generaldata.h"
#include "solver/list.h"
#include "solver/generalsolver.h"
#include "model/generalmodel.h"


#endif


