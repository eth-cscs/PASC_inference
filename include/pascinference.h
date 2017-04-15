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
 #include "external/minlin/pascinference.h"
#endif

/* include Petsc */
#ifdef USE_PETSC
 #include "external/petsc/pascinference.h"
#endif

/* cuda stuff */
#ifdef USE_CUDA
 #include "external/cuda/pascinference.h"
#endif

/* PERMON */
#ifdef USE_PERMON
	#include "fllopqp.h"
#endif

/* boost stuff */
#ifdef USE_BOOST
 #include "external/boost/pascinference.h"
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


