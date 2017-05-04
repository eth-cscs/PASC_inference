#ifndef PASC_PETSCVECTOR_PERMONSOLVER_H
#define	PASC_PETSCVECTOR_PERMONSOLVER_H

#include "general/solver/permonsolver.h"
#include "external/petscvector/solver/qpsolver.h"
#include "external/petscvector/data/qpdata.h"

#ifndef USE_PERMON
	#error 'PERMONSOLVER cannot be used without -DUSE_PERMON=ON'
#else
	#include "fllopqp.h" /* manipulation with quadratic programming problems (QP) */
	#include "fllopqps.h" /* manipulation with solvers (QPS) */

	#include "algebra/feasibleset/simplex_lineqbound.h"
#endif

// TODO:
#ifdef USE_CUDA
    #include <../src/vec/vec/impls/seq/seqcuda/cudavecimpl.h>
#endif


namespace pascinference {
namespace solver {



}
} /* end namespace */

#endif
