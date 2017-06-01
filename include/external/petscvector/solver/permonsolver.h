/** @file permonsolver.h
 *  @brief solve QP with solvers implemented in Permon library
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_PETSCVECTOR_PERMONSOLVER_H
#define	PASC_PETSCVECTOR_PERMONSOLVER_H

#include <iostream>

#include "general/solver/permonsolver.h"

#include "external/petscvector/common/common.h"
#include "external/petscvector/algebra/feasibleset/simplex_lineqbound.h"

//#include "external/petscvector/solver/qpsolver.h"
//#include "external/petscvector/data/qpdata.h"
//#include "external/petscvector/solver/spg_fs.h"
#include "external/petscvector/algebra/matrix/blockgraphsparse.h"

#ifdef USE_PERMON
/* if permon is not defined, then this class does not make any sence */

#include "fllopqp.h" /* manipulation with quadratic programming problems (QP) */
#include "fllopqps.h" /* manipulation with solvers (QPS) */

namespace pascinference {
namespace solver {

/* external-specific stuff */
template<> class PermonSolver<PetscVector>::ExternalContent {
	public:
		QPData<PetscVector> *qpdata;
		QP qp;						/**< Quadratic Programming problem */
		QPS qps;					/**< Quadratic Programming solver */
		
		/** @brief dump data of the solver
		 * 
		 * called before solve()
		 */
		void dump() const; 
};


template<> PermonSolver<PetscVector>::PermonSolver(QPData<PetscVector> &new_qpdata);
template<> void PermonSolver<PetscVector>::solve();



template<> PermonSolver<PetscVector>::ExternalContent * PermonSolver<PetscVector>::get_externalcontent() const;


}
} /* end namespace */

#endif

#endif
