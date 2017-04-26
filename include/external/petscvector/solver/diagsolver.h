#ifndef PASC_PETSCVECTOR_DIAGSOLVER_H
#define	PASC_PETSCVECTOR_DIAGSOLVER_H

#include "general/solver/diagsolver.h"
#include "external/petscvector/algebra/vector/generalvector.h"
//#include "external/petscvector/solver/generalsolver.h"
//#include "external/petscvector/data/diagdata.h"

namespace pascinference {
namespace solver {

template<> void DiagSolver<PetscVector>::solve();

}
} /* end namespace */

#endif
