#ifndef PASC_PETSCVECTOR_CGQPSOLVER_H
#define	PASC_PETSCVECTOR_CGQPSOLVER_H

#include "general/solver/cgqpsolver.h"
//#include "external/petscvector/solver/qpsolver.h"
//#include "external/petscvector/data/qpdata.h"

namespace pascinference {
namespace solver {

template<> void CGQPSolver<PetscVector>::allocate_temp_vectors();
template<> void CGQPSolver<PetscVector>::free_temp_vectors();

}
} /* end namespace */

#endif
