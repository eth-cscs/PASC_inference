#ifndef PASC_PETSCVECTOR_SIMPLEXFEASIBLESET_LOCAL_H
#define	PASC_PETSCVECTOR_SIMPLEXFEASIBLESET_LOCAL_H

#include "external/petscvector/algebra/vector/generalvector.h"
#include "general/algebra/feasibleset/simplex_local.h"

namespace pascinference {
namespace algebra {

template<> void SimplexFeasibleSet_Local<PetscVector>::project(GeneralVector<PetscVector> &x);

}
} /* end of namespace */


#endif
