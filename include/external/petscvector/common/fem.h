#ifndef PASC_PETSCVECTOR_FEM_H
#define	PASC_PETSCVECTOR_FEM_H

#include "external/petscvector/common/decomposition.h"
#include "external/petscvector/algebra/vector/generalvector.h"
#include "general/common/fem.h"

namespace pascinference {
namespace common {

template<> void Fem<PetscVector>::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const;
template<> void Fem<PetscVector>::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const;

}
} /* end of namespace */

#endif
