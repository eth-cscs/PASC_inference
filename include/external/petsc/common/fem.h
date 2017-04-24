#ifndef PASC_FEM_PETSC_H
#define	PASC_FEM_PETSC_H

#include "algebra/vector/generalvector.h"
#include "common/fem.h"

namespace pascinference {
namespace common {

template<> void Fem<PetscVector>::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const;
template<> void Fem<PetscVector>::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const;

}
} /* end of namespace */

#endif
