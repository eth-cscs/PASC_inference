#ifndef PASC_PETSCVECTOR_FEM2D_H
#define	PASC_PETSCVECTOR_FEM2D_H

#include "general/common/fem2D.h"
#include "external/petscvector/common/fem.h"

namespace pascinference {
namespace common {

template<> void Fem2D<PetscVector>::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const;
template<> void Fem2D<PetscVector>::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const;

}
} /* end of namespace */

#endif
