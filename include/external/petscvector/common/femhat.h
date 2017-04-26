#ifndef PASC_PETSCVECTOR_FEMHAT_H
#define	PASC_PETSCVECTOR_FEMHAT_H

#include "general/common/femhat.h"
#include "external/petscvector/common/fem.h"

namespace pascinference {
namespace common {

template<> void FemHat<PetscVector>::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const;
template<> void FemHat<PetscVector>::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const;

}
} /* end of namespace */

#endif
