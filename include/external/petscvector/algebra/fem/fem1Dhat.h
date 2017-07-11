#ifndef PASC_PETSCVECTOR_FEM1DHAT_H
#define	PASC_PETSCVECTOR_FEM1DHAT_H

#include "general/algebra/fem/fem1Dhat.h"

namespace pascinference {
namespace algebra {

template<> void Fem1DHat<PetscVector>::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const;
template<> void Fem1DHat<PetscVector>::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const;

}
} /* end of namespace */

#endif
