#ifndef PASC_PETSCVECTOR_FEM1DSUM_H
#define	PASC_PETSCVECTOR_FEM1DSUM_H

#include "general/algebra/fem/fem1Dsum.h"

namespace pascinference {
namespace algebra {

template<> void Fem1DSum<PetscVector>::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const;
template<> void Fem1DSum<PetscVector>::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const;

}
} /* end of namespace */

#endif
