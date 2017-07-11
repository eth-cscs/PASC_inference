#ifndef PASC_PETSCVECTOR_FEM2DSUM_H
#define	PASC_PETSCVECTOR_FEM2DSUM_H

#include "general/algebra/fem/fem2Dsum.h"

namespace pascinference {
namespace algebra {

template<> void Fem2DSum<PetscVector>::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const;
template<> void Fem2DSum<PetscVector>::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const;

}
} /* end of namespace */

#endif
