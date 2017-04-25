#ifndef PASC_GENERALVECTOR_PETSC_H
#define	PASC_GENERALVECTOR_PETSC_H

#include "algebra/vector/generalvector.h"
#include "external/petsc/algebra/vector/petscvector.h"

typedef petscvector::PetscVector PetscVector;

namespace pascinference {
namespace algebra {

template<> void GeneralVector<PetscVector>::set_random();
template<> void GeneralVector<PetscVector>::set_random2();
template<> std::string GeneralVector<PetscVector>::get_name();


}
} /* end of namespace */

#endif
