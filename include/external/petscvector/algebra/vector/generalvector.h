#ifndef PASC_PETSCVECTOR_GENERALVECTOR_H
#define	PASC_PETSCVECTOR_GENERALVECTOR_H

#include "general/algebra/vector/generalvector.h"
#include "external/petscvector/algebra/vector/petscvector.h"

typedef petscvector::PetscVector PetscVector;

namespace pascinference {
namespace algebra {

template<> void GeneralVector<PetscVector>::set_random();
template<> void GeneralVector<PetscVector>::set_random2();
template<> std::string GeneralVector<PetscVector>::get_name();


}
} /* end of namespace */

#endif
