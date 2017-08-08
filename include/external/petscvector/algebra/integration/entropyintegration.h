#ifndef PASC_PETSCVECTOR_ENTROPYINTEGRATION_H
#define	PASC_PETSCVECTOR_ENTROPYINTEGRATION_H

#include "general/algebra/integration/entropyintegration.h"
#include "external/petscvector/algebra/vector/generalvector.h"

namespace pascinference {
namespace algebra {

template<> void EntropyIntegration<PetscVector>::compute(GeneralVector<PetscVector> &integrals);

}
} /* end of namespace */

#endif
