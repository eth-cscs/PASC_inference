#ifndef PASC_PETSCVECTOR_FEM1DSUM_H
#define	PASC_PETSCVECTOR_FEM1DSUM_H

#include "general/algebra/fem/fem1Dsum.h"

namespace pascinference {
namespace algebra {

/* external-specific stuff */
template<> class Fem1DSum<PetscVector>::ExternalContent {
	public:
};

template<> void Fem1DSum<PetscVector>::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const;
template<> void Fem1DSum<PetscVector>::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const;

template<> Fem1DSum<PetscVector>::ExternalContent * Fem1DSum<PetscVector>::get_externalcontent() const;
template<> void Fem1DSum<PetscVector>::compute_decomposition_reduced();

}
} /* end of namespace */

#endif
