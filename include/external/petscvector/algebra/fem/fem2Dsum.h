#ifndef PASC_PETSCVECTOR_FEM2DSUM_H
#define	PASC_PETSCVECTOR_FEM2DSUM_H

#include "general/algebra/fem/fem2Dsum.h"

namespace pascinference {
namespace algebra {

/* external-specific stuff */
template<> class Fem2DSum<PetscVector>::ExternalContent {
	public:
};

template<> Fem2DSum<PetscVector>::Fem2DSum(Decomposition<PetscVector> *decomposition1, Decomposition<PetscVector> *decomposition2, double fem_reduce);
template<> void Fem2DSum<PetscVector>::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const;
template<> void Fem2DSum<PetscVector>::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const;
template<> void Fem2DSum<PetscVector>::compute_decomposition_reduced();
template<> Fem2DSum<PetscVector>::ExternalContent * Fem2DSum<PetscVector>::get_externalcontent() const;

}
} /* end of namespace */

#endif
