#ifndef PASC_PETSCVECTOR_FEM1DHAT_H
#define	PASC_PETSCVECTOR_FEM1DHAT_H

#include "general/algebra/fem/fem1Dhat.h"

namespace pascinference {
namespace algebra {

/* external-specific stuff */
template<> class Fem1DHat<PetscVector>::ExternalContent {
	public:
};

template<> Fem1DHat<PetscVector>::Fem1DHat(Decomposition<PetscVector> *decomposition1, Decomposition<PetscVector> *decomposition2, double fem_reduce);
template<> void Fem1DHat<PetscVector>::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const;
template<> void Fem1DHat<PetscVector>::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const;
template<> void Fem1DHat<PetscVector>::compute_decomposition_reduced();

}
} /* end of namespace */

#endif
