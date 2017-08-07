#ifndef PASC_PETSCVECTOR_ENTROPYDATA_H
#define	PASC_PETSCVECTOR_ENTROPYDATA_H

#include "general/data/entropydata.h"
#include "external/petscvector/algebra/vector/generalvector.h"

namespace pascinference {
namespace data {

/* external-specific stuff */
template<> class EntropyData<PetscVector>::ExternalContent {
	public:
		Vec *x_powers_Vecs;
};

template<> EntropyData<PetscVector>::EntropyData(Decomposition<PetscVector> *decomposition, int Km);
template<> EntropyData<PetscVector>::~EntropyData(); 

template<> void EntropyData<PetscVector>::set_x(GeneralVector<PetscVector> *x);

template<> void EntropyData<PetscVector>::compute_moments(GeneralVector<PetscVector> *moments);
template<> void EntropyData<PetscVector>::compute_residuum(GeneralVector<PetscVector> *residuum, GeneralVector<PetscVector> *integrals) const;

template<> EntropyData<PetscVector>::ExternalContent * EntropyData<PetscVector>::get_externalcontent() const;


}
} /* end namespace */

#endif
