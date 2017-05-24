#ifndef PASC_PETSCVECTOR_FEMHAT_H
#define	PASC_PETSCVECTOR_FEMHAT_H

#include "general/common/femhat.h"
#include "external/petscvector/common/fem.h"

namespace pascinference {
namespace common {

/* external-specific stuff */
template<> class FemHat<PetscVector>::ExternalContent : Fem<PetscVector>::ExternalContent {
	#ifdef USE_CUDA
		void cuda_reduce_data(Vec &data1, Vec &data2, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int right_t1_idx, int left_t2_idx, int right_t2_idx, double diff);
		void cuda_prolongate_data(Vec &data1, Vec &data2, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int right_t1_idx, int left_t2_idx, int right_t2_idx, double diff);
	#endif
};

template<> FemHat<PetscVector>::FemHat(Decomposition<PetscVector> *decomposition1, Decomposition<PetscVector> *decomposition2, double fem_reduce);
template<> void FemHat<PetscVector>::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const;
template<> void FemHat<PetscVector>::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const;

}
} /* end of namespace */

#endif
