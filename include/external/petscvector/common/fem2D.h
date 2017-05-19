#ifndef PASC_PETSCVECTOR_FEM2D_H
#define	PASC_PETSCVECTOR_FEM2D_H

#include "general/common/fem2D.h"
#include "external/petscvector/common/fem.h"

namespace pascinference {
namespace common {

/* external-specific stuff */
template<> class Fem2D<PetscVector>::ExternalContent : Fem<PetscVector>::ExternalContent {
	#ifdef USE_CUDA
		void cuda_reduce_data(double *data1, double *data2, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int left_t2_idx, double diff);
		void cuda_prolongate_data(double *data1, double *data2, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int left_t2_idx, double diff);

	#endif
};

template<> Fem2D<PetscVector>::Fem2D(Decomposition<PetscVector> *decomposition1, Decomposition<PetscVector> *decomposition2, double fem_reduce);
template<> void Fem2D<PetscVector>::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const;
template<> void Fem2D<PetscVector>::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const;
template<> void Fem2D<PetscVector>::compute_decomposition_reduced();
template<> Fem2D<PetscVector>::ExternalContent * Fem2D<PetscVector>::get_externalcontent() const;


}
} /* end of namespace */

#endif
