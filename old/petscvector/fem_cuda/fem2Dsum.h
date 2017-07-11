#ifndef PASC_PETSCVECTOR_FEM2DSUM_H
#define	PASC_PETSCVECTOR_FEM2DSUM_H

#include "general/algebra/fem/fem2Dsum.h"
#include "external/petscvector/algebra/fem/fem1Dsum.h"

namespace pascinference {
namespace common {

/* external-specific stuff */
template<> class Fem2DSum<PetscVector>::ExternalContent : public Fem1DSum<PetscVector>::ExternalContent {
	public:
	#ifdef USE_CUDA
		void cuda_occupancy();

		void cuda_reduce_data(double *data1_arr, double *data2_arr, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int left_t2_idx, double diff);
		void cuda_prolongate_data(double *data1_arr, double *data2_arr, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int left_t2_idx, double diff);
	#endif
};

template<> Fem2DSum<PetscVector>::Fem2DSum(Decomposition<PetscVector> *decomposition1, Decomposition<PetscVector> *decomposition2, double fem_reduce);
template<> void Fem2DSum<PetscVector>::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const;
template<> void Fem2DSum<PetscVector>::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const;
template<> void Fem2DSum<PetscVector>::compute_decomposition_reduced();
template<> Fem2DSum<PetscVector>::ExternalContent * Fem2DSum<PetscVector>::get_externalcontent() const;


}
} /* end of namespace */

#endif
