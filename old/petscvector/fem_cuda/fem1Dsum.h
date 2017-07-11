#ifndef PASC_PETSCVECTOR_FEM1DSUM_H
#define	PASC_PETSCVECTOR_FEM1DSUM_H

#include "external/petscvector/common/decomposition.h"
#include "external/petscvector/algebra/vector/generalvector.h"
#include "general/algebra/fem/fem1Dsum.h"

namespace pascinference {
namespace common {

/* external-specific stuff */
template<> class Fem1DSum<PetscVector>::ExternalContent {
	public:
	
	#ifdef USE_CUDA
		int blockSize_reduce; /**< block size returned by the launch configurator */
		int minGridSize_reduce; /**< the minimum grid size needed to achieve the maximum occupancy for a full device launch */
		int gridSize_reduce; /**< the actual grid size needed, based on input size */

		int blockSize_prolongate; /**< block size returned by the launch configurator */
		int minGridSize_prolongate; /**< the minimum grid size needed to achieve the maximum occupancy for a full device launch */
		int gridSize_prolongate; /**< the actual grid size needed, based on input size */

		void cuda_occupancy();
		void cuda_reduce_data(double *data1_arr, double *data2_arr, int T1, int T2, int T2local, double diff);
		void cuda_prolongate_data(double *data1_arr, double *data2_arr, int T1, int T2, int T2local, double diff);
	#endif
};


template<> void Fem1DSum<PetscVector>::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const;
template<> void Fem1DSum<PetscVector>::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const;

template<> Fem1DSum<PetscVector>::ExternalContent * Fem1DSum<PetscVector>::get_externalcontent() const;
template<> void Fem1DSum<PetscVector>::compute_decomposition_reduced();

}
} /* end of namespace */

#endif
