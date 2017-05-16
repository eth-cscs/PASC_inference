#ifndef PASC_PETSCVECTOR_FEM2D_H
#define	PASC_PETSCVECTOR_FEM2D_H

#include "general/common/fem2D.h"
#include "external/petscvector/common/fem.h"

namespace pascinference {
namespace common {

/* external-specific stuff */
template<> class Fem2D<PetscVector>::ExternalContent {
	#ifdef USE_CUDA
		int blockSize_reduce; /**< block size returned by the launch configurator */
		int minGridSize_reduce; /**< the minimum grid size needed to achieve the maximum occupancy for a full device launch */
		int gridSize_reduce; /**< the actual grid size needed, based on input size */

		int blockSize_prolongate; /**< block size returned by the launch configurator */
		int minGridSize_prolongate; /**< the minimum grid size needed to achieve the maximum occupancy for a full device launch */
		int gridSize_prolongate; /**< the actual grid size needed, based on input size */
	#endif
};

template<> Fem2D<PetscVector>::Fem2D(Decomposition<PetscVector> *decomposition1, Decomposition<PetscVector> *decomposition2, double fem_reduce);
template<> void Fem2D<PetscVector>::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const;
template<> void Fem2D<PetscVector>::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const;
template<> void Fem2D<PetscVector>::compute_decomposition_reduced();
template<> Fem2D<PetscVector>::ExternalContent * Fem2D<PetscVector>::get_externalcontent() const;

/* cuda kernels cannot be a member of class */
#ifdef USE_CUDA
__global__ void kernel_femhat_reduce_data(double *data1, double *data2, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int left_t2_idx, double diff);
__global__ void kernel_femhat_prolongate_data(double *data1, double *data2, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int left_t2_idx, double diff);
#endif


}
} /* end of namespace */

#endif
