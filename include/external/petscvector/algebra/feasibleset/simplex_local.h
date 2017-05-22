#ifndef PASC_PETSCVECTOR_SIMPLEXFEASIBLESET_LOCAL_H
#define	PASC_PETSCVECTOR_SIMPLEXFEASIBLESET_LOCAL_H

#include "external/petscvector/algebra/vector/generalvector.h"
#include "general/algebra/feasibleset/simplex_local.h"

namespace pascinference {
namespace algebra {

/* external-specific stuff */
template<> class SimplexFeasibleSet_Local<PetscVector>::ExternalContent {
	public:

	#ifdef USE_CUDA
		double *x_sorted; /**< for manipulation with sorted data on GPU */
		int blockSize; /**< block size returned by the launch configurator */
		int minGridSize; /**< the minimum grid size needed to achieve the maximum occupancy for a full device launch */
		int gridSize; /**< the actual grid size needed, based on input size */


		/** @brief kernel projection onto simplex
		*
		* computes a projection of a point onto simplex in nD
		*
		* take K-dimensional vector x[tK,tK+1,...tK+(K-1)] =: p
		* and compute projection
		* P(p) = arg min || p - y ||_2
		* subject to constraints (which define simplex)
		* y_0 + ... y_{K-1} = 1
		* y_i >= 0 for all i=0,...,K-1
		*
		* in practical applications K is much more lower number than T
		* K - number of clusters (2 - 10^2)
		* T - length of time-series (10^5 - 10^9) 
		* 
		* @param x values of whole vector in array
		* @param x_sorted allocated vector for manipulating with sorted x
		* @param Tlocal parameter of vector length
		* @param K parameter of vector length
		*/ 
		void cuda_project(Vec &x_Vec, int T, int K);
		
		void cuda_create(int T, int K);
		void cuda_destroy();
	#endif
};


template<> void SimplexFeasibleSet_Local<PetscVector>::project(GeneralVector<PetscVector> &x);

}
} /* end of namespace */


#endif
