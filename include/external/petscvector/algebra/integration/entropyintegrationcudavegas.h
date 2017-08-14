#ifndef PASC_PETSCVECTOR_ENTROPYINTEGRATIONCUDAVEGAS_H
#define	PASC_PETSCVECTOR_ENTROPYINTEGRATIONCUDAVEGAS_H

#include "general/algebra/integration/entropyintegrationcudavegas.h"
#include "external/petscvector/algebra/vector/generalvector.h"
#include "external/petscvector/algebra/integration/entropyintegration.h"



namespace pascinference {
namespace algebra {

/* external-specific stuff */
template<> class EntropyIntegrationCudaVegas<PetscVector>::ExternalContent {
	private:
		#ifdef USE_CUDA
			__device__ __constant__ int g_ndim;	/**< dimension of integral on CUDA */
			__device__ __constant__ int g_number_of_moments; /**< number of moments on CUDA */
			__device__ __constant__ int g_number_of_integrals; /**< number of integrals on CUDA */
			__device__ __constant__ double *g_lambda;  /**< lagrange multipliers on CUDA */
			__device__ __constant__ int *g_matrix_D_arr;  /**< matrix with powers on CUDA */

		#endif

	public:
		int nBlockSize; /**< number of thread block size */
		int ncall;	/**< number of calls */
		int itmx;	/**< number of max. iterations */
		double acc; /**< precision */
		int ndim; /**< dimension of integral */

		Timer timerVegasCall;
		Timer timerVegasMove;
		Timer timerVegasFill;
		Timer timerVegasRefine;
		
		#ifdef USE_CUDA
			ExternalContent();		
			void cuda_gVegas(double &avgi, double &sd, double &chi2a);
		#endif

};

template<> EntropyIntegrationCudaVegas<PetscVector>::EntropyIntegrationCudaVegas(EntropyData<PetscVector> *entropydata, double new_eps);
template<> EntropyIntegrationCudaVegas<PetscVector>::~EntropyIntegrationCudaVegas();
template<> void EntropyIntegrationCudaVegas<PetscVector>::compute(double *integrals_out, double *lambda, int Km_max);

}
} /* end of namespace */


#endif
