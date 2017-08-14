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
		int xdim;	/**< dimension of integral on CPU */
		int number_of_moments; /**< number of moments on CPU */
		int number_of_integrals; /**< number of integrals on CPU */
		int *matrix_D_arr;  /**< matrix with powers on CPU */

		#ifdef USE_CUDA
			double *g_lambda;  /**< lagrange multipliers on CUDA */
			int *g_matrix_D_arr;  /**< matrix with powers on CUDA */
		#endif

	public:
		int nBlockSize; /**< number of thread block size */
		int ncall;	/**< number of calls */
		int itmx;	/**< number of max. iterations */
		double acc; /**< precision */

		Timer timerVegasCall;
		Timer timerVegasMove;
		Timer timerVegasFill;
		Timer timerVegasRefine;
		
		#ifdef USE_CUDA
			ExternalContent(int xdim, int number_of_moments, int number_of_integrals, int *matrix_D_arr);		
			~ExternalContent();		
			void cuda_gVegas(double &avgi, double &sd, double &chi2a, double *lambda_arr);
		#endif

		int get_xdim() const {
			return this->xdim;
		}

		int get_number_of_moments() const {
			return this->number_of_moments;
		}

		int get_number_of_integrals() const {
			return this->number_of_integrals;
		}

};

template<> EntropyIntegrationCudaVegas<PetscVector>::EntropyIntegrationCudaVegas(EntropyData<PetscVector> *entropydata, double new_eps);
template<> EntropyIntegrationCudaVegas<PetscVector>::~EntropyIntegrationCudaVegas();
template<> void EntropyIntegrationCudaVegas<PetscVector>::compute(double *integrals_out, double *lambda, int Km_max);

}
} /* end of namespace */


#endif
