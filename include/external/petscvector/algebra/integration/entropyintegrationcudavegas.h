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
		int nBlockSize;


	public:

};

template<> EntropyIntegrationCudaVegas<PetscVector>::EntropyIntegrationCudaVegas(EntropyData<PetscVector> *entropydata, double new_eps);
template<> EntropyIntegrationCudaVegas<PetscVector>::~EntropyIntegrationCudaVegas();
template<> void EntropyIntegrationCudaVegas<PetscVector>::compute(double *integrals_out, double *lambda, int Km_max);

}
} /* end of namespace */


#endif
