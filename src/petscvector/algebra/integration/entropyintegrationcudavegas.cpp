#include "external/petscvector/algebra/integration/entropyintegrationcudavegas.h"

namespace pascinference {
namespace algebra {

template<>
void EntropyIntegrationCudaVegas<PetscVector>::compute(double *integrals_arr, double *lambda_arr, int Km_max) {
	LOG_FUNC_BEGIN

	/* call appropriate cuda kernel for integral computation */
	
	LOG_FUNC_END
}



}
}

