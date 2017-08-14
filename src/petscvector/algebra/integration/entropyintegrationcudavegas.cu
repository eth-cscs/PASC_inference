#include "external/petscvector/petscvector.cuh"
#include "external/petscvector/algebra/integration/entropyintegrationcudavegas.h"

namespace pascinference {
namespace algebra {

EntropyIntegrationCudaVegas<PetscVector>::ExternalContent::ExternalContent() {
	LOG_FUNC_BEGIN

	/* restart timers */
	this->timeVegasCall = 0.0;
	this->timeVegasMove = 0.0;
	this->timeVegasFill = 0.0;
	this->timeVegasRefine = 0.0;

	LOG_FUNC_END
}

void EntropyIntegrationCudaVegas<PetscVector>::ExternalContent::cuda_gVegas(double &avgi, double &sd, double &chi2a) {
	LOG_FUNC_BEGIN

	avgi = 11.1;
	sd = 22.2;
	chi2a = 33.33;

	LOG_FUNC_END
}



}
}

