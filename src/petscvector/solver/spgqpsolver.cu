#include "external/petscvector/petscvector.cuh"
#include "external/petscvector/solver/spgqpsolver.h"

namespace pascinference {
namespace solver {

void SPGQPSolver<PetscVector>::ExternalContent::cuda_copytogpu(Vec &x) const {
	LOG_FUNC_BEGIN

	TRYCXX( VecCUDACopyToGPU(x) );

	LOG_FUNC_END
}


}
}
