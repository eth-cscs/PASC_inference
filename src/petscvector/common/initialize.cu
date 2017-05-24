#include "external/petscvector/petscvector.cuh"
#include "external/petscvector/common/initialize.h"

namespace pascinference {
namespace common {

__global__ void kernel_warmup(){
}

void cuda_warmup(){
	kernel_warmup<<<1,1>>>();
	gpuErrchk( cudaDeviceSynchronize() );
}

void cuda_barrier(){
	gpuErrchk( cudaDeviceSynchronize() );
}

void cuda_copytoGPU(Vec &x){
	TRYCXX( VecCUDACopyToGPU(x) );	
}

}
} /* end of namespace */
