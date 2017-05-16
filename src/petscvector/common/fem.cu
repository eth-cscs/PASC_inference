#include "external/petscvector/common/fem.h"

namespace pascinference {
namespace common {

void Fem<PetscVector>::ExternalContent::cuda_Fem_cuda_occupancy(){
	LOG_FUNC_BEGIN

	/* compute optimal kernel calls */
	gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_reduce, &blockSize_reduce, kernel_fem_reduce_data, 0, 0) );
	gridSize_reduce = (decomposition2->get_Tlocal() + blockSize_reduce - 1)/ blockSize_reduce;

	gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_prolongate, &blockSize_prolongate, kernel_fem_prolongate_data, 0, 0) );
	gridSize_prolongate = (decomposition2->get_Tlocal() + blockSize_prolongate - 1)/ blockSize_prolongate;

	LOG_FUNC_END
}

__global__ void kernel_Fem_reduce_data(double *data1, double *data2, int T1, int T2, int T2local, double diff) {
	int t2 = blockIdx.x*blockDim.x + threadIdx.x;

	if(t2 < T2local){
		double mysum = 0.0;
		for(int i=round(t2*diff); i < round((t2+1)*diff);i++){
			mysum += data1[i];
		}

		data2[t2] = mysum;
	}
}


__global__ void kernel_Fem_prolongate_data(double *data1, double *data2, int T1, int T2, int T2local, double diff) {
	int t2 = blockIdx.x*blockDim.x + threadIdx.x;

	if(t2 < T2local){
		for(int i=round(t2*diff); i < round((t2+1)*diff);i++){
			data1[i] = data2[t2];
		}
	}
}


}
} /* end of namespace */
