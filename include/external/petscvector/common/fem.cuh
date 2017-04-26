#ifndef PASC_FEM_CUDA_H
#define	PASC_FEM_CUDA_H

#include "common/fem.h"

namespace pascinference {
namespace common {

extern void cuda_Fem_cuda_occupancy(){
	LOG_FUNC_BEGIN

	/* compute optimal kernel calls */
	gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_reduce, &blockSize_reduce, kernel_fem_reduce_data, 0, 0) );
	gridSize_reduce = (decomposition2->get_Tlocal() + blockSize_reduce - 1)/ blockSize_reduce;

	gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_prolongate, &blockSize_prolongate, kernel_fem_prolongate_data, 0, 0) );
	gridSize_prolongate = (decomposition2->get_Tlocal() + blockSize_prolongate - 1)/ blockSize_prolongate;

	LOG_FUNC_END
}

/* cuda kernels cannot be a member of class */
__global__ void kernel_fem_reduce_data(double *data1, double *data2, int T1, int T2, int T2local, double diff);
__global__ void kernel_fem_prolongate_data(double *data1, double *data2, int T1, int T2, int T2local, double diff);


}
} /* end of namespace */

#endif
