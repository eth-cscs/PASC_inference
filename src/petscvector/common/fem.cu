#include "external/petscvector/common/fem.h"

namespace pascinference {
namespace common {

/* cuda kernels cannot be a member of class */
__global__ void kernel_fem_reduce_data(double *data1, double *data2, int T1, int T2, int T2local, double diff);
__global__ void kernel_fem_prolongate_data(double *data1, double *data2, int T1, int T2, int T2local, double diff);

void Fem<PetscVector>::ExternalContent::cuda_occupancy(){
	LOG_FUNC_BEGIN

	/* compute optimal kernel calls */
	gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_reduce, &blockSize_reduce, kernel_fem_reduce_data, 0, 0) );
	gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_prolongate, &blockSize_prolongate, kernel_fem_prolongate_data, 0, 0) );

	LOG_FUNC_END
}

void Fem<PetscVector>::ExternalContent::cuda_reduce_data(Vec &data1_Vec, Vec &data2_Vec, int T1, int T2, int T2local, double diff){
	LOG_FUNC_BEGIN

	double *data1_arr;
	double *data2_arr;

	/* cuda version */
	TRYCXX( VecCUDAGetArrayReadWrite(data1_Vec,&data1_arr) );
	TRYCXX( VecCUDAGetArrayReadWrite(data2_Vec,&data2_arr) );

	kernel_fem_reduce_data<<<gridSize_reduce, blockSize_reduce>>>(data1_arr, data2_arr, T1, T2, T2local, diff);
	gpuErrchk( cudaDeviceSynchronize() );
	MPI_Barrier( MPI_COMM_WORLD );

	TRYCXX( VecCUDARestoreArrayReadWrite(data1_Vec,&data1_arr) );
	TRYCXX( VecCUDARestoreArrayReadWrite(data2_Vec,&data2_arr) );

	LOG_FUNC_END
}

void Fem<PetscVector>::ExternalContent::cuda_prolongate_data(Vec &data1_Vec, Vec &data2_Vec, int T1, int T2, int T2local, double diff){
	LOG_FUNC_BEGIN

	double *data1_arr;
	double *data2_arr;

	/* cuda version */
	TRYCXX( VecCUDAGetArrayReadWrite(data1_Vec,&data1_arr) );
	TRYCXX( VecCUDAGetArrayReadWrite(data2_Vec,&data2_arr) );

	kernel_fem_prolongate_data<<<gridSize_prolongate, blockSize_prolongate>>>(data1_arr, data2_arr, T1, T2, T2local, diff);
	gpuErrchk( cudaDeviceSynchronize() );
	MPI_Barrier( MPI_COMM_WORLD );

	TRYCXX( VecCUDARestoreArrayReadWrite(data1_Vec,&data1_arr) );
	TRYCXX( VecCUDARestoreArrayReadWrite(data2_Vec,&data2_arr) );

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
