#include "external/petscvector/petscvector.cuh"
#include "external/petscvector/common/femhat.h"

namespace pascinference {
namespace common {

__global__ void kernel_femhat_reduce_data(double *data1, double *data2, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int right_t1_idx, int left_t2_idx, int right_t2_idx, double diff);
__global__ void kernel_femhat_prolongate_data(double *data1, double *data2, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int right_t1_idx, int left_t2_idx, int right_t2_idx, double diff);



void FemHat<PetscVector>::ExternalContent::cuda_occupancy(){
	LOG_FUNC_BEGIN

	/* compute optimal kernel calls */
	gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_reduce, &blockSize_reduce, kernel_femhat_reduce_data, 0, 0) );
	gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize_prolongate, &blockSize_prolongate, kernel_femhat_prolongate_data, 0, 0) );

	LOG_FUNC_END
}

void FemHat<PetscVector>::ExternalContent::cuda_reduce_data(double *data1_arr, double *data2_arr, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int right_t1_idx, int left_t2_idx, int right_t2_idx, double diff){
	LOG_FUNC_BEGIN

	kernel_femhat_reduce_data<<<gridSize_reduce, blockSize_reduce>>>(data1_arr, data2_arr, T1, T2, Tbegin1, Tbegin2, T1local, T2local, left_t1_idx, right_t1_idx, left_t2_idx, right_t2_idx, diff);
	gpuErrchk( cudaDeviceSynchronize() );

	LOG_FUNC_END
}

void FemHat<PetscVector>::ExternalContent::cuda_prolongate_data(double *data1_arr, double *data2_arr, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int right_t1_idx, int left_t2_idx, int right_t2_idx, double diff){
	LOG_FUNC_BEGIN

	kernel_femhat_prolongate_data<<<gridSize_prolongate, blockSize_prolongate>>>(data1_arr, data2_arr, T1, T2, Tbegin1, Tbegin2, T1local, T2local, left_t1_idx, right_t1_idx, left_t2_idx, right_t2_idx, diff);
	gpuErrchk( cudaDeviceSynchronize() );

	LOG_FUNC_END
}


__global__ void kernel_femhat_reduce_data(double *data1, double *data2, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int right_t1_idx, int left_t2_idx, int right_t2_idx, double diff) {
	int t2 = blockIdx.x*blockDim.x + threadIdx.x;

	if(t2 < T2local){
		double center_t1 = (Tbegin2+t2)*diff;
		double left_t1 = (Tbegin2+t2-1)*diff;
		double right_t1 = (Tbegin2+t2+1)*diff;
				
		int id_counter = floor(left_t1) - left_t1_idx; /* first index in provided local t1 array */

		double phi_value; /* value of basis function */

		/* left part of hat function */
		double mysum = 0.0;
		int t1 = floor(left_t1);

		/* compute linear combination with coefficients given by basis functions */
		while(t1 <= center_t1){
			phi_value = (t1 - left_t1)/(center_t1 - left_t1);
			if(id_counter >= 0){
				mysum += phi_value*data1[id_counter];
			}
			t1 += 1;
			id_counter += 1;
		}

		/* right part of hat function */
		while(t1 < right_t1){
			phi_value = (t1 - right_t1)/(center_t1 - right_t1);
			if(id_counter < right_t1_idx - left_t1_idx){
				mysum += phi_value*data1[id_counter];
			}
			t1 += 1;
			id_counter += 1;
		}

		data2[t2] = mysum;
	}
}


__global__ void kernel_femhat_prolongate_data(double *data1, double *data2, int T1, int T2, int Tbegin1, int Tbegin2, int T1local, int T2local, int left_t1_idx, int right_t1_idx, int left_t2_idx, int right_t2_idx, double diff) {
	int t1 = blockIdx.x*blockDim.x + threadIdx.x;

	if(t1 < T1local){
		int t2_left_id_orig = floor((t1 + Tbegin1)/diff);
		int t2_right_id_orig = floor((t1 + Tbegin1)/diff) + 1;

		double t1_left = t2_left_id_orig*diff;
		double t1_right = t2_right_id_orig*diff;

		int t2_left_id = t2_left_id_orig - left_t2_idx;
		int t2_right_id = t2_right_id_orig - left_t2_idx;

		/* value of basis functions */
		double t1_value = 0.0;
		double phi_value_left = (t1 + Tbegin1 - t1_left)/(t1_right - t1_left); 
		t1_value += phi_value_left*data2[t2_right_id];
				
		double phi_value_right = (t1 + Tbegin1 - t1_right)/(t1_left - t1_right); 
		t1_value += phi_value_right*data2[t2_left_id];

		data1[t1] = t1_value;
	}
}


}
} /* end of namespace */

