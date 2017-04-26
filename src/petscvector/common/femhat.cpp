#include "external/petscvector/common/femhat.h"

namespace pascinference {
namespace common {

template <>
void FemHat<PetscVector>::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const {
	LOG_FUNC_BEGIN

	double *gammak1_arr;
	double *gammak2_arr;

	Vec gamma1_Vec = gamma1->get_vector();
	Vec gamma2_Vec = gamma2->get_vector();

	#ifdef USE_CUDA
		TRYCXX( VecCUDACopyToGPU(gamma1_Vec) );
		TRYCXX( VecCUDACopyToGPU(gamma2_Vec) );
	#endif

	Vec gammak1_Vec;
	Vec gammak2_Vec;

	IS gammak1_is;
	IS gammak2_is;

	/* stuff for getting subvector for local computation */
	IS gammak1_sublocal_is;
	Vec gammak1_sublocal_Vec;

	for(int k=0;k<this->decomposition2->get_K();k++){

		/* get gammak */
		this->decomposition1->createIS_gammaK(&gammak1_is, k);
		this->decomposition2->createIS_gammaK(&gammak2_is, k);

		TRYCXX( VecGetSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecGetSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		/* get local necessary part for local computation */
		TRYCXX( ISCreateStride(PETSC_COMM_WORLD, right_t1_idx - left_t1_idx, left_t1_idx, 1, &gammak1_sublocal_is) );
		TRYCXX( VecGetSubVector(gammak1_Vec, gammak1_sublocal_is, &gammak1_sublocal_Vec) );

		#ifndef USE_CUDA
			/* sequential version */
			TRYCXX( VecGetArray(gammak1_sublocal_Vec,&gammak1_arr) );
			TRYCXX( VecGetArray(gammak2_Vec,&gammak2_arr) );

			int Tbegin2 = this->decomposition2->get_Tbegin();

			//TODO: OpenMP?
			for(int t2=0; t2 < this->decomposition2->get_Tlocal(); t2++){
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
						mysum += phi_value*gammak1_arr[id_counter];
					}
					t1 += 1;
					id_counter += 1;
				}

				/* right part of hat function */
				while(t1 < right_t1){
					phi_value = (t1 - right_t1)/(center_t1 - right_t1);
					if(id_counter < right_t1_idx - left_t1_idx){
						mysum += phi_value*gammak1_arr[id_counter];
					}
					t1 += 1;
					id_counter += 1;
				}

				gammak2_arr[t2] = mysum;
			}
		
			TRYCXX( VecRestoreArray(gammak1_sublocal_Vec,&gammak1_arr) );
			TRYCXX( VecRestoreArray(gammak2_Vec,&gammak2_arr) );
			
		#else
			/* cuda version */
			TRYCXX( VecCUDAGetArrayReadWrite(gammak1_sublocal_Vec,&gammak1_arr) );
			TRYCXX( VecCUDAGetArrayReadWrite(gammak2_Vec,&gammak2_arr) );

			kernel_femhat_reduce_data<<<gridSize_reduce, blockSize_reduce>>>(gammak1_arr, gammak2_arr, this->decomposition1->get_T(), this->decomposition2->get_T(), this->decomposition1->get_Tbegin(), this->decomposition2->get_Tbegin(), this->decomposition1->get_Tlocal(), this->decomposition2->get_Tlocal(), left_t1_idx, right_t1_idx, left_t2_idx, right_t2_idx, diff);
			gpuErrchk( cudaDeviceSynchronize() );
			MPI_Barrier( MPI_COMM_WORLD );

			TRYCXX( VecCUDARestoreArrayReadWrite(gammak1_sublocal_Vec,&gammak1_arr) );
			TRYCXX( VecCUDARestoreArrayReadWrite(gammak2_Vec,&gammak2_arr) );			
		#endif

		/* restore local necessary part for local computation */
		TRYCXX( VecRestoreSubVector(gammak1_Vec, gammak1_sublocal_is, &gammak1_sublocal_Vec) );
		TRYCXX( ISDestroy(&gammak1_sublocal_is) );

		TRYCXX( VecRestoreSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecRestoreSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		TRYCXX( ISDestroy(&gammak1_is) );
		TRYCXX( ISDestroy(&gammak2_is) );

	}

	LOG_FUNC_END
}

template<>
void FemHat<PetscVector>::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const {
	LOG_FUNC_BEGIN

	double *gammak1_arr;
	double *gammak2_arr;

	Vec gamma1_Vec = gamma1->get_vector();
	Vec gamma2_Vec = gamma2->get_vector();

	#ifdef USE_CUDA
		TRYCXX( VecCUDACopyToGPU(gamma1_Vec) );
		TRYCXX( VecCUDACopyToGPU(gamma2_Vec) );
	#endif

	Vec gammak1_Vec;
	Vec gammak2_Vec;
	
	IS gammak1_is;
	IS gammak2_is;

	/* stuff for getting subvector for local computation */
	IS gammak2_sublocal_is;
	Vec gammak2_sublocal_Vec;

	for(int k=0;k<this->decomposition2->get_K();k++){

		/* get gammak */
		this->decomposition1->createIS_gammaK(&gammak1_is, k);
		this->decomposition2->createIS_gammaK(&gammak2_is, k);

		TRYCXX( VecGetSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecGetSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		TRYCXX( ISCreateStride(PETSC_COMM_WORLD, right_t2_idx - left_t2_idx + 1, left_t2_idx, 1, &gammak2_sublocal_is) );
		TRYCXX( VecGetSubVector(gammak2_Vec, gammak2_sublocal_is, &gammak2_sublocal_Vec) );

		#ifndef USE_CUDA
			TRYCXX( VecGetArray(gammak1_Vec,&gammak1_arr) );
			TRYCXX( VecGetArray(gammak2_sublocal_Vec,&gammak2_arr) );

			int Tbegin1 = this->decomposition1->get_Tbegin();

			//TODO: OpenMP?
			for(int t1=0; t1 < this->decomposition1->get_Tlocal(); t1++){
				int t2_left_id_orig = floor((t1 + Tbegin1)/diff);
				int t2_right_id_orig = floor((t1 + Tbegin1)/diff) + 1;

				double t1_left = t2_left_id_orig*diff;
				double t1_right = t2_right_id_orig*diff;

				int t2_left_id = t2_left_id_orig - left_t2_idx;
				int t2_right_id = t2_right_id_orig - left_t2_idx;

				/* value of basis functions */
				double t1_value = 0.0;
				double phi_value_left = (t1 + Tbegin1 - t1_left)/(t1_right - t1_left); 
				t1_value += phi_value_left*gammak2_arr[t2_right_id];
				
				double phi_value_right = (t1 + Tbegin1 - t1_right)/(t1_left - t1_right); 
				t1_value += phi_value_right*gammak2_arr[t2_left_id];

				gammak1_arr[t1] = t1_value;
			}

			TRYCXX( VecRestoreArray(gammak1_Vec,&gammak1_arr) );
			TRYCXX( VecRestoreArray(gammak2_sublocal_Vec,&gammak2_arr) );
		#else
			/* cuda version */
			TRYCXX( VecCUDAGetArrayReadWrite(gammak1_Vec,&gammak1_arr) );
			TRYCXX( VecCUDAGetArrayReadWrite(gammak2_sublocal_Vec,&gammak2_arr) );

			kernel_femhat_prolongate_data<<<gridSize_prolongate, blockSize_prolongate>>>(gammak1_arr, gammak2_arr, this->decomposition1->get_T(), this->decomposition2->get_T(), this->decomposition1->get_Tbegin(), this->decomposition2->get_Tbegin(), this->decomposition1->get_Tlocal(), this->decomposition2->get_Tlocal(), left_t1_idx, right_t1_idx, left_t2_idx, right_t2_idx, diff);
			gpuErrchk( cudaDeviceSynchronize() );
			MPI_Barrier( MPI_COMM_WORLD );

			TRYCXX( VecCUDARestoreArrayReadWrite(gammak1_Vec,&gammak1_arr) );
			TRYCXX( VecCUDARestoreArrayReadWrite(gammak2_sublocal_Vec,&gammak2_arr) );			

		#endif

		/* restore local necessary part for local computation */
		TRYCXX( VecRestoreSubVector(gammak2_Vec, gammak2_sublocal_is, &gammak2_sublocal_Vec) );
		TRYCXX( ISDestroy(&gammak2_sublocal_is) );

		TRYCXX( VecRestoreSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecRestoreSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		TRYCXX( ISDestroy(&gammak1_is) );
		TRYCXX( ISDestroy(&gammak2_is) );
	}

	LOG_FUNC_END
}



#ifdef USE_CUDA
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

#endif


}
} /* end of namespace */

