#include "external/petscvector/algebra/fem/fem1Dhat.h"

namespace pascinference {
namespace common {

template<>
Fem1DHat<PetscVector>::Fem1DHat(Decomposition<PetscVector> *decomposition1, Decomposition<PetscVector> *decomposition2, double fem_reduce) : Fem<PetscVector>(decomposition1, decomposition2, fem_reduce){
	LOG_FUNC_BEGIN

	#ifdef USE_CUDA
		/* compute optimal kernel calls */
		externalcontent->cuda_occupancy();
		externalcontent->gridSize_reduce = (decomposition2->get_Tlocal() + externalcontent->blockSize_reduce - 1)/ externalcontent->blockSize_reduce;
		externalcontent->gridSize_prolongate = (decomposition1->get_Tlocal() + externalcontent->blockSize_prolongate - 1)/ externalcontent->blockSize_prolongate;
	#endif

	this->diff = (decomposition1->get_T() - 1)/(double)(decomposition2->get_T() - 1);

	compute_overlaps();

	LOG_FUNC_END
}



template<>
void Fem1DHat<PetscVector>::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const {
	LOG_FUNC_BEGIN

	double *gammak1_arr;
	double *gammak2_arr;

	Vec gamma1_Vec = gamma1->get_vector();
	Vec gamma2_Vec = gamma2->get_vector();

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

			externalcontent->cuda_reduce_data(gammak1_arr, gammak2_arr, this->decomposition1->get_T(), this->decomposition2->get_T(), this->decomposition1->get_Tbegin(), this->decomposition2->get_Tbegin(), this->decomposition1->get_Tlocal(), this->decomposition2->get_Tlocal(), left_t1_idx, right_t1_idx, left_t2_idx, right_t2_idx, diff);

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
void Fem1DHat<PetscVector>::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const {
	LOG_FUNC_BEGIN

	double *gammak1_arr;
	double *gammak2_arr;

	Vec gamma1_Vec = gamma1->get_vector();
	Vec gamma2_Vec = gamma2->get_vector();

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
			
			externalcontent->cuda_prolongate_data(gammak1_arr, gammak2_arr, this->decomposition1->get_T(), this->decomposition2->get_T(), this->decomposition1->get_Tbegin(), this->decomposition2->get_Tbegin(), this->decomposition1->get_Tlocal(), this->decomposition2->get_Tlocal(), left_t1_idx, right_t1_idx, left_t2_idx, right_t2_idx, diff);

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

template<>
void Fem1DHat<PetscVector>::compute_decomposition_reduced() {
	LOG_FUNC_BEGIN
	
	if(this->is_reduced()){
		int T_reduced = ceil(this->decomposition1->get_T()*this->fem_reduce);
		
		/* compute new decomposition */
		this->decomposition2 = new Decomposition<PetscVector>(T_reduced, 
				*(this->decomposition1->get_graph()), 
				this->decomposition1->get_K(), 
				this->decomposition1->get_xdim(), 
				this->decomposition1->get_DDT_size(), 
				this->decomposition1->get_DDR_size());

	} else {
		/* there is not reduction of the data, we can reuse the decomposition */
		this->decomposition2 = this->decomposition1;
	}

	#ifdef USE_CUDA
		/* compute optimal kernel calls */
		externalcontent->cuda_occupancy();
		externalcontent->gridSize_reduce = (this->decomposition2->get_Tlocal() + externalcontent->blockSize_reduce - 1)/ externalcontent->blockSize_reduce;
		externalcontent->gridSize_prolongate = (this->decomposition1->get_Tlocal() + externalcontent->blockSize_prolongate - 1)/ externalcontent->blockSize_prolongate;
	#endif

	this->diff = (this->decomposition1->get_T() - 1)/(double)(this->decomposition2->get_T() - 1);

	compute_overlaps();
	
	LOG_FUNC_END
}

}
} /* end of namespace */

