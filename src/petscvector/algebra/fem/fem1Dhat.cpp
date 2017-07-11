#include "external/petscvector/algebra/fem/fem1Dhat.h"

namespace pascinference {
namespace algebra {

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

		/* sequential version */
		TRYCXX( VecGetArray(gammak1_sublocal_Vec,&gammak1_arr) );
		TRYCXX( VecGetArray(gammak2_Vec,&gammak2_arr) );

		int Tbegin2 = this->decomposition2->get_Tbegin();

		//TODO: OpenMP? GPU?
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

		TRYCXX( VecGetArray(gammak1_Vec,&gammak1_arr) );
		TRYCXX( VecGetArray(gammak2_sublocal_Vec,&gammak2_arr) );

		int Tbegin1 = this->decomposition1->get_Tbegin();

		//TODO: OpenMP? GPU?
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


}
} /* end of namespace */

