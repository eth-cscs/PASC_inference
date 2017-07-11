#include "external/petscvector/algebra/fem/fem1Dsum.h"

namespace pascinference {
namespace algebra {

template<>
Fem1DSum<PetscVector>::Fem1DSum(Decomposition<PetscVector> *decomposition1, Decomposition<PetscVector> *decomposition2, double fem_reduce) : Fem<PetscVector>(decomposition1, decomposition2, fem_reduce){
	LOG_FUNC_BEGIN

	diff = (decomposition1->get_T())/(double)(decomposition2->get_T());
	externalcontent = new ExternalContent();

	LOG_FUNC_END
}

template<>
void Fem1DSum<PetscVector>::reduce_gamma(GeneralVector<PetscVector> *gamma1, GeneralVector<PetscVector> *gamma2) const {
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

	for(int k=0;k<decomposition2->get_K();k++){

		/* get gammak */
		decomposition1->createIS_gammaK(&gammak1_is, k);
		decomposition2->createIS_gammaK(&gammak2_is, k);

		TRYCXX( VecGetSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecGetSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		/* get local necessary part for local computation */
		TRYCXX( ISCreateStride(PETSC_COMM_WORLD, round(decomposition2->get_Tlocal()*diff), round(decomposition2->get_Tbegin()*diff), 1, &gammak1_sublocal_is) );
		TRYCXX( VecGetSubVector(gammak1_Vec, gammak1_sublocal_is, &gammak1_sublocal_Vec) );

		/* sequential version */
		TRYCXX( VecGetArray(gammak1_sublocal_Vec,&gammak1_arr) );
		TRYCXX( VecGetArray(gammak2_Vec,&gammak2_arr) );

		//TODO: OpenMP? GPU?
		for(int t2=0; t2 < decomposition2->get_Tlocal(); t2++){
			double mysum = 0.0;
			for(int i=round(t2*diff); i < round((t2+1)*diff);i++){
				mysum += gammak1_arr[i];
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
void Fem1DSum<PetscVector>::prolongate_gamma(GeneralVector<PetscVector> *gamma2, GeneralVector<PetscVector> *gamma1) const {
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

	for(int k=0;k<decomposition2->get_K();k++){

		/* get gammak */
		decomposition1->createIS_gammaK(&gammak1_is, k);
		decomposition2->createIS_gammaK(&gammak2_is, k);

		TRYCXX( VecGetSubVector(gamma1_Vec, gammak1_is, &gammak1_Vec) );
		TRYCXX( VecGetSubVector(gamma2_Vec, gammak2_is, &gammak2_Vec) );

		/* get local necessary part for local computation */
		TRYCXX( ISCreateStride(PETSC_COMM_WORLD, round(decomposition2->get_Tlocal()*diff), round(decomposition2->get_Tbegin()*diff), 1, &gammak1_sublocal_is) );
		TRYCXX( VecGetSubVector(gammak1_Vec, gammak1_sublocal_is, &gammak1_sublocal_Vec) );

		TRYCXX( VecGetArray(gammak1_sublocal_Vec,&gammak1_arr) );
		TRYCXX( VecGetArray(gammak2_Vec,&gammak2_arr) );

		//TODO: OpenMP? GPU?
		for(int t2=0; t2 < decomposition2->get_Tlocal(); t2++){
			for(int i=round(t2*diff); i < round((t2+1)*diff);i++){
				gammak1_arr[i] = gammak2_arr[t2];
			}
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
void Fem1DSum<PetscVector>::compute_decomposition_reduced() {
	LOG_FUNC_BEGIN
	
	if(is_reduced()){
		int T_reduced = ceil(decomposition1->get_T()*fem_reduce);
		
		/* compute new decomposition */
		decomposition2 = new Decomposition<PetscVector>(T_reduced, 
				*(decomposition1->get_graph()), 
				decomposition1->get_K(), 
				decomposition1->get_xdim(), 
				decomposition1->get_DDT_size(), 
				decomposition1->get_DDR_size());

	} else {
		/* there is not reduction of the data, we can reuse the decomposition */
		this->set_decomposition_reduced(decomposition1);
	}

	diff = (decomposition1->get_T())/(double)(decomposition2->get_T());
	
	LOG_FUNC_END
}

template<> Fem1DSum<PetscVector>::ExternalContent * Fem1DSum<PetscVector>::get_externalcontent() const {
	return externalcontent;
}


}
} /* end of namespace */

