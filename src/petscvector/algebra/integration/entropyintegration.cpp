#include "external/petscvector/algebra/integration/entropyintegration.h"

namespace pascinference {
namespace algebra {

template<>
void EntropyIntegration<PetscVector>::compute(GeneralVector<PetscVector> &integrals) {
	LOG_FUNC_BEGIN

	int K = this->entropydata->get_decomposition()->get_K();
	int number_of_moments = this->entropydata->get_number_of_moments();
	int n = number_of_moments-1;
	int number_of_integrals = this->get_number_of_integrals();

	Vec lambda_Vec = this->entropydata->get_lambda()->get_vector();
	Vec integrals_Vec = integrals.get_vector();
	
	IS lambdak_is;
	IS integralsk_is;
	Vec lambdak_Vec;
	Vec integralsk_Vec;
	double *lambdak_arr;
	double *integralsk_arr;

	/* through all clusters */
	for(int k = 0; k < K; k++){

		/* prepare index set to get subvectors from moments, x, g, s, y */
		TRYCXX( ISCreateStride(PETSC_COMM_SELF, n, k*number_of_moments+1, 1, &lambdak_is) );
		TRYCXX( ISCreateStride(PETSC_COMM_SELF, number_of_integrals, k*number_of_integrals, 1, &integralsk_is) );

		/* get subvectors for this cluster */
		TRYCXX( VecGetSubVector(lambda_Vec, lambdak_is, &lambdak_Vec) );
		TRYCXX( VecGetSubVector(integrals_Vec, integralsk_is, &integralsk_Vec) );

		/* compute integrals */
		TRYCXX( VecGetArray(lambdak_Vec,&lambdak_arr));
		TRYCXX( VecGetArray(integralsk_Vec,&integralsk_arr));
		
 		 compute(integralsk_arr, lambdak_arr, number_of_integrals);

		TRYCXX( VecRestoreArray(lambdak_Vec,&lambdak_arr));
		TRYCXX( VecRestoreArray(integralsk_Vec,&integralsk_arr));

		/* restore subvectors */
		TRYCXX( VecRestoreSubVector(lambda_Vec, lambdak_is, &lambdak_Vec) );
		TRYCXX( VecRestoreSubVector(integrals_Vec, integralsk_is, &integralsk_Vec) );

		/* destroy index set */
		TRYCXX( ISDestroy(&lambdak_is) );
		TRYCXX( ISDestroy(&integralsk_is) );

	} /* endfor through clusters */
	
	LOG_FUNC_END
}

}
}

