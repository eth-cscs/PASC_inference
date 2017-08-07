#include "external/petscvector/data/entropydata.h"

namespace pascinference {
namespace data {
	
/* constructor */
template<>
EntropyData<PetscVector>::EntropyData(Decomposition<PetscVector> *decomposition, int Km){
	LOG_FUNC_BEGIN

	/* set initial content */
	this->lambda = NULL;
	this->x = NULL;
	this->gamma = NULL;

	this->decomposition = decomposition;
	this->Km = Km;

	this->number_of_moments = compute_number_of_moments(get_xdim(), this->Km);
	this->matrix_D = new int[get_xdim()*this->number_of_moments];
	this->prepare_matrix_D();

	/* prepare external content */
	externalcontent = new ExternalContent();
	externalcontent->x_powers_Vecs = new Vec[Km+1]; /* 0 = x^0, 1 = x^1 ... Km = x^Km */

	LOG_FUNC_END
}

/* destructor */
template<>
EntropyData<PetscVector>::~EntropyData(){
	LOG_FUNC_BEGIN
	
	free(this->matrix_D);

	/* destroy auxiliary vectors of powers */
	free(externalcontent->x_powers_Vecs);
	//TODO: VecDestroy to x_powers_Vecs?

	/* destroy external content */
	free(externalcontent);
		
	LOG_FUNC_END
}

template<> 
void EntropyData<PetscVector>::set_x(GeneralVector<PetscVector> *x) {
	LOG_FUNC_BEGIN

	this->x = x;

	/* prepare and compute auxiliary array of powers */
	// TODO: aaaah! bottleneck, maybe can be performed in different way (I hope so)
	Vec x_Vec = this->x->get_vector();
	TRYCXX( VecDuplicate(x_Vec, &(externalcontent->x_powers_Vecs[0])) );
	TRYCXX( VecSet(externalcontent->x_powers_Vecs[0], 1.0) );
	TRYCXX( VecAssemblyBegin(externalcontent->x_powers_Vecs[0]));
	TRYCXX( VecAssemblyEnd(externalcontent->x_powers_Vecs[0]));

	for(int km = 1; km <= get_Km();km++){
		TRYCXX( VecDuplicate(x_Vec, &(externalcontent->x_powers_Vecs[km])) );
		TRYCXX( VecPointwiseMult(externalcontent->x_powers_Vecs[km], externalcontent->x_powers_Vecs[km-1], x_Vec) );
		TRYCXX( VecAssemblyBegin(externalcontent->x_powers_Vecs[km]));
		TRYCXX( VecAssemblyEnd(externalcontent->x_powers_Vecs[km]));
	}

	LOG_FUNC_END
}

template<>
void EntropyData<PetscVector>::compute_moments(GeneralVector<PetscVector> *moments) {
	LOG_FUNC_BEGIN

	int Tlocal = get_decomposition()->get_Tlocal();
	int Rlocal = get_decomposition()->get_Rlocal();

	int T = get_decomposition()->get_T();
	int R = get_decomposition()->get_R();

	/* I assume that externalcontent->x_powers_Vecs is computed and constant */
	/* I assume that D_matrix is computed and prepared */

	int xdim = get_xdim();
	int number_of_moments = get_number_of_moments();

	Vec x_Vec = get_x()->get_vector();
	Vec gamma_Vec = get_gamma()->get_vector();
	Vec moments_Vec = moments->get_vector();

	Vec gammak_Vec;
	IS gammak_is;

	/* temp = (x_1^D*x_2^D*...) */
	Vec temp_Vec;
	TRYCXX( VecCreate(PETSC_COMM_WORLD,&temp_Vec) );
	#ifdef USE_CUDA
		TRYCXX(VecSetType(temp_Vec, VECMPICUDA));
	#else
		TRYCXX(VecSetType(temp_Vec, VECMPI));
	#endif
	TRYCXX( VecSetSizes(temp_Vec, Tlocal*Rlocal, T*R) );
	TRYCXX( VecSetFromOptions(temp_Vec) );

	/* temp2 = x_n^D */
	Vec temp2_Vec;
	TRYCXX( VecDuplicate(temp_Vec, &temp2_Vec));

	/* temp3 = gammak*temp; */
	Vec temp3_Vec;
	TRYCXX( VecDuplicate(temp_Vec, &temp3_Vec)); //TODO: I cannot reuse temp2 ?

	/* mom = sum(temp3) */
	IS xn_is;

	double *moments_arr;
	TRYCXX( VecGetArray(moments_Vec, &moments_arr) );

	int *matrix_D_arr = get_matrix_D();

	double mysum, gammaksum;

	int D_value;

	for(int D_row_idx=0; D_row_idx < number_of_moments; D_row_idx++){ /* go through all rows of matrix D */

		/* temp = 1 */
		TRYCXX( VecSet(temp_Vec, 1.0) );
		TRYCXX( VecAssemblyBegin(temp_Vec));
		TRYCXX( VecAssemblyEnd(temp_Vec));

		/* throught columns of D */
		for(int D_col_idx=0; D_col_idx < xdim; D_col_idx++){
			D_value = (int)matrix_D_arr[D_row_idx*xdim + D_col_idx];

			/* get x_n^D */
			get_decomposition()->createIS_datan(&xn_is, D_col_idx);

			TRYCXX( VecGetSubVector(externalcontent->x_powers_Vecs[D_value], xn_is, &temp2_Vec) );

			/* compute temp *= x_n^D */
			TRYCXX( VecPointwiseMult(temp_Vec, temp_Vec, temp2_Vec) );
			TRYCXX( VecAssemblyBegin(temp_Vec));
			TRYCXX( VecAssemblyEnd(temp_Vec));

			TRYCXX( VecRestoreSubVector(externalcontent->x_powers_Vecs[D_value], xn_is, &temp2_Vec) );
			TRYCXX( ISDestroy(&xn_is) );
		}

		/* go throught clusters and multiply with coefficients */
		for(int k=0;k<get_K();k++){
			/* get gammak */
			get_decomposition()->createIS_gammaK(&gammak_is, k);
			TRYCXX( VecGetSubVector(gamma_Vec, gammak_is, &gammak_Vec) );

			/* compute temp_Vec*gammak_Vec */
			TRYCXX( VecPointwiseMult(temp3_Vec, gammak_Vec, temp_Vec) ); /* x_power_gammak = x_power.*gammak */

			/* compute gammaksum */
			TRYCXX( VecSum(gammak_Vec, &gammaksum) );
			TRYCXX( VecSum(temp3_Vec, &mysum) );

			/* store computed moment */
			if(gammaksum != 0){
				moments_arr[k*number_of_moments + D_row_idx] = mysum/gammaksum;
			} else {
				coutMaster << "ERROR: norm(gammak) = 0" << std::endl;

				moments_arr[k*get_number_of_moments() + D_row_idx] = 0.0;
			}

			TRYCXX( VecRestoreSubVector(gamma_Vec, gammak_is, &gammak_Vec) );
			TRYCXX( ISDestroy(&gammak_is) );
		}

	}

	TRYCXX( VecRestoreArray(moments_Vec, &moments_arr) );

	LOG_FUNC_END
}

template<>
void EntropyData<PetscVector>::compute_residuum(GeneralVector<PetscVector> *residuum, GeneralVector<PetscVector> *integrals) const {
	LOG_FUNC_BEGIN

	int T = get_T();
	int Tlocal = get_decomposition()->get_Tlocal();
	int K = get_K();
	int Km = get_Km();
	int xdim = get_xdim();
	int number_of_moments = get_number_of_moments();

	int n = number_of_moments-1;
	int number_of_integrals = 1 + n + (int)(0.5*n*(n+1));

	/* update gamma_solver data - prepare new linear term */
	/* theta includes all moments */
	double *lambda_arr;
	TRYCXX( VecGetArray(get_lambda()->get_vector(), &lambda_arr) );

	double *integrals_arr;
	TRYCXX( VecGetArray(integrals->get_vector(), &integrals_arr) );

	/* mom_powers - exponents */
	int *matrix_D_arr = get_matrix_D();

	/* index set of the data for given dimension component 1,...,xdim */
	IS xn_is;

	/* indexes of appropriate components in residuum 1,...,K */
	IS gammak_is;
	Vec gammak_Vec;

	/* temp = (x_1^D*x_2^D*...) */
	Vec temp_Vec;
	TRYCXX( VecCreate(PETSC_COMM_WORLD,&temp_Vec) );
	#ifdef USE_CUDA
		TRYCXX(VecSetType(temp_Vec, VECMPICUDA));
	#else
		TRYCXX(VecSetType(temp_Vec, VECMPI));
	#endif
	TRYCXX( VecSetSizes(temp_Vec, Tlocal, T) );
	TRYCXX( VecSetFromOptions(temp_Vec) );

	/* temp2 = x_n^D */
	Vec temp2_Vec;
	TRYCXX( VecDuplicate(temp_Vec, &temp2_Vec));

	/* vector of residuum */
	Vec residuum_Vec = residuum->get_vector();
	TRYCXX( VecSet(residuum_Vec, 0.0) ); /* residuum = 0 */
	/* add log part to residuum */
	double logF;
	for(int k=0; k<K;k++){
		logF = log(integrals_arr[k*number_of_integrals]);

		get_decomposition()->createIS_gammaK(&gammak_is, k);
		TRYCXX( VecGetSubVector(residuum_Vec, gammak_is, &gammak_Vec) );

		/* multiply with correspoinding computed lagrange multiplier and add it to residuum */
		TRYCXX( VecSet(gammak_Vec, logF) );

		TRYCXX( VecRestoreSubVector(residuum_Vec, gammak_is, &gammak_Vec) );
		TRYCXX( ISDestroy(&gammak_is) );
	}

	int D_value;
	for(int D_row_idx=1; D_row_idx < number_of_moments; D_row_idx++){ /* go through all rows of matrix D */

		/* temp = 1 */
		TRYCXX( VecSet(temp_Vec, 1.0) );

		/* throught columns of D */
		for(int D_col_idx=0; D_col_idx < xdim; D_col_idx++){
			D_value = (int)matrix_D_arr[D_row_idx*xdim + D_col_idx];

			/* get x_n^D */
			get_decomposition()->createIS_datan(&xn_is, D_col_idx);
			TRYCXX( VecGetSubVector(externalcontent->x_powers_Vecs[D_value], xn_is, &temp2_Vec) );

			/* compute temp *= x_n^D */
			TRYCXX( VecPointwiseMult(temp_Vec, temp_Vec, temp2_Vec) );

			TRYCXX( VecRestoreSubVector(externalcontent->x_powers_Vecs[D_value], xn_is, &temp2_Vec) );
			TRYCXX( ISDestroy(&xn_is) );
		}

		/* go throught clusters and multiply with lambda */
		for(int k=0;k<get_K();k++){
			/* get gammak */
			get_decomposition()->createIS_gammaK(&gammak_is, k);
			TRYCXX( VecGetSubVector(residuum_Vec, gammak_is, &gammak_Vec) );

			/* multiply with correspoinding computed lagrange multiplier and add it to residuum */
			TRYCXX( VecAXPY(gammak_Vec, lambda_arr[k*number_of_moments+D_row_idx], temp_Vec) );

			TRYCXX( VecRestoreSubVector(residuum_Vec, gammak_is, &gammak_Vec) );
			TRYCXX( ISDestroy(&gammak_is) );
		}

	}

	/* restore arrays */
	TRYCXX( VecRestoreArray(get_lambda()->get_vector(), &lambda_arr) );
	TRYCXX( VecRestoreArray(integrals->get_vector(), &integrals_arr) );

	LOG_FUNC_END
}

template<> EntropyData<PetscVector>::ExternalContent * EntropyData<PetscVector>::get_externalcontent() const {
	return externalcontent;
}

}
}
