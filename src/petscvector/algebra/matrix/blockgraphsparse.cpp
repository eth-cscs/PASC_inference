#include "external/petscvector/algebra/matrix/blockgraphsparse.h"

namespace pascinference {
namespace algebra {

template<>
BlockGraphSparseMatrix<PetscVector>::BlockGraphSparseMatrix(Decomposition<PetscVector> &new_decomposition, double alpha, GeneralVector<PetscVector> *new_coeffs){
	LOG_FUNC_BEGIN

	this->decomposition = &new_decomposition;

	this->alpha = alpha;
	this->coeffs = new_coeffs;

	int K = get_K();

	int T = get_T();
	int Tlocal = get_Tlocal();
	int Tbegin = decomposition->get_Tbegin();
	int Tend = decomposition->get_Tend();
	int R = get_R();
	int Rlocal = get_Rlocal();
	int Rbegin = decomposition->get_Rbegin();
	int Rend = decomposition->get_Rend();

	int* neighbor_nmbs = decomposition->get_graph()->get_neighbor_nmbs();
	int **neightbor_ids = decomposition->get_graph()->get_neighbor_ids();

	/* prepare external content with PETSc stuff */
	externalcontent = new ExternalContent();

	/* create matrix */
	TRYCXX( MatCreate(PETSC_COMM_WORLD,&(externalcontent->A_petsc)) );
	TRYCXX( MatSetSizes(externalcontent->A_petsc,K*Rlocal*Tlocal,K*Rlocal*Tlocal,K*R*T,K*R*T) );

	#ifndef USE_CUDA
		TRYCXX( MatSetType(externalcontent->A_petsc,MATMPIAIJ) );
	#else
		TRYCXX( MatSetType(externalcontent->A_petsc,MATAIJCUSPARSE) );
	#endif

	/* compute preallocation of number of non-zero elements in matrix */
	TRYCXX( MatMPIAIJSetPreallocation(externalcontent->A_petsc,3*(1+2*decomposition->get_graph()->get_m_max()),NULL,2*(decomposition->get_graph()->get_m_max()+1),NULL) );
	TRYCXX( MatSeqAIJSetPreallocation(externalcontent->A_petsc,3*(1+2*decomposition->get_graph()->get_m_max()),NULL) );

	TRYCXX( MatSetFromOptions(externalcontent->A_petsc) );

	double coeff = 1.0;

////	#pragma omp parallel for
	for(int k=0; k < K; k++){
		for(int r=Rbegin; r < Rend; r++){
			for(int t=Tbegin;t < Tend;t++){
				int diag_idx = this->decomposition->get_pdTRb_idx(t,this->decomposition->get_DDR_permutation(r),K,k);

				int Wsum;

				/* compute sum of W entries in row */
				if(t == 0 || t == T-1){
					if(T > 1){
						Wsum = 2*neighbor_nmbs[r]+2; /* +1 for diagonal block */
//						Wsum = 2*neighbor_nmbs[r_orig]; /* +1 for diagonal block */
					} else {
						Wsum = neighbor_nmbs[r];
					}
				} else {
					if(T > 1){
						Wsum = 3*neighbor_nmbs[r]+4; /* +2 for diagonal block */
//						Wsum = 3*neighbor_nmbs[r_orig]; /* +2 for diagonal block */
					} else {
						Wsum = (neighbor_nmbs[r])+1; /* +1 for diagonal block */
					}
				}

				/* diagonal entry */
				TRYCXX( MatSetValue(externalcontent->A_petsc, diag_idx, diag_idx, coeff*Wsum, INSERT_VALUES) );

				/* my nondiagonal entries */
				if(T>1){
					if(t > 0) {
						TRYCXX( MatSetValue(externalcontent->A_petsc, diag_idx, this->decomposition->get_pdTRb_idx(t-1,this->decomposition->get_DDR_permutation(r),K,k), -2*coeff, INSERT_VALUES) );
					}
					if(t < T-1) {
						TRYCXX( MatSetValue(externalcontent->A_petsc, diag_idx, this->decomposition->get_pdTRb_idx(t+1,this->decomposition->get_DDR_permutation(r),K,k), -2*coeff, INSERT_VALUES) );
					}
				}

				/* non-diagonal neighbor entries */
				for(int neighbor=0;neighbor<neighbor_nmbs[r];neighbor++){
					int r_neighbor = this->decomposition->get_DDR_permutation(neightbor_ids[r][neighbor]);
					int idx2 = this->decomposition->get_pdTRb_idx(t,r_neighbor,K,k);
					TRYCXX( MatSetValue(externalcontent->A_petsc, diag_idx, idx2, -coeff, INSERT_VALUES) );

					if(t > 0) {
						TRYCXX( MatSetValue(externalcontent->A_petsc, diag_idx, this->decomposition->get_pdTRb_idx(t-1,r_neighbor,K,k), -coeff, INSERT_VALUES) );
					}
					if(t < T-1) {
						TRYCXX( MatSetValue(externalcontent->A_petsc, diag_idx, this->decomposition->get_pdTRb_idx(t+1,r_neighbor,K,k), -coeff, INSERT_VALUES) );
					}
				}
			} /* end T */

		} /* end R */

	} /* end K */

	/* finish all writting in matrix */
	TRYCXX( PetscBarrier(NULL) );

	/* assemble matrix */
	TRYCXX( MatAssemblyBegin(externalcontent->A_petsc,MAT_FINAL_ASSEMBLY) );
	TRYCXX( MatAssemblyEnd(externalcontent->A_petsc,MAT_FINAL_ASSEMBLY) );
	TRYCXX( PetscObjectSetName((PetscObject)(externalcontent->A_petsc),"Regularization matrix") );

	LOG_FUNC_END
}


template<>
BlockGraphSparseMatrix<PetscVector>::~BlockGraphSparseMatrix(){
	LOG_FUNC_BEGIN

	if(petscvector::PETSC_INITIALIZED){ /* maybe Petsc was already finalized and there is nothing to destroy */
		TRYCXX( MatDestroy(&(externalcontent->A_petsc)) );
	}

	LOG_FUNC_END
}

/* print matrix */
template<>
void BlockGraphSparseMatrix<PetscVector>::printcontent(ConsoleOutput &output) const
{
	LOG_FUNC_BEGIN

	output << "BlockGraphMatrix (MatView from Petsc follows):" << std::endl;
	output << "----------------------------------------------------------" << std::endl;

	TRYCXX( MatView(externalcontent->A_petsc, PETSC_VIEWER_STDOUT_WORLD) );

	output << "----------------------------------------------------------" << std::endl;

	LOG_FUNC_END
}

/* matrix-vector multiplication */
template<>
void BlockGraphSparseMatrix<PetscVector>::matmult(PetscVector &y, const PetscVector &x) const {
	LOG_FUNC_BEGIN

	// TODO: maybe y is not initialized, who knows

	/* multiply with constant part of matrix */
	TRYCXX( MatMult(externalcontent->A_petsc, x.get_vector(), y.get_vector()) );

	/* multiply with coeffs */
	if(coeffs){
		int K = decomposition->get_K();

		double *coeffs_arr;
		TRYCXX( VecGetArray(coeffs->get_vector(),&coeffs_arr) );

		Vec xk_Vec;
		Vec x_Vec = y.get_vector();
		IS xk_is;
		double coeff;

		/* get vector corresponding to coeff */
		for(int k=0;k<K;k++){
			this->decomposition->createIS_gammaK(&xk_is, k);
			TRYCXX( VecGetSubVector(x_Vec, xk_is, &xk_Vec) );

//			coeff = coeffs_arr[k]*coeffs_arr[k];
			coeff = alpha*coeffs_arr[k]*coeffs_arr[k];
			TRYCXX( VecScale(xk_Vec, coeff) );

			TRYCXX( VecRestoreSubVector(x_Vec, xk_is, &xk_Vec) );
			TRYCXX( ISDestroy(&xk_is) );
		}

		TRYCXX( VecRestoreArray(coeffs->get_vector(),&coeffs_arr) );
	} else {
	    TRYCXX( VecScale(y.get_vector(), this->alpha) );
	}

	TRYCXX( VecAssemblyBegin(y.get_vector()) );
	TRYCXX( VecAssemblyEnd(y.get_vector()) );

	LOG_FUNC_END
}

template<>
BlockGraphSparseMatrix<PetscVector>::ExternalContent * BlockGraphSparseMatrix<PetscVector>::get_externalcontent() const {
	return this->externalcontent;
}

}
} /* end of namespace */

