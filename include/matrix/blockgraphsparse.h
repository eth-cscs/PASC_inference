/** @file blockgraphsparse.h
 *  @brief block graph matrix with tridiag blocks used in GRAPHH1FEM model with matrix-vector multiplication implemented as sparse-matrix operation
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_BLOCKGRAPHSPARSEMATRIX_H
#define	PASC_BLOCKGRAPHSPARSEMATRIX_H

#include "pascinference.h"
#include "common/bgmgraph.h"

#ifndef USE_PETSCVECTOR
 #error 'BLOCKGRAPHSPARSEMATRIX is for PETSCVECTOR only, sorry'
#endif

#ifdef USE_CUDA
    #include <../src/vec/vec/impls/seq/seqcuda/cudavecimpl.h>
#endif

namespace pascinference {

typedef petscvector::PetscVector PetscVector;
typedef Mat PetscMatrix;

template<class VectorBase>
class BlockGraphSparseMatrix: public GeneralMatrix<VectorBase> {
	private:
		Decomposition *decomposition;

		double alpha; /**< general matrix multiplicator */
		
		#ifdef USE_PETSCVECTOR
			/* Petsc stuff */ 
			PetscMatrix A_petsc;
		#endif

		GeneralVector<VectorBase> *coeffs; /**< vector of coefficient for each block */

	public:
		BlockGraphSparseMatrix(Decomposition &decomposition, double alpha=1.0, GeneralVector<VectorBase> *new_coeffs=NULL);
		~BlockGraphSparseMatrix(); /* destructor - destroy inner matrix */

		void print(ConsoleOutput &output) const; /* print matrix */
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		
		void printcontent(ConsoleOutput &output) const;
		std::string get_name() const;
		
		void matmult(VectorBase &y, const VectorBase &x) const; /* y = A*x */

		int get_R() const;
		int get_Rlocal() const;
		int get_K() const;
		int get_T() const;
		int get_Tlocal() const;
		double get_alpha() const;

};

/* ---------------- IMPLEMENTATION -------------------- */

template<class VectorBase>
std::string BlockGraphSparseMatrix<VectorBase>::get_name() const {
	return "BlockGraphSparseMatrix";
}


template<class VectorBase>
BlockGraphSparseMatrix<VectorBase>::BlockGraphSparseMatrix(Decomposition &new_decomposition, double alpha, GeneralVector<VectorBase> *new_coeffs){
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

	/* create matrix */
	TRY( MatCreate(PETSC_COMM_WORLD,&A_petsc) );
	TRY( MatSetSizes(A_petsc,K*Rlocal*Tlocal,K*Rlocal*Tlocal,K*R*T,K*R*T) );

	#ifndef USE_GPU
		TRY( MatSetType(A_petsc,MATMPIAIJ) ); 
	#else
		TRY( MatSetType(A_petsc,MATAIJCUSPARSE) ); 
	#endif

	/* compute preallocation of number of non-zero elements in matrix */
	TRY( MatMPIAIJSetPreallocation(A_petsc,3*(1+2*decomposition->get_graph()->get_m_max()),NULL,2*(decomposition->get_graph()->get_m_max()+1),NULL) ); 
	TRY( MatSeqAIJSetPreallocation(A_petsc,3*(1+2*decomposition->get_graph()->get_m_max()),NULL) );

	TRY( MatSetFromOptions(A_petsc) ); 
	
	double coeff = 1.0;
	
////	#pragma omp parallel for
	for(int k=0; k < K; k++){
		for(int r=Rbegin; r < Rend; r++){
			for(int t=Tbegin;t < Tend;t++){
				int diag_idx = t*R*K + r*K + k;
				int r_orig = decomposition->get_invPr(r);
				
				int Wsum;

				/* compute sum of W entries in row */
				if(t == 0 || t == T-1){
					if(T > 1){
						Wsum = 2*neighbor_nmbs[r_orig]+1; /* +1 for diagonal block */
					} else {
						Wsum = neighbor_nmbs[r_orig];
					}
				} else {
					if(T > 1){
						Wsum = 3*neighbor_nmbs[r_orig]+2; /* +2 for diagonal block */
					} else {
						Wsum = (neighbor_nmbs[r_orig])+1; /* +1 for diagonal block */
					}
				}
				
				/* diagonal entry */
				TRY( MatSetValue(A_petsc, diag_idx, diag_idx, coeff*Wsum, INSERT_VALUES) );

				/* my nondiagonal entries */
				if(T>1){
					if(t > 0) {
						/* t - 1 */
						TRY( MatSetValue(A_petsc, diag_idx, diag_idx-R*K, -coeff, INSERT_VALUES) );
					}
					if(t < T-1) {
						/* t + 1 */
						TRY( MatSetValue(A_petsc, diag_idx, diag_idx+R*K, -coeff, INSERT_VALUES) );
					}
				}

				/* non-diagonal neighbor entries */
				for(int neighbor=0;neighbor<neighbor_nmbs[r_orig];neighbor++){
					int r_new = decomposition->get_Pr(neightbor_ids[r_orig][neighbor]);
					int idx2 = t*R*K + r_new*K + k;

					TRY( MatSetValue(A_petsc, diag_idx, idx2, -coeff, INSERT_VALUES) );
//					TRY( MatSetValue(A_petsc, idx2, diag_idx, -coeff, INSERT_VALUES) );
					if(t > 0) {
						TRY( MatSetValue(A_petsc, diag_idx, idx2-R*K, -coeff, INSERT_VALUES) );
					}
					if(t < T-1) {
						TRY( MatSetValue(A_petsc, diag_idx, idx2+R*K, -coeff, INSERT_VALUES) );
					}
				}
			} /* end T */

		} /* end R */

	} /* end K */

	/* finish all writting in matrix */
	TRY( PetscBarrier(NULL) );

	/* assemble matrix */
	TRY( MatAssemblyBegin(A_petsc,MAT_FINAL_ASSEMBLY) );
	TRY( MatAssemblyEnd(A_petsc,MAT_FINAL_ASSEMBLY) );

	LOG_FUNC_END
}	


template<class VectorBase>
BlockGraphSparseMatrix<VectorBase>::~BlockGraphSparseMatrix(){
	LOG_FUNC_BEGIN
	
	if(petscvector::PETSC_INITIALIZED){ /* maybe Petsc was already finalized and there is nothing to destroy */
		TRY( MatDestroy(&A_petsc) );
	}
	
	LOG_FUNC_END	
}

/* print matrix */
template<class VectorBase>
void BlockGraphSparseMatrix<VectorBase>::print(ConsoleOutput &output) const		
{
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	output << " - T:     " << get_T() << std::endl;
	output << " - R:     " << get_R() << std::endl;
	output << " - K:     " << get_K() << std::endl;
	output << " - size:  " << get_T()*get_R()*get_K() << std::endl;
	output << " - alpha: " << alpha << std::endl;

	if(coeffs){
		output << " - coeffs: " << *coeffs << std::endl;
	}
	
	LOG_FUNC_END	
}

template<class VectorBase>
void BlockGraphSparseMatrix<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global << this->get_name() << std::endl;
	output_global << " - T:     " << get_T() << std::endl;
	output_global.push();
		output_local << " - Tlocal:  " << get_Tlocal() << " (" << decomposition->get_Tbegin() << "," << decomposition->get_Tend() << ")" << std::endl;
		output_local.synchronize();
	output_global.pop();
	
	output_global << " - R:     " << get_R() << std::endl;
	output_global.push();
		output_local << " - Rlocal:  " << get_Rlocal() << " (" << decomposition->get_Rbegin() << "," << decomposition->get_Rend() << ")" << std::endl;
		output_local.synchronize();
	output_global.pop();

	output_global << " - K:     " << get_K() << std::endl;
	output_global << " - size:  " << get_T()*get_R()*get_K() << std::endl;
	output_global.push();
		output_local << " - sizelocal:  " << get_Tlocal()*get_Rlocal()*get_K() << std::endl;
		output_local.synchronize();
	output_global.pop();

	output_global << " - alpha: " << alpha << std::endl;

	if(coeffs){
		output_local << " - coeffs: " << *coeffs << std::endl;
		output_local.synchronize();
	}

	LOG_FUNC_END	
}

/* print matrix */
template<class VectorBase>
void BlockGraphSparseMatrix<VectorBase>::printcontent(ConsoleOutput &output) const		
{
	LOG_FUNC_BEGIN
	
	output << "BlockGraphMatrix (sorry, 'only' MatView from Petsc follows):" << std::endl;
	output << "----------------------------------------------------------" << std::endl;
	
	TRY( MatView(A_petsc, PETSC_VIEWER_STDOUT_WORLD) );

	output << "----------------------------------------------------------" << std::endl;
	
	LOG_FUNC_END
}

/* matrix-vector multiplication */
template<class VectorBase>
void BlockGraphSparseMatrix<VectorBase>::matmult(VectorBase &y, const VectorBase &x) const { 
	LOG_FUNC_BEGIN
	
	// TODO: maybe y is not initialized, who knows

	/* multiply with constant part of matrix */
	TRY( MatMult(A_petsc, x.get_vector(), y.get_vector()) );

	/* multiply with coeffs */
	if(coeffs){
		int K = decomposition->get_K();
		
		double *coeffs_arr;
		TRY( VecGetArray(coeffs->get_vector(),&coeffs_arr) );

		Vec xk_Vec;
		Vec x_Vec = y.get_vector();
		IS xk_is;
		double coeff;
		
		/* get vector corresponding to coeff */
		for(int k=0;k<K;k++){
			this->decomposition->createIS_gammaK(&xk_is, k);
			TRY( VecGetSubVector(x_Vec, xk_is, &xk_Vec) );
			
//			coeff = coeffs_arr[k]*coeffs_arr[k];
			coeff = alpha*coeffs_arr[k]*coeffs_arr[k];
			TRY( VecScale(xk_Vec, coeff) );

			TRY( VecRestoreSubVector(x_Vec, xk_is, &xk_Vec) );
			TRY( ISDestroy(&xk_is) );
		}

		TRY( VecRestoreArray(coeffs->get_vector(),&coeffs_arr) );
	} else {
	    TRY( VecScale(y.get_vector(), this->alpha) );
	}

	TRY( VecAssemblyBegin(y.get_vector()) );
	TRY( VecAssemblyEnd(y.get_vector()) );

	LOG_FUNC_END
}

template<class VectorBase>
int BlockGraphSparseMatrix<VectorBase>::get_R() const { 
	return decomposition->get_R();
}

template<class VectorBase>
int BlockGraphSparseMatrix<VectorBase>::get_K() const { 
	return decomposition->get_K();
}

template<class VectorBase>
int BlockGraphSparseMatrix<VectorBase>::get_T() const { 
	return decomposition->get_T();
}

template<class VectorBase>
int BlockGraphSparseMatrix<VectorBase>::get_Tlocal() const { 
	return decomposition->get_Tlocal();
}

template<class VectorBase>
int BlockGraphSparseMatrix<VectorBase>::get_Rlocal() const { 
	return decomposition->get_Rlocal();
}

template<class VectorBase>
double BlockGraphSparseMatrix<VectorBase>::get_alpha() const { 
	return this->alpha;
}


} /* end of namespace */


#endif
