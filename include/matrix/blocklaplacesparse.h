#ifndef BLOCKLAPLACESPARSEMATRIX_H
#define	BLOCKLAPLACESPARSEMATRIX_H

#include "pascinference.h"

#ifndef USE_PETSCVECTOR
 #error 'BLOCKLAPLACESPARSEMATRIX is for PETSCVECTOR only, sorry'
#endif

#ifdef USE_CUDA
    #include <../src/vec/vec/impls/seq/seqcuda/cudavecimpl.h>
#endif

typedef petscvector::PetscVector PetscVector;
typedef Mat PetscMatrix;

namespace pascinference {

template<class VectorBase>
class BlockLaplaceSparseMatrix: public GeneralMatrix<VectorBase> {
	private:
		Decomposition *decomposition;
		
		double alpha; /**< general matrix multiplicator */
		
		#ifdef USE_PETSCVECTOR
			/* Petsc stuff */ 
			PetscMatrix A_petsc;
		#endif
		
	public:
		BlockLaplaceSparseMatrix(Decomposition &decomposition, double alpha=1.0);
		~BlockLaplaceSparseMatrix(); /* destructor - destroy inner matrix */

		void print(ConsoleOutput &output) const; /* print matrix */
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		
		void printcontent(ConsoleOutput &output) const;
		std::string get_name() const;
		
		void matmult(VectorBase &y, const VectorBase &x) const; /* y = A*x */

		int get_K() const;
		int get_T() const;
		int get_Tlocal() const;
		double get_alpha() const;

};

/* --------------- IMPLEMENTATION -------------------- */
// TODO: move to implementation
template<class VectorBase>
std::string BlockLaplaceSparseMatrix<VectorBase>::get_name() const {
	return "BlockLaplaceSparseMatrix";
}

template<>
BlockLaplaceSparseMatrix<PetscVector>::BlockLaplaceSparseMatrix(Decomposition &new_decomposition, double alpha){
	LOG_FUNC_BEGIN

	this->decomposition = &new_decomposition;

	this->alpha = alpha;

	int K = get_K();
	
	int T = get_T();
	int Tlocal = get_Tlocal();
	int Tbegin = decomposition->get_Tbegin();
	int Tend = decomposition->get_Tend();

	/* create matrix */
	TRY( MatCreate(PETSC_COMM_WORLD,&A_petsc) );
	TRY( MatSetSizes(A_petsc,K*Tlocal,K*Tlocal,K*T,K*T) );

	#ifndef USE_CUDA
		TRY( MatSetType(A_petsc,MATMPIAIJ) ); 
	#else
		TRY( MatSetType(A_petsc,MATAIJCUSPARSE) ); 
	#endif

	/* compute preallocation of number of non-zero elements in matrix */
	TRY( MatMPIAIJSetPreallocation(A_petsc,3,NULL,2,NULL) ); 
	TRY( MatSeqAIJSetPreallocation(A_petsc,3,NULL) );

	TRY( MatSetFromOptions(A_petsc) ); 
	
////	#pragma omp parallel for
	for(int k=0; k < K; k++){
		for(int t=Tbegin; t < Tend; t++){
			int diag_idx = t*K + k;

			int Wsum;

			/* compute sum of W entries in row */
			if(T > 1){
				if(t == 0 || t == T-1){
					Wsum = 1;
				} else {
					Wsum = 2;
				}
			} else {
				Wsum = 1;
			}
				
			/* diagonal entry */
			TRY( MatSetValue(A_petsc, diag_idx, diag_idx, Wsum, INSERT_VALUES) );

			/* my nondiagonal entries */
			if(T>1){
				if(t > 0) {
					/* t - 1 */
					TRY( MatSetValue(A_petsc, diag_idx, diag_idx-K, -1.0, INSERT_VALUES) );
				}
				if(t < T-1) {
					/* t + 1 */
					TRY( MatSetValue(A_petsc, diag_idx, diag_idx+K, -1.0, INSERT_VALUES) );
				}
			}

		} /* end T */

	} /* end K */

	/* finish all writting in matrix */
	TRY( PetscBarrier(NULL) );

	/* assemble matrix */
	TRY( MatAssemblyBegin(A_petsc,MAT_FINAL_ASSEMBLY) );
	TRY( MatAssemblyEnd(A_petsc,MAT_FINAL_ASSEMBLY) );

	LOG_FUNC_END
}	

template<class VectorBase>
BlockLaplaceSparseMatrix<VectorBase>::~BlockLaplaceSparseMatrix(){
	LOG_FUNC_BEGIN
	
	if(petscvector::PETSC_INITIALIZED){ /* maybe Petsc was already finalized and there is nothing to destroy */
		TRY( MatDestroy(&A_petsc) );
	}
	
	LOG_FUNC_END	
}

/* print matrix */
template<class VectorBase>
void BlockLaplaceSparseMatrix<VectorBase>::print(ConsoleOutput &output) const		
{
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	output << " - T:     " << get_T() << std::endl;
	output << " - K:     " << get_K() << std::endl;
	output << " - size:  " << get_T()*get_K() << std::endl;

	output << " - alpha: " << alpha << std::endl;

	LOG_FUNC_END	
}

template<class VectorBase>
void BlockLaplaceSparseMatrix<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global << this->get_name() << std::endl;
	output_global << " - T:     " << get_T() << std::endl;
	output_global.push();
		output_local << " - Tlocal:  " << get_Tlocal() << " (" << decomposition->get_Tbegin() << "," << decomposition->get_Tend() << ")" << std::endl;
		output_local.synchronize();
	output_global.pop();
	
	output_global << " - K:     " << get_K() << std::endl;
	output_global << " - size:  " << get_T()*get_K() << std::endl;
	output_global.push();
		output_local << " - sizelocal:  " << get_Tlocal()*get_K() << std::endl;
		output_local.synchronize();
	output_global.pop();

	output_global << " - alpha: " << alpha << std::endl;

	LOG_FUNC_END	
}

/* print matrix */
template<class VectorBase>
void BlockLaplaceSparseMatrix<VectorBase>::printcontent(ConsoleOutput &output) const		
{
	LOG_FUNC_BEGIN

	output << "BlockLaplaceMatrix (sorry, 'only' MatView from Petsc follows):" << std::endl;
	output << "----------------------------------------------------------" << std::endl;
	
	TRY( MatView(A_petsc, PETSC_VIEWER_STDOUT_WORLD) );

	output << "----------------------------------------------------------" << std::endl;	
	
	LOG_FUNC_END
}

template<class VectorBase>
int BlockLaplaceSparseMatrix<VectorBase>::get_K() const { 
	return decomposition->get_K();
}

template<class VectorBase>
int BlockLaplaceSparseMatrix<VectorBase>::get_T() const { 
	return decomposition->get_T();
}

template<class VectorBase>
int BlockLaplaceSparseMatrix<VectorBase>::get_Tlocal() const { 
	return decomposition->get_Tlocal();
}

template<class VectorBase>
double BlockLaplaceSparseMatrix<VectorBase>::get_alpha() const { 
	return this->alpha;
}

/* matrix-vector multiplication */
template<class VectorBase>
void BlockLaplaceSparseMatrix<VectorBase>::matmult(VectorBase &y, const VectorBase &x) const { 
	LOG_FUNC_BEGIN
	
	// TODO: maybe y is not initialized, who knows

	/* multiply with constant part of matrix */
	TRY( MatMult(A_petsc, x.get_vector(), y.get_vector()) );
	TRY( VecScale(y.get_vector(),this->alpha) );

	LOG_FUNC_END	
}


} /* end of namespace */


#endif
