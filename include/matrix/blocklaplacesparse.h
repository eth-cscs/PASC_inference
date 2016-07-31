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
		int T; /**< dimension of each block */
		int Tlocal; /**< local dimension of each block */
		int Tbegin; /**< ownership begin */
		int Tend; /**< ownership end */
		
		int K; /**< number of diagonal blocks */
		double alpha; /**< general matrix multiplicator */
		
		#ifdef USE_PETSCVECTOR
			/* Petsc stuff */ 
			PetscMatrix A_petsc;
		#endif
		
	public:
		BlockLaplaceSparseMatrix(VectorBase &x, int K, double alpha=1.0);
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
BlockLaplaceSparseMatrix<PetscVector>::BlockLaplaceSparseMatrix(PetscVector &x, int K, double alpha){
	LOG_FUNC_BEGIN

	this->K = K;
	this->alpha = alpha;
	
	/* get informations from given vector */
	int size, size_local, low, high;
	TRY( VecGetSize(x.get_vector(), &size) );
	TRY( VecGetLocalSize(x.get_vector(), &size_local) );
	TRY( VecGetOwnershipRange(x.get_vector(), &low, &high) );

	this->T = size/(double)K;
	this->Tlocal = size_local/(double)K;
	this->Tbegin = low/(double)K;
	this->Tend = high/(double)K;

	/* get ranges of all processors - necessary to compute overlaping indexes */
	const int *ranges;
	ranges = (int*)malloc((GlobalManager.get_size()+1)*sizeof(int));
    TRY( VecGetOwnershipRanges(x.get_vector(),&ranges) );
	int myrank = GlobalManager.get_rank();

	int N, n;
	N = x.size(); /* length of whole matrix N = K*T */
	n = x.local_size();

	/* create matrix */
	TRY( MatCreate(PETSC_COMM_WORLD,&A_petsc) );
	TRY( MatSetSizes(A_petsc,n,n,N,N) );

	#ifndef USE_GPU
		TRY( MatSetType(A_petsc,MATMPIAIJ) ); 
	#else
		TRY( MatSetType(A_petsc,MATAIJCUSPARSE) ); 
	#endif

	/* compute preallocation of number of non-zero elements in matrix */
	TRY( MatMPIAIJSetPreallocation(A_petsc,3,NULL,2,NULL) ); 
	TRY( MatSeqAIJSetPreallocation(A_petsc,3,NULL) );

	TRY( MatSetFromOptions(A_petsc) ); 
	
//	#pragma omp parallel for
	for(int k=0; k < K; k++){
		for(int t=0;t<Tlocal;t++){
			int tglobal = this->Tbegin+t;
			int diag_idx = tglobal*K + k;
			int row_idx =  tglobal*K + k;

			int Wsum;

			/* compute sum of W entries in row */
			if(T > 1){
				if(tglobal > 0 && tglobal < T-1){
					Wsum = 2;
				} else {
					Wsum = 1;
				}
			} else {
				Wsum = 1;
			}
			
			/* diagonal entry */
			TRY( MatSetValue(A_petsc, row_idx, diag_idx, Wsum, INSERT_VALUES) );

			/* my nondiagonal entries */
			if(T>1){
				if(tglobal > 0) {
					/* t - 1 */
					TRY( MatSetValue(A_petsc, row_idx, diag_idx-K, -1.0, INSERT_VALUES) );
				}
				if(tglobal < T-1) {
					/* t + 1 */
					TRY( MatSetValue(A_petsc, row_idx, diag_idx+K, -1.0, INSERT_VALUES) );
				}
			}

		}

	}

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
	output << " - T:     " << T << std::endl;
	output << " - K:     " << K << std::endl;
	output << " - size:  " << T*K << std::endl;

	output << " - alpha: " << alpha << std::endl;

	LOG_FUNC_END	
}

template<class VectorBase>
void BlockLaplaceSparseMatrix<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global << this->get_name() << std::endl;
	output_global << " - T:     " << T << std::endl;
	output_global.push();
		output_local << " - Tlocal:  " << Tlocal << " (" << Tbegin << "," << Tend << ")" << std::endl;
		output_local.synchronize();
	output_global.pop();
	
	output_global << " - K:     " << K << std::endl;
	output_global << " - size:  " << T*K << std::endl;
	output_global.push();
		output_local << " - sizelocal:  " << Tlocal*K << std::endl;
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
	return this->K;
}

template<class VectorBase>
int BlockLaplaceSparseMatrix<VectorBase>::get_T() const { 
	return this->T;
}

template<class VectorBase>
int BlockLaplaceSparseMatrix<VectorBase>::get_Tlocal() const { 
	return this->Tlocal;
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
