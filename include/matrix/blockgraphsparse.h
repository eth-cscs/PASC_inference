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
		int T; /**< dimension of each block */
		int R; /**< number of vertices = number of blocks in row,col */
		int K; /**< number of diagonal blocks */

		/* T decomposition */
		int DD_T_size; /**< time - number of domains for decomposition */
		int DD_T_ranges; /**< time - local ownership */

		/* R decomposition */
		int DD_R_size; /**< space - number of domains for decomposition */
		int *DD_R_affiliation; /**< space - domain affiliation of vertices */
		int *DD_R_local_idx; /**< space - local indexes of vertices */
		int *DD_R_lengths; /**< space - array of local lengths */


		int Tlocal; /**< local dimension of each block */
		int Tbegin; /**< ownership begin */
		int Tend; /**< ownership end */
		
		double alpha; /**< general matrix multiplicator */
		
		BGMGraph *graph; /**< graph with stucture of the matrix */
		
		#ifdef USE_PETSCVECTOR
			/* Petsc stuff */ 
			PetscMatrix A_petsc;
		#endif

		GeneralVector<VectorBase> *coeffs; /**< vector of coefficient for each block */
		
	public:
		BlockGraphSparseMatrix(const VectorBase &x, BGMGraph &new_graph, int K, double alpha=1.0, GeneralVector<VectorBase> *new_coeffs=NULL);
		~BlockGraphSparseMatrix(); /* destructor - destroy inner matrix */

		void print(ConsoleOutput &output) const; /* print matrix */
		void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		
		void printcontent(ConsoleOutput &output) const;
		std::string get_name() const;
		
		void matmult(VectorBase &y, const VectorBase &x) const; /* y = A*x */

		int get_R() const;
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
BlockGraphSparseMatrix<VectorBase>::BlockGraphSparseMatrix(const VectorBase &x, BGMGraph &new_graph, int K, double alpha, GeneralVector<VectorBase> *new_coeffs){
	LOG_FUNC_BEGIN

	this->graph = &new_graph;
	this->R = this->graph->get_n();

	this->K = K;
	this->alpha = alpha;
	
	this->coeffs = new_coeffs;
	
	/* get informations from given vector */
	int N, n, low, high;
	TRY( VecGetSize(x.get_vector(), &N) );
	TRY( VecGetLocalSize(x.get_vector(), &n) );
	TRY( VecGetOwnershipRange(x.get_vector(), &low, &high) );

	this->T = N/(double)(R*K);
	this->Tlocal = n/(double)(R*K);
	this->Tbegin = low/(double)(R*K);
	this->Tend = high/(double)(R*K);

	int* neighbor_nmbs = graph->get_neighbor_nmbs();
	int **neightbor_ids = graph->get_neighbor_ids();

	/* get ranges of all processors - necessary to compute overlaping indexes */
	const int *ranges;
	ranges = (int*)malloc((GlobalManager.get_size()+1)*sizeof(int));
	TRY( VecGetOwnershipRanges(x.get_vector(),&ranges) );
	int myrank = GlobalManager.get_rank();

	/* create matrix */
	TRY( MatCreate(PETSC_COMM_WORLD,&A_petsc) );
	TRY( MatSetSizes(A_petsc,n,n,N,N) );

	#ifndef USE_GPU
		TRY( MatSetType(A_petsc,MATMPIAIJ) ); 
	#else
		TRY( MatSetType(A_petsc,MATAIJCUSPARSE) ); 
	#endif

	/* compute preallocation of number of non-zero elements in matrix */
	TRY( MatMPIAIJSetPreallocation(A_petsc,3*(1+2*this->graph->get_m_max()),NULL,2*(this->graph->get_m_max()+1),NULL) ); 
	TRY( MatSeqAIJSetPreallocation(A_petsc,3*(1+2*this->graph->get_m_max()),NULL) );

	TRY( MatSetFromOptions(A_petsc) ); 
	
//	#pragma omp parallel for
	for(int k=0; k < K; k++){
		for(int r=0; r < R; r++){
			for(int t=0;t<Tlocal;t++){
				int tglobal = this->Tbegin+t;
				int diag_idx = tglobal*R*K + r*K + k;
				int row_idx =  tglobal*R*K + r*K + k;

				int Wsum;

				/* compute sum of W entries in row */
				if(tglobal == 0 || tglobal == T-1){
					if(T > 1){
						Wsum = 2*neighbor_nmbs[r]+1; /* +1 for diagonal block */
					} else {
						Wsum = neighbor_nmbs[r];
					}
				} else {
					if(T > 1){
						Wsum = 3*neighbor_nmbs[r]+2; /* +2 for diagonal block */
					} else {
						Wsum = (neighbor_nmbs[r])+1; /* +1 for diagonal block */
					}
				}
				
				/* diagonal entry */
				TRY( MatSetValue(A_petsc, row_idx, diag_idx, Wsum, INSERT_VALUES) );

				/* my nondiagonal entries */
				if(T>1){
					if(tglobal > 0) {
						/* t - 1 */
						TRY( MatSetValue(A_petsc, row_idx, diag_idx-R*K, -1.0, INSERT_VALUES) );
					}
					if(tglobal < T-1) {
						/* t + 1 */
						TRY( MatSetValue(A_petsc, row_idx, diag_idx+R*K, -1.0, INSERT_VALUES) );
					}
				}

				/* non-diagonal neighbor entries */
				int neighbor;
				for(neighbor=0;neighbor<neighbor_nmbs[r];neighbor++){
//					int idx2 = Tbegin*K*R + k*Tlocal*R + *Tlocal + t;
					int idx2 = tglobal*R*K + (neightbor_ids[r][neighbor])*K + k;

					TRY( MatSetValue(A_petsc, row_idx, idx2, -1.0, INSERT_VALUES) );
					if(tglobal > 0) {
						TRY( MatSetValue(A_petsc, row_idx, idx2-R*K, -1.0, INSERT_VALUES) );
					}
					if(tglobal < T-1) {
						TRY( MatSetValue(A_petsc, row_idx, idx2+R*K, -1.0, INSERT_VALUES) );
					}
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
	output << " - T:     " << T << std::endl;
	output << " - R:     " << R << std::endl;
	output << " - K:     " << K << std::endl;
	output << " - size:  " << T*R*K << std::endl;

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
	output_global << " - T:     " << T << std::endl;
	output_global.push();
		output_local << " - Tlocal:  " << Tlocal << " (" << Tbegin << "," << Tend << ")" << std::endl;
		output_local.synchronize();
	output_global.pop();
	
	output_global << " - R:     " << R << std::endl;
	output_global << " - K:     " << K << std::endl;
	output_global << " - size:  " << T*R*K << std::endl;
	output_global.push();
		output_local << " - sizelocal:  " << Tlocal*R*K << std::endl;
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
		double *coeffs_arr;
		TRY( VecGetArray(coeffs->get_vector(),&coeffs_arr) );

		Vec xk_Vec;
		Vec x_Vec = y.get_vector();
		IS xk_is;
		double coeff;
		
		/* get vector corresponding to coeff */
		for(int k=0;k<K;k++){
			TRY( ISCreateStride(PETSC_COMM_WORLD, R*Tlocal, Tbegin*K*R + k, K, &xk_is) );
			TRY( VecGetSubVector(x_Vec, xk_is, &xk_Vec) );
			
			coeff = alpha*coeffs_arr[k]*coeffs_arr[k];
			TRY( VecScale(xk_Vec, coeff) );
			
			TRY( VecRestoreSubVector(x_Vec, xk_is, &xk_Vec) );
			TRY( ISDestroy(&xk_is) );
		}

		TRY( VecRestoreArray(coeffs->get_vector(),&coeffs_arr) );
	} else {
	    TRY( VecScale(y.get_vector(), this->alpha) );
	}

	LOG_FUNC_END
}

template<class VectorBase>
int BlockGraphSparseMatrix<VectorBase>::get_R() const { 
	return this->R;
}

template<class VectorBase>
int BlockGraphSparseMatrix<VectorBase>::get_K() const { 
	return this->K;
}

template<class VectorBase>
int BlockGraphSparseMatrix<VectorBase>::get_T() const { 
	return this->T;
}

template<class VectorBase>
int BlockGraphSparseMatrix<VectorBase>::get_Tlocal() const { 
	return this->Tlocal;
}

template<class VectorBase>
double BlockGraphSparseMatrix<VectorBase>::get_alpha() const { 
	return this->alpha;
}


} /* end of namespace */


#endif
