#ifndef BLOCKLAPLACEFREEMATRIX_H
#define	BLOCKLAPLACEFREEMATRIX_H

#include "pascinference.h"

#ifndef USE_PETSCVECTOR
 #error 'BLOCKLAPLACEFREEMATRIX is for PETSCVECTOR only, sorry'
#endif

#ifdef USE_CUDA
    #include <../src/vec/vec/impls/seq/seqcuda/cudavecimpl.h>
#endif

typedef petscvector::PetscVector PetscVector;

namespace pascinference {
namespace algebra {

template<class VectorBase>
class BlockLaplaceFreeMatrix: public GeneralMatrix<VectorBase> {
	private:
		Decomposition *decomposition;

		double alpha; /**< general matrix multiplicator */
		
		#ifdef USE_CUDA
			int blockSize; /**< block size returned by the launch configurator */
			int minGridSize; /**< the minimum grid size needed to achieve the maximum occupancy for a full device launch */
			int gridSize; /**< the actual grid size needed, based on input size */
		#endif	

		/* stuff for 3diag mult */
		#ifdef USE_PETSCVECTOR
			int left_overlap, right_overlap;
			IS left_overlap_is;
			IS right_overlap_is;
		#endif
		
	public:
		BlockLaplaceFreeMatrix(Decomposition &decomposition, double alpha=1.0);
		~BlockLaplaceFreeMatrix(); /* destructor - destroy inner matrix */

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

#ifdef USE_CUDA
__global__ void kernel_mult(double* y_arr, double* x_arr, double *left_overlap_arr, double *right_overlap_arr, int T, int Tlocal, int Tbegin, int K, double alpha);
#endif


/* --------------- IMPLEMENTATION -------------------- */
// TODO: move to implementation
template<class VectorBase>
std::string BlockLaplaceFreeMatrix<VectorBase>::get_name() const {
	return "BlockLaplaceFreeMatrix";
}

template<>
BlockLaplaceFreeMatrix<PetscVector>::BlockLaplaceFreeMatrix(Decomposition &new_decomposition, double alpha){
	LOG_FUNC_BEGIN

	this->decomposition = &new_decomposition;
	this->alpha = alpha;
	
	#ifdef USE_CUDA
		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize,kernel_mult, 0, 0) );
		gridSize = (get_K()*get_Tlocal() + blockSize - 1)/ blockSize;
	#endif

	/* create IS for 3diag mult */
	if(decomposition->get_Tbegin() > 0){
		left_overlap = 1;
		TRY( ISCreateStride(PETSC_COMM_SELF, get_K(), (decomposition->get_Tbegin()-1)*get_K() , 1, &left_overlap_is) );
	} else {
		left_overlap = 0;
		TRY( ISCreateStride(PETSC_COMM_SELF, 0, 0, 0, &left_overlap_is) );
	}
	if(decomposition->get_Tend() < get_T()){
		right_overlap = 1;
		TRY( ISCreateStride(PETSC_COMM_SELF, get_K(), decomposition->get_Tend()*get_K(), 1, &right_overlap_is) );
	} else {
		right_overlap = 0;
		TRY( ISCreateStride(PETSC_COMM_SELF, 0, 0, 0, &right_overlap_is) );
	}

	LOG_FUNC_END
}	

template<class VectorBase>
BlockLaplaceFreeMatrix<VectorBase>::~BlockLaplaceFreeMatrix(){
	LOG_FUNC_BEGIN
	
	/* destroy aux vector */
//	TRY( VecDestroy(&x_aux) );
	
	/* destroy IS for 3diag mult */
//	TRY( ISDestroy(&left_overlap_is) );
//	TRY( ISDestroy(&right_overlap_is) );
	
	LOG_FUNC_END	
}

/* print matrix */
template<class VectorBase>
void BlockLaplaceFreeMatrix<VectorBase>::print(ConsoleOutput &output) const		
{
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	output << " - T:     " << get_T() << std::endl;
	output << " - K:     " << get_K() << std::endl;
	output << " - size:  " << get_T()*get_K() << std::endl;

	output << " - alpha: " << alpha << std::endl;

	#ifdef USE_CUDA
		output <<  " - blockSize:   " << blockSize << std::endl;
		output <<  " - gridSize:    " << gridSize << std::endl;
		output <<  " - minGridSize: " << minGridSize << std::endl;
	#endif
	
	LOG_FUNC_END	
}

template<class VectorBase>
void BlockLaplaceFreeMatrix<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
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

	#ifdef USE_CUDA
		output_global <<  " - blockSize:   " << blockSize << std::endl;
		output_global <<  " - gridSize:    " << gridSize << std::endl;
		output_global <<  " - minGridSize: " << minGridSize << std::endl;
	#endif

	LOG_FUNC_END	
}

/* print matrix */
template<class VectorBase>
void BlockLaplaceFreeMatrix<VectorBase>::printcontent(ConsoleOutput &output) const		
{
	LOG_FUNC_BEGIN
	
	
	LOG_FUNC_END
}

template<class VectorBase>
int BlockLaplaceFreeMatrix<VectorBase>::get_K() const { 
	return decomposition->get_K();
}

template<class VectorBase>
int BlockLaplaceFreeMatrix<VectorBase>::get_T() const { 
	return decomposition->get_T();
}

template<class VectorBase>
int BlockLaplaceFreeMatrix<VectorBase>::get_Tlocal() const { 
	return decomposition->get_Tlocal();
}

template<class VectorBase>
double BlockLaplaceFreeMatrix<VectorBase>::get_alpha() const { 
	return this->alpha;
}

#ifndef USE_CUDA
/* A*x using openmp */
template<>
void BlockLaplaceFreeMatrix<PetscVector>::matmult(PetscVector &y, const PetscVector &x) const { 
	double *y_arr;
	const double *x_arr;
	
	Vec left_overlap_vec;
	Vec right_overlap_vec;
	const double *left_overlap_arr;
	const double *right_overlap_arr;

	/* get subvector and array */
	TRY( VecGetSubVector(x.get_vector(),left_overlap_is,&left_overlap_vec) );
	TRY( VecGetArrayRead(left_overlap_vec,&left_overlap_arr) );

	TRY( VecGetSubVector(x.get_vector(),right_overlap_is,&right_overlap_vec) );
	TRY( VecGetArrayRead(right_overlap_vec,&right_overlap_arr) );
	
	TRY( VecGetArrayRead(x.get_vector(),&x_arr) );
	TRY( VecGetArray(y.get_vector(),&y_arr) );

	int T = get_T();
	int Tlocal = get_Tlocal();
	int Tbegin = decomposition->get_Tbegin();
	int K = get_K();

	/* use openmp */
	#pragma omp parallel for
	for(int y_arr_idx=0;y_arr_idx<Tlocal*K;y_arr_idx++){
		int tlocal = floor(y_arr_idx/(double)(K));
		int k = y_arr_idx - tlocal*K;
		int tglobal = Tbegin+tlocal;
		int x_arr_idx = tlocal*K + k;
		
		int overlap_id = k;

		double value;
		double right_value;
		double left_value;

		/* first row */
		if(tglobal == 0){
			if(tlocal+1 >= Tlocal){
				right_value = right_overlap_arr[overlap_id];
			} else {
				right_value = x_arr[x_arr_idx+K];
			}
			value = alpha*x_arr[x_arr_idx] - alpha*right_value;
		}
		/* common row */
		if(tglobal > 0 && tglobal < T-1){
			if(tlocal+1 >= Tlocal){
				right_value = right_overlap_arr[overlap_id];
			} else {
				right_value = x_arr[x_arr_idx+K];
			}
			if(tlocal-1 < 0){
				left_value = left_overlap_arr[overlap_id];
			} else {
				left_value = x_arr[x_arr_idx-K];
			}
			value = -alpha*left_value + 2*alpha*x_arr[x_arr_idx] - alpha*right_value;
		}
		/* last row */
		if(tglobal == T-1){
			if(tlocal-1 < 0){
				left_value = left_overlap_arr[overlap_id];
			} else {
				left_value = x_arr[x_arr_idx-K];
			}
			value = -alpha*left_value + alpha*x_arr[x_arr_idx];
		}

		y_arr[y_arr_idx] = value;
//		y_arr[x_arr_idx] = k;

	}

	/* restore subvector and array */
	TRY( VecRestoreArrayRead(x.get_vector(),&x_arr) );
	TRY( VecRestoreArray(y.get_vector(),&y_arr) );

	TRY( VecRestoreArrayRead(left_overlap_vec,&left_overlap_arr) );
	TRY( VecRestoreSubVector(x.get_vector(),left_overlap_is,&left_overlap_vec) );

	TRY( VecRestoreArrayRead(right_overlap_vec,&right_overlap_arr) );
	TRY( VecRestoreSubVector(x.get_vector(),right_overlap_is,&right_overlap_vec) );

}

#else

/* A*x using CUDA kernel */
__global__ void kernel_mult(double* y_arr, double* x_arr, double *left_overlap_arr, double *right_overlap_arr, int T, int Tlocal, int Tbegin, int K, double alpha)
{
	/* compute my id */
	int y_arr_idx = blockIdx.x*blockDim.x + threadIdx.x;

	if(y_arr_idx < K*Tlocal){
		int tlocal = floor(y_arr_idx/(double)(K));
		int k = y_arr_idx - tlocal*K;
		int tglobal = Tbegin+tlocal;
		int x_arr_idx = tlocal*K + k;		
		int overlap_id = k;

		double value;
		double right_value;
		double left_value;

		/* first row */
		if(tglobal == 0){
			if(tlocal+1 >= Tlocal){
				right_value = right_overlap_arr[overlap_id];
			} else {
				right_value = x_arr[x_arr_idx+K];
			}
			value = alpha*x_arr[x_arr_idx] - alpha*right_value;
		}
		/* common row */
		if(tglobal > 0 && tglobal < T-1){
			if(tlocal+1 >= Tlocal){
				right_value = right_overlap_arr[overlap_id];
			} else {
				right_value = x_arr[x_arr_idx+K];
			}
			if(tlocal-1 < 0){
				left_value = left_overlap_arr[overlap_id];
			} else {
				left_value = x_arr[x_arr_idx-K];
			}
			value = -alpha*left_value + 2*alpha*x_arr[x_arr_idx] - alpha*right_value;
		}
		/* last row */
		if(tglobal == T-1){
			if(tlocal-1 < 0){
				left_value = left_overlap_arr[overlap_id];
			} else {
				left_value = x_arr[x_arr_idx-K];
			}
			value = -alpha*left_value + alpha*x_arr[x_arr_idx];
		}
		
		y_arr[y_arr_idx] = value; 
	}

	/* if id_row >= K*T then relax and do nothing */	

}

template<>
void BlockLaplaceFreeMatrix<PetscVector>::matmult(PetscVector &y, const PetscVector &x) const {
	double *y_arr;
	double *x_arr;
	
	Vec left_overlap_vec;
	Vec right_overlap_vec;
	double *left_overlap_arr;
	double *right_overlap_arr;

	int T = get_T();
	int Tlocal = get_Tlocal();
	int Tbegin = get_Tbegin();
	int K = get_K();

	/* get subvector and array */
	TRY( VecGetSubVector(x.get_vector(),left_overlap_is,&left_overlap_vec) );
	TRY( VecCUDAGetArrayReadWrite(left_overlap_vec,&left_overlap_arr) );

	TRY( VecGetSubVector(x.get_vector(),right_overlap_is,&right_overlap_vec) );
	TRY( VecCUDAGetArrayReadWrite(right_overlap_vec,&right_overlap_arr) );
	
	TRY( VecCUDAGetArrayReadWrite(x.get_vector(),&x_arr) );
	TRY( VecCUDAGetArrayReadWrite(y.get_vector(),&y_arr) );

	/* call kernel */
	kernel_mult<<<gridSize, blockSize>>>(y_arr,x_arr,left_overlap_arr,right_overlap_arr,T,Tlocal,Tbegin,K,alpha);
	gpuErrchk( cudaDeviceSynchronize() );
	
	/* restore subvector and array */
	TRY( VecCUDARestoreArrayReadWrite(x.get_vector(),&x_arr) );
	TRY( VecCUDARestoreArrayReadWrite(y.get_vector(),&y_arr) );

	TRY( VecCUDARestoreArrayReadWrite(left_overlap_vec,&left_overlap_arr) );
	TRY( VecRestoreSubVector(x.get_vector(),left_overlap_is,&left_overlap_vec) );

	TRY( VecCUDARestoreArrayReadWrite(right_overlap_vec,&right_overlap_arr) );
	TRY( VecRestoreSubVector(x.get_vector(),right_overlap_is,&right_overlap_vec) );
}

#endif


}
} /* end of namespace */


#endif
