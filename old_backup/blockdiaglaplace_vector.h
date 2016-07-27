#ifndef BLOCKDIAGLAPLACEVECTORMATRIX_H
#define	BLOCKDIAGLAPLACEVECTORMATRIX_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include "pascinference.h"

#ifdef USE_PETSCVECTOR
	typedef petscvector::PetscVector PetscVector;
#endif

#ifdef USE_MINLIN
	typedef minlin::threx::HostVector<double> MinlinHostVector;
	typedef minlin::threx::DeviceVector<double> MinlinDeviceVector;
#endif

#ifdef USE_CUDA
    #include <../src/vec/vec/impls/seq/seqcuda/cudavecimpl.h>
#endif

namespace pascinference {

/* laplace matrix */ 
template<class VectorBase>
class BlockDiagLaplaceVectorMatrix: public GeneralMatrix<VectorBase> {
	private:
		int K; /**< number of block */
		int N; /**< size of the matrix */
		int T; /**< size of each block N = K*T */
		double alpha; /**< scale of whole matrix alpha*A */

		#ifdef USE_GPU
			int blockSize; /**< block size returned by the launch configurator */
			int minGridSize; /**< the minimum grid size needed to achieve the maximum occupancy for a full device launch */
			int gridSize; /**< the actual grid size needed, based on input size */
		#endif	
		
	public:
		BlockDiagLaplaceVectorMatrix(const VectorBase &x, int K, int T, double alpha = 1.0); /* constructor from vector and number of blocks */

		~BlockDiagLaplaceVectorMatrix(); /* destructor - destroy inner matrix */

		void print(ConsoleOutput &output) const; /* print matrix */
		std::string get_name() const;

		void matmult(VectorBase &y, const VectorBase &x) const; /* y = A*x */

};

#ifdef USE_GPU
__global__ void kernel_mult(double* Axp, double* xp, int T, int K, double alpha);
#endif



template<class VectorBase>
std::string BlockDiagLaplaceVectorMatrix<VectorBase>::get_name() const {
	return "BlockDiagLaplaceVectorMatrix";
}

template<class VectorBase>
BlockDiagLaplaceVectorMatrix<VectorBase>::BlockDiagLaplaceVectorMatrix(const VectorBase &x, int newK, int newT, double newalpha){
	if(DEBUG_MODE >= 100) coutMaster << "(BlockDiagLaplaceVectorMatrix)CONSTRUCTOR: from PetscVector" << std::endl;

	/* get informations from given vector */
	alpha = newalpha;
	K = newK; /* number of blocks */
	T = newT; /* size of each block */
	N = K*T; /* length of whole matrix N = K*T */
	
	#ifdef USE_GPU
		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize,kernel_mult, 0, 0) );
		gridSize = (K*T + blockSize - 1)/ blockSize;
	#endif
}

template<class VectorBase>
BlockDiagLaplaceVectorMatrix<VectorBase>::~BlockDiagLaplaceVectorMatrix(){
	if(DEBUG_MODE >= 100) coutMaster << "(BlockDiagLaplaceVectorMatrix)DESTRUCTOR" << std::endl;

}

template<class VectorBase>
void BlockDiagLaplaceVectorMatrix<VectorBase>::print(ConsoleOutput &output) const		
{
	if(DEBUG_MODE >= 100) coutMaster << "(BlockDiagLaplaceVectorMatrix)OPERATOR: << print" << std::endl;

	output << "Laplace matrix (";
	output << "K = " << K << ", ";
	output << "T = " << T << ", ";
	output << "N = " << N << ", ";
	output << "alpha = " << alpha << ")" << std::endl;

	#ifdef USE_GPU
		output <<  " - blockSize:   " << blockSize << std::endl;
		output <<  " - gridSize:    " << gridSize << std::endl;
		output <<  " - minGridSize: " << minGridSize << std::endl;
	#endif

}

template<class VectorBase>
void BlockDiagLaplaceVectorMatrix<VectorBase>::matmult(VectorBase &y, const VectorBase &x) const { 
	if(DEBUG_MODE >= 100) coutMaster << "(BlockDiagLaplaceVectorMatrix)FUNCTION: matmult" << std::endl;

	// TODO: maybe y is not initialized, who knows

	y(1,N-2) = alpha*2*x(1,N-2) - alpha*x(0,N-3) - alpha*x(2,N-1);
	
	/* begin and end of each block */
	for(int k=0;k<K;k++){
		y(k*T) = alpha*x(k*T) - alpha*x(k*T+1);
		y((k+1)*T-1) = alpha*x((k+1)*T-1) - alpha*x((k+1)*T-2);
	}
	
}

#ifdef USE_PETSCVECTOR

typedef petscvector::PetscVector PetscVector;

#ifndef USE_GPU
template<>
void BlockDiagLaplaceVectorMatrix<PetscVector>::matmult(PetscVector &y, const PetscVector &x) const { 
	if(DEBUG_MODE >= 100) coutMaster << "(BlockDiagLaplaceVectorMatrix)FUNCTION: matmult" << std::endl;

	// TODO: maybe y is not initialized, who knows
	double *y_arr;
	const double *x_arr;
	TRY( VecGetArray(y.get_vector(),&y_arr) );
	TRY( VecGetArrayRead(x.get_vector(),&x_arr) );

//	int k,t,id_row;
	
	/* use openmp */
	#pragma omp parallel for
	for(int id_row=0;id_row<T*K;id_row++){
		int k = (int)(id_row/(double)T);
		int t = id_row - k*T;

		double value;

		/* first row */
		if(t == 0){
			value = alpha*x_arr[id_row] - alpha*x_arr[id_row+1];
		}
		/* common row */
		if(t > 0 && t < T-1){
			value = -alpha*x_arr[id_row-1] + 2.0*alpha*x_arr[id_row] - alpha*x_arr[id_row+1];
		}
		/* last row */
		if(t == T-1){
			value = -alpha*x_arr[id_row-1] + alpha*x_arr[id_row];
		}
		
		y_arr[id_row] = value;
	}

	TRY( VecRestoreArray(y.get_vector(),&y_arr) );
	TRY( VecRestoreArrayRead(x.get_vector(),&x_arr) );

}
#else

/* A*x using CUDA kernel */
__global__ void kernel_mult(double* y, double* x, int T, int K, double alpha)
{
	/* compute my id */
	int id_row = blockIdx.x*blockDim.x + threadIdx.x;

	if(id_row < K*T){
		/* compute id of cluster */
		int k = (int)(id_row/T);
	
		/* compute id_row in local block */
		int t = id_row-k*T;
		double value;

		/* first row */
		if(t == 0){
			value = x[id_row] - x[id_row+1];
		}
		/* common row */
		if(t > 0 && t < T-1){
			value = -x[id_row-1] + 2.0*x[id_row] - x[id_row+1];
		}
		/* last row */
		if(t == T-1){
			value = -x[id_row-1] + x[id_row];
		}
		
		y[id_row] = alpha*value;
	}

	/* if id_row >= K*T then relax and do nothing */	

}

template<>
void BlockDiagLaplaceVectorMatrix<PetscVector>::matmult(PetscVector &y, const PetscVector &x) const { 
	if(DEBUG_MODE >= 100) coutMaster << "(BlockDiagLaplaceVectorMatrix)FUNCTION: matmult" << std::endl;

	double *y_arr;
	double *x_arr;
	TRY( VecCUDAGetArrayReadWrite(y.get_vector(),&y_arr) );
	TRY( VecCUDAGetArrayReadWrite(x.get_vector(),&x_arr) );

	/* call kernel */
	kernel_mult<<<gridSize, blockSize>>>(y_arr,x_arr,T,K,alpha);
	gpuErrchk( cudaDeviceSynchronize() );

	TRY( VecCUDARestoreArrayReadWrite(y.get_vector(),&y_arr) );
	TRY( VecCUDARestoreArrayReadWrite(x.get_vector(),&x_arr) );

}

#endif


#endif



} /* end of namespace */


#endif