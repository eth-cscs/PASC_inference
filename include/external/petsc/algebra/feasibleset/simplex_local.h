/** @file simplex_local.h
 *  @brief simplex feasible set for petscvector with local problems
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_SIMPLEXFEASIBLESET_LOCAL_PETSC_H
#define	PASC_SIMPLEXFEASIBLESET_LOCAL_PETSC_H

#include "pascinference.h"
#include "algebra/feasibleset/simplex_local.h"

typedef petscvector::PetscVector PetscVector;

namespace pascinference {
namespace algebra {

/* constructor */
template<class VectorBase>
SimplexFeasibleSet_Local<VectorBase>::SimplexFeasibleSet_Local(int T, int K){
	LOG_FUNC_BEGIN

	/* set initial content */
	this->T = T;
	this->K = K;

	#ifdef USE_CUDA
		/* allocate space for sorting */
		gpuErrchk( cudaMalloc((void **)&x_sorted,K*T*sizeof(double)) );
		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize,kernel_project, 0, 0) );
		gridSize = (T + blockSize - 1)/ blockSize;
	#endif

	LOG_FUNC_END
}

/* general destructor */
template<class VectorBase>
SimplexFeasibleSet_Local<VectorBase>::~SimplexFeasibleSet_Local(){
	LOG_FUNC_BEGIN
	
	#ifdef USE_CUDA
		/* destroy space for sorting */
		gpuErrchk( cudaFree(&x_sorted) );
	#endif	
	
	LOG_FUNC_END	
}

/* print info about feasible set */
template<class VectorBase>
void SimplexFeasibleSet_Local<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN
	
	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - nmb of subsets:     " << T << std::endl;
	output <<  " - size of subset:     " << K << std::endl;

	#ifdef USE_CUDA
		output <<  " - blockSize:   " << blockSize << std::endl;
		output <<  " - gridSize:    " << gridSize << std::endl;
		output <<  " - minGridSize: " << minGridSize << std::endl;
	#endif

	LOG_FUNC_END
}

template<class VectorBase>
std::string SimplexFeasibleSet_Local<VectorBase>::get_name() const {
	return "SimplexFeasibleSet_Local";
}

template<>
void SimplexFeasibleSet_Local<PetscVector>::project(GeneralVector<PetscVector> &x) {
	LOG_FUNC_BEGIN
	
	/* get local array */
	double *x_arr;
	
	#ifdef USE_CUDA
		TRYCXX( VecCUDAGetArrayReadWrite(x.get_vector(),&x_arr) );

		/* use kernel to compute projection */
		//TODO: here should be actually the comparison of Vec type! not simple use_gpu
		kernel_project<<<gridSize, blockSize>>>(x_arr,x_sorted,T,K);
		gpuErrchk( cudaDeviceSynchronize() );
		MPI_Barrier( MPI_COMM_WORLD );

		TRYCXX( VecCUDARestoreArrayReadWrite(x.get_vector(),&x_arr) );
	#else
		TRYCXX( VecGetArray(x.get_vector(),&x_arr) );
	
		/* use openmp */
		#pragma omp parallel for
		for(int t=0;t<T;t++){
			project_sub(x_arr,t,T,K);
		}

		TRYCXX( VecRestoreArray(x.get_vector(),&x_arr) );
	#endif

	TRYCXX( PetscBarrier(NULL) );

	LOG_FUNC_END
}

template<class VectorBase>
void SimplexFeasibleSet_Local<VectorBase>::sort_bubble(double *x, int n){
	int i;
	int m = n;
	int mnew;
	double swap;

	while(m > 0){
		/* Iterate through x */
		mnew = 0;
		for(i=1;i<m;i++){
			/* Swap elements in wrong order */
			if (x[i] < x[i - 1]){
				swap = x[i];
				x[i] = x[i-1];
				x[i-1] = swap;
				mnew = i;
			}
        }
		m = mnew;
    }
}

template<class VectorBase>
void SimplexFeasibleSet_Local<VectorBase>::project_sub(double *x, int t, int T, int K){
	if(t<T){ /* maybe we call more than T kernels */
		int k;

		bool is_inside = true;
		double sum = 0.0;
	
		/* control inequality constraints */
		for(k = 0; k < K; k++){ // TODO: could be performed parallely  
			if(x[t*K+k] < 0.0){
				is_inside = false;
			}
			sum += x[t*K+k];
		}

		/* control equality constraints */
		if(sum != 1){ 
			is_inside = false;
		}

		/* if given point is not inside the feasible domain, then do projection */
		if(!is_inside){
			int j,i;
			/* compute sorted x_sub */
			double *y = new double[K];
			double sum_y;
			for(k=0;k<K;k++){
				y[k] = x[t*K+k]; 
			}
			sort_bubble(y,K);

			/* now perform analytical solution of projection problem */	
			double t_hat = 0.0;
			i = K - 1;
			double ti;

			while(i >= 1){
				/* compute sum(y) */
				sum_y = 0.0;
				for(j=i;j<K;j++){ /* sum(y(i,n-1)) */
					sum_y += y[j];
				}
				
				ti = (sum_y - 1.0)/(double)(K-i);
				if(ti >= y[i-1]){
					t_hat = ti;
					i = -1; /* break */
				} else {
					i = i - 1;
				}
			}

			if(i == 0){
				t_hat = (sum-1.0)/(double)K; /* uses sum=sum(x_sub) */
			}
    
			for(k = 0; k < K; k++){ // TODO: could be performed parallely  
				/* (*x_sub)(i) = max(*x_sub-t_hat,0); */
				ti = x[t*K+k] - t_hat;	
				if(ti > 0.0){
					x[t*K+k] = ti;
				} else {
					x[t*K+k] = 0.0;
				}
			}
			
			delete y;
		}
		
	}

	/* if t >= T then relax and do nothing */	
}


}
} /* end of namespace */


#endif
