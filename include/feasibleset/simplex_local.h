/** @file simplex_local.h
 *  @brief simplex feasible set for petscvector with local problems
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_SIMPLEXFEASIBLESET_LOCAL_H
#define	PASC_SIMPLEXFEASIBLESET_LOCAL_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include "pascinference.h"

#ifndef USE_PETSCVECTOR
 #error 'SIMPLEXFEASIBLESET_LOCAL is for PETSCVECTOR'
#endif

typedef petscvector::PetscVector PetscVector;

#ifdef USE_CUDA
    #include <../src/vec/vec/impls/seq/seqcuda/cudavecimpl.h>
#endif

namespace pascinference {


/** \class SimplexFeasibleSet_Local
 *  \brief class for manipulation with simplex set
 *
 *  Provides structure for manipulation with simplex feasible set, i.e.
 * \f[
 * \Omega = 
 *  \left\lbrace x \in R^{KT}: 
 *    \forall t = 0,\dots,T-1: \sum\limits_{k=0}^{K-1} x_{t+kT} = 1,  
 *    x \geq 0
 *  \right\rbrace
 *	\f]
 * 
*/
class SimplexFeasibleSet_Local: public GeneralFeasibleSet<PetscVector> {
	private:
		
		/** @brief sort array using bubble sort
		 * 
		 * @param x array with values
		 * @param n size of array
		*/ 		
		void sort_bubble(double *x, int n);

		/** @brief projection onto simplex
		 *
		 * computes a projection of a point onto simplex in nD
		 *
		 * take K-dimensional vector x[t,t+T,t+2T,...t+(K-1)T] =: p
		 * and compute projection
		 * P(p) = arg min || p - y ||_2
		 * subject to constraints (which define simplex)
		 * y_0 + ... y_{K-1} = 1
		 * y_i >= 0 for all i=0,...,K-1
		 *
		 * in practical applications K is much more lower number than T
		 * K - number of clusters (2 - 10^2)
		 * T - length of time-series (10^5 - 10^9) 
		 * 
		 * @param t where my subvector starts
		 * @param x values of whole vector in array
		 * @param K size of subvector
		*/ 		
		void project_sub(double *x, int t, int T, int K);

		#ifdef USE_GPU
			double *x_sorted; /**< for manipulation with sorted data on GPU */
			int blockSize; /**< block size returned by the launch configurator */
			int minGridSize; /**< the minimum grid size needed to achieve the maximum occupancy for a full device launch */
			int gridSize; /**< the actual grid size needed, based on input size */
		#endif

		int T; /**< number of disjoint simplex subsets */
		int K_local; /**< size of each simplex subset */
				
	public:
		/** @brief default constructor
		*/ 	
		SimplexFeasibleSet_Local(int T, int K_local);
		
		/** @brief default destructor
		 */ 
		~SimplexFeasibleSet_Local();

		/** @brief print properties of this feasible set
		 * 
		 * @param output where to print
		 */ 
		void print(ConsoleOutput &output) const;

		/** @brief get name of this feasible set
		 */
		virtual std::string get_name() const;

		/** @brief compute projection onto feasible set
		 * 
		 * @param x point which will be projected
		 */		
		void project(GeneralVector<PetscVector> &x);
		
};

/* cuda kernels cannot be a member of class */
#ifdef USE_CUDA
__device__ void device_sort_bubble(double *x_sorted, int t, int T, int K);

/** @brief kernel projection onto simplex
 *
 * computes a projection of a point onto simplex in nD
 *
 * take K-dimensional vector x[t,t+T,t+2T,...t+(K-1)T] =: p
 * and compute projection
 * P(p) = arg min || p - y ||_2
 * subject to constraints (which define simplex)
 * y_0 + ... y_{K-1} = 1
 * y_i >= 0 for all i=0,...,K-1
 *
 * in practical applications K is much more lower number than T
 * K - number of clusters (2 - 10^2)
 * T - length of time-series (10^5 - 10^9) 
 * 
 * @param x values of whole vector in array
 * @param x_sorted allocated vector for manipulating with sorted x
 * @param T parameter of vector length
 * @param K parameter of vector length
 */ 
__global__ void kernel_project(double *x, double *x_sorted, int T, int K);
#endif

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
SimplexFeasibleSet_Local::SimplexFeasibleSet_Local(int Tnew, int Knew){
	LOG_FUNC_BEGIN

	/* set initial content */
	this->T = Tnew;
	this->K_local = Knew;

	#ifdef USE_GPU
		/* allocate space for sorting */
		gpuErrchk( cudaMalloc((void **)&x_sorted,K_local*T*sizeof(double)) );
		gpuErrchk( cudaOccupancyMaxPotentialBlockSize( &minGridSize, &blockSize,kernel_project, 0, 0) );
		gridSize = (T + blockSize - 1)/ blockSize;
	#endif

	LOG_FUNC_END
}

/* general destructor */
SimplexFeasibleSet_Local::~SimplexFeasibleSet_Local(){
	LOG_FUNC_BEGIN
	
	#ifdef USE_GPU
		/* destroy space for sorting */
		gpuErrchk( cudaFree(&x_sorted) );
	#endif	
	
	LOG_FUNC_END	
}

/* print info about feasible set */
void SimplexFeasibleSet_Local::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN
	
	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - T:           " << T << std::endl;
	output <<  " - K_local:     " << K_local << std::endl;

	#ifdef USE_GPU
		output <<  " - blockSize:   " << blockSize << std::endl;
		output <<  " - gridSize:    " << gridSize << std::endl;
		output <<  " - minGridSize: " << minGridSize << std::endl;
	#endif

	LOG_FUNC_END
}

std::string SimplexFeasibleSet_Local::get_name() const {
	return "SimplexFeasibleSet_Local";
}

void SimplexFeasibleSet_Local::project(GeneralVector<PetscVector> &x) {
	LOG_FUNC_BEGIN
	
	/* get local array */
	double *x_arr;
	
	#ifdef USE_GPU
		TRY( VecCUDAGetArrayReadWrite(x.get_vector(),&x_arr) );

		/* use kernel to compute projection */
		//TODO: here should be actually the comparison of Vec type! not simple use_gpu
		kernel_project<<<gridSize, blockSize>>>(x_arr,x_sorted,T,K_local);
		gpuErrchk( cudaDeviceSynchronize() );

		TRY( VecCUDARestoreArrayReadWrite(x.get_vector(),&x_arr) );
	#else
		TRY( VecGetArray(x.get_vector(),&x_arr) );
	
		/* use openmp */
		#pragma omp parallel for
		for(int t=0;t<T;t++){
			project_sub(x_arr,t,T,K_local);
		}

		TRY( VecRestoreArray(x.get_vector(),&x_arr) );
	#endif

	LOG_FUNC_END
}


void SimplexFeasibleSet_Local::sort_bubble(double *x, int n){
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

void SimplexFeasibleSet_Local::project_sub(double *x, int t, int T, int K){
	if(t<T){ /* maybe we call more than T kernels */
		int k;

		bool is_inside = true;
		double sum = 0.0;
	
		/* control inequality constraints */
		for(k = 0; k < K; k++){ // TODO: could be performed parallely  
			if(x[k*T+t] < 0.0){
				is_inside = false;
			}
			sum += x[k*T+t];
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
				y[k] = x[k*T+t]; 
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
				ti = x[k*T+t] - t_hat;	
				if(ti > 0.0){
					x[k*T+t] = ti;
				} else {
					x[k*T+t] = 0.0;
				}
			}
			
			delete y;
		}
		
	}

	/* if t >= T then relax and do nothing */	
}




/* kernels in cuda */
#ifdef USE_GPU

__device__ void device_sort_bubble(double *x_sorted, int t, int T, int K){
	int i;
	int m=K;
	int mnew;
	double swap;

	while(m > 0){
		/* Iterate through x */
		mnew = 0;
		for(i=1;i<m;i++){
			/* Swap elements in wrong order */
			if (x_sorted[i*T+t] < x_sorted[(i - 1)*T+t]){
				swap = x_sorted[i*T+t];
				x_sorted[i*T+t] = x_sorted[(i-1)*T+t];
				x_sorted[(i-1)*T+t] = swap;
				mnew = i;
			}
	        }
		m = mnew;
	}
}

__global__ void kernel_project(double *x, double *x_sorted, int T, int K){
	int t = blockIdx.x*blockDim.x + threadIdx.x; /* thread t */

	int k;
	if(t<T){
		bool is_inside = true;
		double sum = 0.0;
	
		/* control inequality constraints */
		for(k = 0; k < K; k++){ // TODO: could be performed parallely  
			if(x[k*T+t] < 0.0){
				is_inside = false;
			}
			sum += x[k*T+t];
			
			x_sorted[k*T+t] = x[k*T+t];
		}

		/* control equality constraints */
		if(sum != 1){ 
			is_inside = false;
		}

		/* if given point is not inside the feasible domain, then do projection */
		if(!is_inside){
			int j,i;
			/* compute sorted x_sub */
			double sum_y;
			device_sort_bubble(x_sorted,t,T,K);

			/* now perform analytical solution of projection problem */	
			double t_hat = 0.0;
			i = K - 1;
			double ti;

			while(i >= 1){
				/* compute sum(y) */
				sum_y = 0.0;
				for(j=i;j<K;j++){ /* sum(y(i,n-1)) */
					sum_y += x_sorted[j*T+t];
				}
				
				ti = (sum_y - 1.0)/(double)(K-i);
				if(ti >= x_sorted[(i-1)*T+t]){
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
				ti = x[k*T+t] - t_hat;	
				if(ti > 0.0){
					x[k*T+t] = ti;
				} else {
					x[k*T+t] = 0.0;
				}
			}
		}
		
	}

	/* if t >= T then relax and do nothing */	
}

#endif



} /* end namespace */

#endif
