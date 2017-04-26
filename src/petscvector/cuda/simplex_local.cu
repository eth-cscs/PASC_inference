/** @file simplex_local.cu
 *  @brief simplex feasible set for petscvector with local problems
 *
 *  cuda implementation
 * 
 *  @author Lukas Pospisil
 */

#include "algebra/feasibleset/simplex_local.h"

namespace pascinference {
namespace algebra {


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
			if (x_sorted[t*K+i] < x_sorted[t*K + (i - 1)]){
				swap = x_sorted[t*K + i];
				x_sorted[t*K + i] = x_sorted[t*K + (i - 1)];
				x_sorted[t*K + (i - 1)] = swap;
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
			if(x[t*K+k] < 0.0){
				is_inside = false;
			}
			sum += x[t*K + k];
			
			x_sorted[t*K + k] = x[t*K + k];
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
					sum_y += x_sorted[t*K + j];
				}
				
				ti = (sum_y - 1.0)/(double)(K-i);
				if(ti >= x_sorted[t*K + (i-1)]){
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
				ti = x[t*K + k] - t_hat;	
				if(ti > 0.0){
					x[t*K + k] = ti;
				} else {
					x[t*K + k] = 0.0;
				}
			}
		}
		
	}

	/* if t >= T then relax and do nothing */	
}



}
} /* end namespace */

