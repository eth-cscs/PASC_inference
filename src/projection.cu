#include "projection.h"


/* -------- HostVector ---------- */

void get_projection(HostVector<Scalar> & x, int K){

	int N = x.size();
	int T = N/K; /* length of vectors */

	int t,k;
	Scalar x_sub[K];  /* GammaVector x_sub(K); */

	#pragma omp parallel for private(t)	
	for(t=0;t<T;t++){
		/* cut x_sub from x */
		for(k=0;k<K;k++){
			x_sub[k] = x(k*T+t);
		}
		
		/* compute subprojection */
		get_projection_sub(x_sub, K);

		/* add x_sub back to x */
		for(k=0;k<K;k++){
			x(k*T+t) = x_sub[k];
		}
	}
}


/* project x_sub to feasible set defined by equality and inequality constraints
 * sum(x_sub) = 1
 * x_sub >= 0
 */
void get_projection_sub(Scalar *x_sub, int n){
	int i;

	bool is_inside = true;
	Scalar sum = 0.0;
	
	/* control inequality constraints */
	for(i = 0; i < n; i++){ // TODO: could be performed parallely  
		if(x_sub[i] < 0.0){
			is_inside = false;
		}
		sum += x_sub[i];
	}

	/* control equality constraints */
	if(sum != 1){ 
		is_inside = false;
	}


	/* if given point is not inside the feasible domain, then do projection */
	if(!is_inside){
		int j;
		/* compute sorted x_sub */
		Scalar y[n], sum_y;
		for(i=0;i<n;i++){
			y[i] = x_sub[i]; 
		}
		sort_bubble(y,n);

		/* now perform analytical solution of projection problem */	
		Scalar t_hat = 0.0;
		i = n - 1;
		Scalar ti;

		while(i >= 1){
			/* compute sum(y) */
			sum_y = 0.0;
			for(j=i;j<n;j++){ /* sum(y(i,n-1)) */
				sum_y += y[j];
			}
				
			ti = (sum_y - 1.0)/(Scalar)(n-i);
			if(ti >= y[i-1]){
				t_hat = ti;
				i = -1; /* break */
			} else {
				i = i - 1;
			}
		}

		if(i == 0){
			t_hat = (sum-1.0)/(Scalar)n; /* uses sum=sum(x_sub) */
		}
    
		for(i = 0; i < n; i++){ // TODO: could be performed parallely  
			/* (*x_sub)(i) = max(*x_sub-t_hat,0); */
			ti = x_sub[i] - t_hat;	
			if(ti > 0.0){
				x_sub[i] = ti;
			} else {
				x_sub[i] = 0.0;
			}
		}
	}
}

/* sort values of scalar vector */
void sort_bubble(Scalar *x, int n){
	int i;
	int m = n;
	int mnew;
	Scalar swap;

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


/* -------- DeviceVector ---------- */
#ifdef USE_GPU

void get_projection(DeviceVector<Scalar> & x, int K){

	int N = x.size();
	int T = N/K; /* length of vectors */

	/* call projection using kernel */
	Scalar *xp = x.pointer();
	
	// TODO: compute optimal nmb of threads/kernels
	kernel_get_projection_sub<<<T, 1>>>(xp,T,K);
	
	/* synchronize kernels, if there is an error with cuda, then it will appear here */ 
	gpuErrchk( cudaDeviceSynchronize() );	
	

}


__device__
void device_sort_bubble(Scalar *x, int n){
	int i;
	int m = n;
	int mnew;
	Scalar swap;

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

__global__
void kernel_get_projection_sub(Scalar *x, int T, const int K){
	/* compute my id */
	int t = blockIdx.x*blockDim.x + threadIdx.x;

	if(t<T){
		int k;

		bool is_inside = true;
		Scalar sum = 0.0;
	
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
			Scalar *y = new Scalar[K];
			Scalar sum_y;
			for(k=0;k<K;k++){
				y[k] = x[k*T+t]; 
			}
			device_sort_bubble(y,K);

			/* now perform analytical solution of projection problem */	
			Scalar t_hat = 0.0;
			i = K - 1;
			Scalar ti;

			while(i >= 1){
				/* compute sum(y) */
				sum_y = 0.0;
				for(j=i;j<K;j++){ /* sum(y(i,n-1)) */
					sum_y += y[j];
				}
				
				ti = (sum_y - 1.0)/(Scalar)(K-i);
				if(ti >= y[i-1]){
					t_hat = ti;
					i = -1; /* break */
				} else {
					i = i - 1;
				}
			}

			if(i == 0){
				t_hat = (sum-1.0)/(Scalar)K; /* uses sum=sum(x_sub) */
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
	

	/* if t >= N then relax and do nothing */	

}

#endif

/* -------- PetscVector ---------- */
#ifdef USE_PETSC

bool petsc_projection_init = false; /* if the initialization of projection was not performed, then = false */
IS petsc_projection_is; 
int petsc_projection_Townership_low, petsc_projection_Townership_high;

void get_projection(PetscVector & x, int K){
	
	int T = x.size()/(double)K; /* length of time-series */

	/* initialization - how much I will compute? */
	if(!petsc_projection_init){
		if(DEBUG_MODE >= 5){
			Message("  - initialization of projection");
		}
		
		petsc_projection_init = true;
		
		/* try to make a global vector of length T and then get the indexes of begin and end of local portion */
		PetscVector TVector(T);

		/* get the ownership range - now I know how much I will calculate from the time-series */
		TVector.get_ownership(&petsc_projection_Townership_low,&petsc_projection_Townership_high);

		/* destroy testing vector - it is useless now */
//		~TVector();

	}

	if(DEBUG_MODE >= 5){
		std::cout << "my ownership: [" << petsc_projection_Townership_low << ", " << petsc_projection_Townership_high << "]" << std::endl;
	}

	int t;
	PetscVector x_sub;
	Scalar *x_sub_arr;
	/* go throught local portion of time-serie and perform the projection */
	for(t = petsc_projection_Townership_low; t < petsc_projection_Townership_high; t++){
		/* prepare index set [low, low + T, ... , low + (K-1)*T ] */
		ISCreateStride(PETSC_COMM_SELF, K, petsc_projection_Townership_low + t, T, &petsc_projection_is);

		/* get the subvector from global vector */
		x_sub = x(petsc_projection_is);

		/* get the array */
		x_sub.get_array(&x_sub_arr);

		/* perform the projection on this subvector */
		get_projection_sub(x_sub_arr, K);

		/* print the array of subvector */
		if(DEBUG_MODE >= 5){
			int i;
			std::cout << "xsub_" << t << " = [ ";
			for(i=0;i<K;i++){
				std::cout << x_sub_arr[i];
				if(i < K-1) std::cout << ", ";
			}
			std::cout << " ]" << std::endl;
		}

		/* restore the array */
		x_sub.restore_array(&x_sub_arr);

		x(petsc_projection_is) = x_sub;

	}

}

#endif
