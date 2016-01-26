#include "projection.h"

/* projection with timer */
void get_projection(GammaVector<Scalar> *x, int K, double *time_to_add){
	timer.start(); /* add to projection time */
	
	get_projection(x, K);
	
	(*time_to_add) += timer.stop();
}

void get_projection(GammaVector<Scalar> *x, int K){

	int t,k;
	int N = (*x).size();
	int T = N/K; /* length of vectors */
	Scalar x_sub[K];  /* GammaVector<Scalar> x_sub(K); */

#ifdef USE_GPU
	
	
#else
	#pragma omp parallel for private(t)	
	for(t=0;t<T;t++){
		/* cut x_sub from x */
		for(k=0;k<K;k++){
			x_sub[k] = (*x)(k*T+t);
		}
		
		/* compute subprojection */
		get_projection_sub(x_sub, K);

		/* add x_sub back to x */
		for(k=0;k<K;k++){
			(*x)(k*T+t) = x_sub[k];
		}
	}
#endif
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
