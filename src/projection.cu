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
	GammaVector<Scalar> x_sub(K);

	for(t=0;t<T;t++){	// TODO: this is the place, where the parallel impementation should make a point
		/* cut x_sub from x */
		for(k=0;k<K;k++){
			x_sub(k) = (*x)(k*T+t);
		}
		
		/* compute subprojection */
		get_projection_sub(&x_sub);

		/* add x_sub back to x */
		for(k=0;k<K;k++){
			(*x)(k*T+t) = x_sub(k);
		}
	}
}

/* project x_sub to feasible set defined by equality and inequality constraints
 * sum(x_sub) = 1
 * x_sub >= 0
 */
void get_projection_sub(GammaVector<Scalar> *x_sub){
	int n = x_sub->size(); /* nmb of components of x_sub */
	int i;

	bool is_inside = true;
	
	/* control equality constraints */
	if(sum(*x_sub) != 1){ 
		is_inside = false;
	}
	
	/* control inequality constraints */
	for(i = 0; i < n; i++){ // TODO: could be performed parallely  
		if((*x_sub)(i) < 0.0){
			is_inside = false;
		}
	}

	/* if given point is not inside the feasible domain, then do projection */
	if(!is_inside){
		/* compute sorted x_sub */
		GammaVector<Scalar> y(n);
		for(i=0;i<n;i++){ // TODO: it is really necessary?
			y(i) = (*x_sub)(i); 
		}
		sort_bubble(&y);

		/* now perform analytical solution of projection problem */	
		Scalar t_hat = 0.0;
		i = n - 1;
		Scalar ti;

		while(i >= 1){
			ti = (sum(y(i,n-1)) - 1.0)/(Scalar)(n-i);
			if(ti >= y(i-1)){
				t_hat = ti;
				i = -1; /* break */
			} else {
				i = i - 1;
			}
		}

		if(i == 0){
			t_hat = (sum(y)-1.0)/(Scalar)n;
		}
    
		for(i = 0; i < n; i++){ // TODO: could be performed parallely  
			/* (*x_sub)(i) = max(*x_sub-t_hat,0); */
			ti = (*x_sub)(i) - t_hat;	
			if(ti > 0.0){
				(*x_sub)(i) = ti;
			} else {
				(*x_sub)(i) = 0.0;
			}
		}
	}
}

/* sort values of scalar vector */
void sort_bubble(GammaVector<Scalar> *x){
	int n = x->size();
	int i;
	int nnew;
	Scalar swap;

	while(n > 0){
		/* Iterate through x */
		nnew = 0;
		for(i=1;i<n;i++){
			/* Swap elements in wrong order */
			if ((*x)(i) < (*x)(i - 1)){
				swap = (*x)(i);
				(*x)(i) = (*x)(i-1);
				(*x)(i-1) = swap;
				nnew = i;
			}
        }
		n = nnew;
    }
}
