#include "operations.h"

/* main functions */
void get_dot(Scalar *xx, GammaVector<Scalar> x, GammaVector<Scalar> y){
	*xx = dot(x,y);
}

void get_Ax_laplace(GammaVector<Scalar> *Ax, GammaVector<Scalar> x){
	int N = x.size();
	int t;

	for(t=0;t<N;t++){
		/* first row */
		if(t == 0){
			(*Ax)(t) = x(t) - x(t+1);
		}
		/* common row */
		if(t > 0 && t < N-1){
			(*Ax)(t) = -x(t-1) + 2.0*x(t) - x(t+1);
		}
		/* last row */
		if(t == N-1){
			(*Ax)(t) = -x(t-1) + x(t);
		}
	}
}

/* overloaded functions */
void get_dot(Scalar *xx, GammaVector<Scalar> x, GammaVector<Scalar> y, double *time_to_add){
	timer.start(); 

	*xx = get_dot(x,y);

	(*time_to_add) += timer.stop();	
}

void get_Ax_laplace(GammaVector<Scalar> *Ax, GammaVector<Scalar> x, double *time_to_add){
	timer.start(); 

	get_Ax_laplace(Ax,x);

	(*time_to_add) += timer.stop();	
}


Scalar get_dot(GammaVector<Scalar> x, GammaVector<Scalar> y){
	Scalar xx;

	get_dot(&xx, x, y);

	return xx;
}

Scalar get_dot(GammaVector<Scalar> x, GammaVector<Scalar> y, double *time_to_add){
	Scalar xx;

	get_dot(&xx, x, y,time_to_add);

	return xx;
}
