#include "operations.h"

/* main functions */
void get_dot(Scalar *xx, GammaVector x, GammaVector y){
	*xx = dot(x,y);
}

Scalar get_dot(GammaVector x, GammaVector y){
	Scalar xx;

	get_dot(&xx, x, y);

	return xx;
}

void get_Ax_laplace(GammaVector& Ax, GammaVector x, int K){
	get_Ax_laplace(Ax,x,K,1.0);
}

void get_Ax_laplace(GammaVector& Ax, GammaVector x, int K, Scalar alpha){
	int N = x.size();
	int T = N/(double)K;
	int k;

	Ax(1,N-2) = 2*x(1,N-2) - x(0,N-3) - x(2,N-1);
	
	/* first and last in each block */
	for(k=0;k<K;k++){
		Ax(k*T) = x(k*T) - x(k*T+1);
		Ax((k+1)*T-1) = x((k+1)*T-1) - x((k+1)*T-2);
	}
	
	Ax *= alpha;
}
