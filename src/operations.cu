#include "operations.h"

/* main functions */
void get_dot(Scalar *xx, GammaVector<Scalar> x, GammaVector<Scalar> y){
	*xx = dot(x,y);
}

Scalar get_dot(GammaVector<Scalar> x, GammaVector<Scalar> y){
	Scalar xx;

	get_dot(&xx, x, y);

	return xx;
}

void get_Ax_laplace(GammaVector<Scalar>& Ax, GammaVector<Scalar> x){
	get_Ax_laplace(Ax,x,1.0);
}

void get_Ax_laplace(GammaVector<Scalar>& Ax, GammaVector<Scalar> x, Scalar alpha){
	int N = x.size();

	Ax(1,N-2) = 2*x(1,N-2) - x(0,N-3) - x(2,N-1);
	
	/* first and last */
	Ax(0) = x(0) - x(1);
	Ax(N-1) = x(N-1) - x(N-2);

	Ax *= alpha;
}
