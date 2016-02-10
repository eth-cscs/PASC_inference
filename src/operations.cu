#include "operations.h"

/* ----------------- HostVector ------------------- */

/* return the dot product */
void get_dot(Scalar *xx, HostVector<Scalar> x, HostVector<Scalar> y){
	*xx = dot(x,y); /* Minlin operation */
}

Scalar get_dot(HostVector<Scalar> x, HostVector<Scalar> y){
	Scalar xx;

	get_dot(&xx, x, y);

	return xx;
}

void get_Ax_laplace(HostVector<Scalar> &Ax, HostVector<Scalar> x, int K){
	get_Ax_laplace(Ax,x,K,1.0);
}

/* return matrix-vector multiplication - MinLin operations */
void get_Ax_laplace(HostVector<Scalar>& Ax, HostVector<Scalar> x, int K, Scalar alpha){
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


/* ----------------- DeviceVector ------------------- */

#ifdef USE_GPU

	/* return the dot product */
	void get_dot(Scalar *xx, DeviceVector<Scalar> x, DeviceVector<Scalar> y){
		*xx = dot(x,y); /* Minlin operation */
	}

	Scalar get_dot(DeviceVector<Scalar> x, DeviceVector<Scalar> y){
		Scalar xx;

		get_dot(&xx, x, y);

		return xx;
	}

	void get_Ax_laplace(DeviceVector<Scalar> &Ax, DeviceVector<Scalar> x, int K){
		get_Ax_laplace(Ax,x,K,1.0);
	}

	/* return matrix-vector multiplication - MinLin operations */
	void get_Ax_laplace(DeviceVector<Scalar>& Ax, DeviceVector<Scalar> x, int K, Scalar alpha){
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

#endif

/* ----------------- PetscVector ------------------- */

#ifdef USE_PETSC

	/* return the dot product */
	void get_dot(Scalar *xx, PetscVector x, PetscVector y){
		*xx = dot(x,y); /* Petsc operation */
	}

	Scalar get_dot(PetscVector x, PetscVector y){
		Scalar xx;

		get_dot(&xx, x, y);

		return xx;
	}

	void get_Ax_laplace(PetscVector &Ax, PetscVector x, int K){
		get_Ax_laplace(Ax,x,K,1.0);
	}

	/* return matrix-vector multiplication - Petsc operations */
	void get_Ax_laplace(PetscVector& Ax, PetscVector x, int K, Scalar alpha){
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

#endif


