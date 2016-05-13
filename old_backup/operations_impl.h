
namespace pascinference {

/* ----------------- minlin::threx::HostVector ------------------- */

/* return the dot product */
void get_dot(Scalar *xx, minlin::threx::HostVector<Scalar> x, minlin::threx::HostVector<Scalar> y){
	*xx = dot(x,y); /* Minlin operation */
}

Scalar get_dot(minlin::threx::HostVector<Scalar> x, minlin::threx::HostVector<Scalar> y){
	Scalar xx;

	get_dot(&xx, x, y);

	return xx;
}

void get_Ax_laplace(minlin::threx::HostVector<Scalar> &Ax, minlin::threx::HostVector<Scalar> x, int K){
	get_Ax_laplace(Ax,x,K,1.0);
}

/* return matrix-vector multiplication - MinLin operations */
void get_Ax_laplace(minlin::threx::HostVector<Scalar>& Ax, minlin::threx::HostVector<Scalar> x, int K, Scalar alpha){
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


/* ----------------- minlin::threx::DeviceVector ------------------- */

#ifdef USE_GPU

	/* return the dot product */
	void get_dot(Scalar *xx, minlin::threx::DeviceVector<Scalar> x, minlin::threx::DeviceVector<Scalar> y){
		*xx = dot(x,y); /* Minlin operation */
	}

	Scalar get_dot(minlin::threx::DeviceVector<Scalar> x, minlin::threx::DeviceVector<Scalar> y){
		Scalar xx;

		get_dot(&xx, x, y);

		return xx;
	}

	void get_Ax_laplace(minlin::threx::DeviceVector<Scalar> &Ax, minlin::threx::DeviceVector<Scalar> x, int K){
		get_Ax_laplace(Ax,x,K,1.0);
	}

	/* return matrix-vector multiplication - MinLin operations */
	void get_Ax_laplace(minlin::threx::DeviceVector<Scalar>& Ax, minlin::threx::DeviceVector<Scalar> x, int K, Scalar alpha){
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

/* ----------------- petscvector::PetscVector ------------------- */

#ifdef USE_PETSC

	/* return the dot product */
	void get_dot(Scalar *xx, petscvector::PetscVector &x, petscvector::PetscVector &y){
		*xx = dot(x,y); /* Petsc operation */
	}

	Scalar get_dot(petscvector::PetscVector &x, petscvector::PetscVector &y){
		Scalar xx;

		get_dot(&xx, x, y);

		return xx;
	}

	void get_Ax_laplace(petscvector::PetscVector &Ax, petscvector::PetscVector &x, int K){
		get_Ax_laplace(Ax,x,K,1.0);
	}

	/* return matrix-vector multiplication - Petsc operations */
	void get_Ax_laplace(petscvector::PetscVector& Ax, petscvector::PetscVector &x, int K, Scalar alpha){
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

} /* end of namespace */
