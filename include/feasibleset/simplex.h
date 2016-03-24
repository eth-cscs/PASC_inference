/** @file simplex.h
 *  @brief simplex feasible set
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_SIMPLEXFEASIBLESET_H
#define	PASC_SIMPLEXFEASIBLESET_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "generalfeasibleset.h"
#include "algebra.h"

#ifdef USE_PETSCVECTOR
	typedef petscvector::PetscVector GlobalPetscVector;
#endif

#ifdef USE_MINLIN
	typedef minlin::threx::HostVector<double> HostMinLinVector;
	typedef minlin::threx::DeviceVector<double> DeviceMinLinVector;
#endif

namespace pascinference {


/** \class SimplexFeasibleSet
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
template<class VectorBase>
class SimplexFeasibleSet: public GeneralFeasibleSet<VectorBase> {
	private:
		
		/** @brief sort array using bubble sort
		 * 
		 * @param x array with values
		 * @param n size of array
		*/ 		
		void sort_bubble(double *x, int n);

		/** @brief compute projection to one simplex subset
		 * 
		 * @param x_sub values of subvector in array
		 * @param n size of subvector
		*/ 		
		void get_projection_sub(double *x_sub, int n);
		
		#ifdef USE_PETSCVECTOR
			bool petsc_projection_init; /**< if the initialization of projection was not performed, then = false */
			IS *petsc_projection_is; /**< array of indexsets with coeficient indexes */
			int petsc_projection_Townership_low, petsc_projection_Townership_high;
		#endif

	public:
		/** @brief default constructor
		*/ 	
		SimplexFeasibleSet(int T, int K);
		
		/** @brief default destructor
		 */ 
		~SimplexFeasibleSet();

		/** @brief print properties of this feasible set
		 * 
		 * @param output where to print
		 */ 
		void print(std::ostream &output) const;

		/** @brief get name of this feasible set
		 */
		virtual std::string get_name() const;

		/* variables */
		int T; /**< number of disjoint simplex subsets */
		int K; /**< size of each simplex subset */
		
		/** @brief compute projection onto feasible set
		 * 
		 * @param x point which will be projected
		 */		
		void project(GeneralVector<VectorBase> &x);


};

/* sorry, __device__ and __global__ functions cannot be members of template class */
#ifdef USE_GPU
__device__ void SimplexFeasibleSet_device_sort_bubble(double *x, int n);
__global__ void SimplexFeasibleSet_kernel_get_projection_sub(double *x, int T, int K);
#endif

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
SimplexFeasibleSet<VectorBase>::SimplexFeasibleSet(int Tnew, int Knew){
	if(DEBUG_MODE >= 100) coutMaster << "(SimplexFeasibleSet)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->T = Tnew;
	this->K = Knew;
	
	this->petsc_projection_init = false;
}

/* general destructor */
template<class VectorBase>
SimplexFeasibleSet<VectorBase>::~SimplexFeasibleSet(){
	if(DEBUG_MODE >= 100) coutMaster << "(SimplexFeasibleSet)DESTRUCTOR" << std::endl;
	
}

/* print info about feasible set */
template<class VectorBase>
void SimplexFeasibleSet<VectorBase>::print(std::ostream &output) const {
	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - T:     " << T << std::endl;
	output <<  " - K:     " << K << std::endl;
		
}

template<class VectorBase>
std::string SimplexFeasibleSet<VectorBase>::get_name() const {
	return "SimplexFeasibleSet";
}



/* -------- minlin::threx::HostVector ---------- */
#ifdef USE_MINLIN
template<>
void SimplexFeasibleSet<HostMinLinVector>::project(GeneralVector<HostMinLinVector> &x) {
	if(DEBUG_MODE >= 100) coutMaster << "(SimplexFeasibleSet)FUNCTION: project HostMinLinVector" << std::endl;

	int t,k;
	double x_sub[K];  /* GammaVector x_sub(K); */

	for(t=0;t<T;t++){
		/* cut x_sub from x */
		for(k=0;k<K;k++){
			x_sub[k] = x(k*T+t);
		}
		
		/* compute subprojection */
		get_projection_sub(x_sub, this->K);

		/* add x_sub back to x */
		for(k=0;k<K;k++){
			x(k*T+t) = x_sub[k];
		}
	}

}
#endif

/* project x_sub to feasible set defined by equality and inequality constraints
 * sum(x_sub) = 1
 * x_sub >= 0
 */
template<class VectorBase>
void SimplexFeasibleSet<VectorBase>::get_projection_sub(double *x_sub, int n){
	int i;

	bool is_inside = true;
	double sum = 0.0;
	
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
		double y[n], sum_y;
		for(i=0;i<n;i++){
			y[i] = x_sub[i]; 
		}
		sort_bubble(y,n);

		/* now perform analytical solution of projection problem */	
		double t_hat = 0.0;
		i = n - 1;
		double ti;

		while(i >= 1){
			/* compute sum(y) */
			sum_y = 0.0;
			for(j=i;j<n;j++){ /* sum(y(i,n-1)) */
				sum_y += y[j];
			}
				
			ti = (sum_y - 1.0)/(double)(n-i);
			if(ti >= y[i-1]){
				t_hat = ti;
				i = -1; /* break */
			} else {
				i = i - 1;
			}
		}

		if(i == 0){
			t_hat = (sum-1.0)/(double)n; /* uses sum=sum(x_sub) */
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
template<class VectorBase>
void SimplexFeasibleSet<VectorBase>::sort_bubble(double *x, int n){
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



/* -------- minlin::threx::DeviceVector ---------- */
#ifdef USE_GPU
template<>
void SimplexFeasibleSet<DeviceMinLinVector>::project(GeneralVector<DeviceMinLinVector> &x) {

	/* call projection using kernel */
	double *xp = x.pointer();
	
	// TODO: compute optimal nmb of threads/kernels
	SimplexFeasibleSet_kernel_get_projection_sub<<<this->T, 1>>>(xp,this->T,this->K);
	
	/* synchronize kernels, if there is an error with cuda, then it will appear here */ 
	gpuErrchk( cudaDeviceSynchronize() );	
	

}


__device__
void SimplexFeasibleSet_device_sort_bubble(double *x, int n){
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

__global__
void SimplexFeasibleSet_kernel_get_projection_sub(double *x, int T, int K){
	/* compute my id */
	int t = blockIdx.x*blockDim.x + threadIdx.x;

	if(t<T){
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
			SimplexFeasibleSet_device_sort_bubble(y,K);

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

	/* if t >= N then relax and do nothing */	

}

#endif



/* -------- petscvector::PetscVector ---------- */
#ifdef USE_PETSCVECTOR
template<>
void SimplexFeasibleSet<GlobalPetscVector>::project(GeneralVector<GlobalPetscVector> &x) {

	int i;

	/* initialization - how much I will compute? */
	if(!this->petsc_projection_init){
		
		Vec layout;

		/* try to make a global vector of length T and then get the indexes of begin and end of local portion */
		TRY( VecCreate(PETSC_COMM_WORLD,&layout) );
		TRY( VecSetSizes(layout,PETSC_DECIDE,this->T) );
		TRY( VecSetFromOptions(layout) );

		/* get the ownership range - now I know how much I will calculate from the time-series */
		TRY( VecGetOwnershipRange(layout, &(this->petsc_projection_Townership_low), &(this->petsc_projection_Townership_high)) );

		/* destroy testing vector - it is useless now */
		TRY(VecDestroy(&layout));

		/* create array index set of "my" indeces */
		TRY( PetscMalloc((this->petsc_projection_Townership_high - this->petsc_projection_Townership_low)*sizeof(IS),&this->petsc_projection_is) );
		for(i =0; i<this->petsc_projection_Townership_high - this->petsc_projection_Townership_low;i++){
			/* prepare index set [low, low + T, ... , low + (K-1)*T ] */
			TRY( ISCreateStride(PETSC_COMM_SELF, this->K, this->petsc_projection_Townership_low + i, this->T, &(this->petsc_projection_is[i])) );
		}

		/* projection was initialized */
		this->petsc_projection_init = true;

	}

	if(DEBUG_MODE >= 100){
		coutAll << " my ownership: [" << this->petsc_projection_Townership_low << ", " << this->petsc_projection_Townership_high << "]" << std::endl;
	}

	Vec x_sub;
	double *x_sub_arr;
	
	Vec x_vec = x.get_vector();

	/* go throught local portion of time-serie and perform the projection */
	for(i = 0; i < this->petsc_projection_Townership_high - this->petsc_projection_Townership_low; i++){

		/* get the subvector from global vector */
		TRY( VecGetSubVector(x_vec,petsc_projection_is[i], &x_sub) );

		/* get the array */
		TRY( VecGetArray(x_sub, &x_sub_arr) );

		/* perform the projection on this subvector */
		get_projection_sub(x_sub_arr, this->K);

		/* print the array of subvector */
		if(DEBUG_MODE >= 100){
			int j;
			coutAll << " xsub_" << this->petsc_projection_Townership_low+i << " = [ ";
			for(j=0;j<this->K;j++){
				coutAll << x_sub_arr[j];
				if(j < this->K-1) coutAll << ", ";
			}
			coutAll << " ]" << std::endl;
		}

		/* restore the array */
		TRY( VecRestoreArray(x_sub, &x_sub_arr) );

		TRY( VecRestoreSubVector(x_vec,petsc_projection_is[i], &x_sub) );

		//TODO: deal with x_sub

	}

	x.valuesUpdate();
	
}

/* petsc specific destructor */
template<>
SimplexFeasibleSet<GlobalPetscVector>::~SimplexFeasibleSet(){
	if(DEBUG_MODE >= 100) coutMaster << "(SimplexFeasibleSet)DESTRUCTOR" << std::endl;

	/* destruction of index sets */
	if(this->petsc_projection_init){
		int i;
		for(i =0; i<this->petsc_projection_Townership_high - this->petsc_projection_Townership_low;i++){
			/* prepare index set [low, low + T, ... , low + (K-1)*T ] */
			TRY( ISDestroy(&(this->petsc_projection_is[i])) );
		}

		/* projection was initialized */
		this->petsc_projection_init = false;

	}
	
}

#endif



} /* end namespace */

#endif
