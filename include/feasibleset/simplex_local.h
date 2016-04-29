/** @file simplex_global.h
 *  @brief simplex feasible set for petscvector with local problems
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_SIMPLEXFEASIBLESET_GLOBAL_H
#define	PASC_SIMPLEXFEASIBLESET_GLOBAL_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "generalfeasibleset.h"
#include "algebra.h"

#ifndef USE_PETSCVECTOR
 #error 'TSSOLVER_GLOBAL is for PETSCVECTOR'
#endif

typedef petscvector::PetscVector PetscVector;


namespace pascinference {


/** \class SimplexFeasibleSet_Local
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
class SimplexFeasibleSet_Local: public GeneralFeasibleSet<PetscVector> {
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

	public:
		/** @brief default constructor
		*/ 	
		SimplexFeasibleSet_Local(int T, int K_local);
		
		/** @brief default destructor
		 */ 
		~SimplexFeasibleSet_Local();

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
		int K_local; /**< size of each simplex subset */
		
		/** @brief compute projection onto feasible set
		 * 
		 * @param x point which will be projected
		 */		
		void project(GeneralVector<PetscVector> &x);

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
SimplexFeasibleSet_Local::SimplexFeasibleSet_Local(int Tnew, int Knew){
	if(DEBUG_MODE >= 100) coutMaster << "(SimplexFeasibleSet_Local)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->T = Tnew;
	this->K_local = Knew;

}

/* general destructor */
SimplexFeasibleSet_Local::~SimplexFeasibleSet_Local(){
	if(DEBUG_MODE >= 100) coutMaster << "(SimplexFeasibleSet_Local)DESTRUCTOR" << std::endl;
	
}

/* print info about feasible set */
void SimplexFeasibleSet_Local::print(std::ostream &output) const {
	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - T:       " << T << std::endl;
	output <<  " - K_local: " << K_local << std::endl;
		
}

std::string SimplexFeasibleSet_Local::get_name() const {
	return "SimplexFeasibleSet_Local";
}



void SimplexFeasibleSet_Local::project(GeneralVector<PetscVector> &x) {
	if(DEBUG_MODE >= 100) coutMaster << "(SimplexFeasibleSet_Local)FUNCTION: project" << std::endl;

	int t,k;
	double x_sub[K_local];  /* GammaVector x_sub(K); */

	/* get local array */
	double *x_arr;
	
	TRY( VecGetArray(x.get_vector(),&x_arr) );

	for(t=0;t<T;t++){
		/* cut x_sub from x */
		for(k=0;k<K_local;k++){
			x_sub[k] = x_arr[k*T+t];
		}
		
		/* compute subprojection */
		get_projection_sub(x_sub, this->K_local);

		/* add x_sub back to x */
		for(k=0;k<K_local;k++){
			x_arr[k*T+t] = x_sub[k];
		}
	}

	TRY( VecRestoreArray(x.get_vector(),&x_arr) );


}


/* project x_sub to feasible set defined by equality and inequality constraints
 * sum(x_sub) = 1
 * x_sub >= 0
 */
void SimplexFeasibleSet_Local::get_projection_sub(double *x_sub, int n){
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
void SimplexFeasibleSet_Local::sort_bubble(double *x, int n){
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



} /* end namespace */

#endif
