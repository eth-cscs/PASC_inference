/** @file simplex_local.h
 *  @brief simplex feasible set for vector with local problems
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_SIMPLEXFEASIBLESET_LOCAL_H
#define	PASC_SIMPLEXFEASIBLESET_LOCAL_H

#include "general/algebra/feasibleset/generalfeasibleset.h"

namespace pascinference {
namespace algebra {

/** \class SimplexFeasibleSet_Local
 *  \brief class for manipulation with simplex set
 *
 *  Provides structure for manipulation with simplex feasible set, i.e.
 * \f[
 * \Omega = 
 *  \left\lbrace x \in R^{KT}: 
 *    \forall t = 0,\dots,T-1: \sum\limits_{k=0}^{K-1} x_{tK+k} = 1,  
 *    x \geq 0
 *  \right\rbrace
 *	\f]
 * 
*/
template<class VectorBase>
class SimplexFeasibleSet_Local: public GeneralFeasibleSet<VectorBase> {
	public:
		class ExternalContent;

	protected:
		friend class ExternalContent;
		ExternalContent *externalcontent;			/**< for manipulation with external-specific stuff */

		/** @brief sort array using bubble sort
		 * 
		 * @param x array with values
		 * @param n size of array
		*/ 		
		void sort_bubble(double *x, int n);

		/** @brief projection onto simplex
		 *
		 * computes a projection of a point onto simplex in nD
		 *
		 * take K-dimensional vector x[tK,tK+1,...tK+K-1] =: p
		 * and compute projection
		 * P(p) = arg min || p - y ||_2
		 * subject to constraints (which define simplex)
		 * y_0 + ... y_{K-1} = 1
		 * y_i >= 0 for all i=0,...,K-1
		 *
		 * in practical applications K is much more lower number than T
		 * K - number of clusters (2 - 10^2)
		 * T - length of time-series (10^5 - 10^9) 
		 * 
		 * @param x values of whole vector in array
		 * @param t where my subvector starts
		 * @param T number of local disjoint simplex subsets
		 * @param K size of subvector
		*/ 		
		void project_sub(double *x, int t, int T, int K);

		int T; /**< number of local disjoint simplex subsets */
		int K; /**< size of each simplex subset */
				
	public:
		/** @brief default constructor
		*/ 	
		SimplexFeasibleSet_Local(int T, int K);
		
		/** @brief default destructor
		 */ 
		~SimplexFeasibleSet_Local();

		/** @brief print properties of this feasible set
		 * 
		 * @param output where to print
		 */ 
		void print(ConsoleOutput &output) const;

		/** @brief get name of this feasible set
		 */
		virtual std::string get_name() const;

		/** @brief compute projection onto feasible set
		 * 
		 * @param x point which will be projected
		 */		
		void project(GeneralVector<VectorBase> &x);

		ExternalContent *get_externalcontent() const;		
	
};


}
} /* end of namespace */


/* ---------------- implementation ---------- */
namespace pascinference {
namespace algebra {

/* constructor */
template<class VectorBase>
SimplexFeasibleSet_Local<VectorBase>::SimplexFeasibleSet_Local(int T, int K){
	LOG_FUNC_BEGIN

	/* set initial content */
	this->T = T;
	this->K = K;

	LOG_FUNC_END
}

/* general destructor */
template<class VectorBase>
SimplexFeasibleSet_Local<VectorBase>::~SimplexFeasibleSet_Local(){
	LOG_FUNC_BEGIN
	
	LOG_FUNC_END	
}

/* print info about feasible set */
template<class VectorBase>
void SimplexFeasibleSet_Local<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN
	
	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - nmb of subsets:     " << T << std::endl;
	output <<  " - size of subset:     " << K << std::endl;

	LOG_FUNC_END
}

template<class VectorBase>
std::string SimplexFeasibleSet_Local<VectorBase>::get_name() const {
	return "SimplexFeasibleSet_Local";
}

template<class VectorBase>
void SimplexFeasibleSet_Local<VectorBase>::project(GeneralVector<VectorBase> &x) {
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

template<class VectorBase>
void SimplexFeasibleSet_Local<VectorBase>::sort_bubble(double *x, int n){
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

template<class VectorBase>
void SimplexFeasibleSet_Local<VectorBase>::project_sub(double *x, int t, int T, int K){
	if(t<T){ /* maybe we call more than T kernels */
		int k;

		bool is_inside = true;
		double sum = 0.0;
	
		/* control inequality constraints */
		for(k = 0; k < K; k++){ // TODO: could be performed parallely  
			if(x[t*K+k] < 0.0){
				is_inside = false;
			}
			sum += x[t*K+k];
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
				y[k] = x[t*K+k]; 
			}
			sort_bubble(y,K);

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
				ti = x[t*K+k] - t_hat;	
				if(ti > 0.0){
					x[t*K+k] = ti;
				} else {
					x[t*K+k] = 0.0;
				}
			}
			
			delete y;
		}
		
	}

	/* if t >= T then relax and do nothing */	
}


}
} /* end of namespace */



#endif
