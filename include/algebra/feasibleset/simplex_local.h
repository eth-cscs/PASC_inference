/** @file simplex_local.h
 *  @brief simplex feasible set for petscvector with local problems
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_SIMPLEXFEASIBLESET_LOCAL_H
#define	PASC_SIMPLEXFEASIBLESET_LOCAL_H

#include "pascinference.h"

#ifndef USE_PETSC
 #error 'SIMPLEXFEASIBLESET_LOCAL is for PETSC'
#endif

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
	private:
		
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

		#ifdef USE_CUDA
			double *x_sorted; /**< for manipulation with sorted data on GPU */
			int blockSize; /**< block size returned by the launch configurator */
			int minGridSize; /**< the minimum grid size needed to achieve the maximum occupancy for a full device launch */
			int gridSize; /**< the actual grid size needed, based on input size */
		#endif

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
		void project(GeneralVector<PetscVector> &x);
		
};

/* cuda kernels cannot be a member of class */
#ifdef USE_CUDA
__device__ void device_sort_bubble(double *x_sorted, int t, int T, int K);

/** @brief kernel projection onto simplex
 *
 * computes a projection of a point onto simplex in nD
 *
 * take K-dimensional vector x[tK,tK+1,...tK+(K-1)] =: p
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
 * @param x_sorted allocated vector for manipulating with sorted x
 * @param Tlocal parameter of vector length
 * @param K parameter of vector length
 */ 
__global__ void kernel_project(double *x, double *x_sorted, int T, int K);
#endif


}
} /* end of namespace */


#ifdef USE_PETSC
 #include "external/petsc/algebra/feasibleset/simplex_local.h"
#endif



#endif
