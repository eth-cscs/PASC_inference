/** @file simplex_lineqbound.h
 *  @brief simplex feasible set defined by linear equality constraints and bound constraints
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_SIMPLEX_LINEQBOUNDFEASIBLESET_H
#define	PASC_SIMPLEX_LINEQBOUNDFEASIBLESET_H

#include "algebra/feasibleset/generalfeasibleset.h"


namespace pascinference {
namespace algebra {

/** \class SimplexFeasibleSet_LinEqBound
 *  \brief simplex feasible set defined by linear equality constraints and bound constraints
 *
 *  Provides structure for manipulation with feasible set in form
 * \f[
 * \Omega = 
 *  \left\lbrace x \in R^{KT}: 
 *    \forall t = 0,\dots,T-1: \sum\limits_{k=0}^{K-1} x_{tK+k} = 1,  
 *    x \geq 0
 *  \right\rbrace
 *	\f] =
 *  \left\lbrace x \in R^{KT}: 
 *    Bx = c \wedge x \geq 0
 *  \right\rbrace
 * 
*/
template<class VectorBase>
class SimplexFeasibleSet_LinEqBound: public GeneralFeasibleSet<VectorBase> {
	private:
		
		int T; /**< number of global disjoint simplex subsets */
		int Tlocal; /**< number of local disjoint simplex subsets */
		int K; /**< size of each simplex subset */

		Mat B; /**< matrix of equality constraints Bx=c */
		Vec c; /**< vector of equality constraints Bx=c */
		Vec lb; /**< vector of lower bounds */
		
	public:
		/** @brief default constructor
		*/ 	
		SimplexFeasibleSet_LinEqBound(int T, int Tlocal, int K);
		
		/** @brief default destructor
		 */ 
		~SimplexFeasibleSet_LinEqBound();

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
		
		Mat get_B() const;
		Vec get_c() const;
		Vec get_lb() const;
		
};


}
} /* end of namespace */

#ifdef USE_PETSC
 #include "external/petsc/algebra/feasibleset/simplex_lineqbound.h"
#endif

#endif
