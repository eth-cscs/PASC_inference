/** @file simplex_lineqbound.h
 *  @brief simplex feasible set defined by linear equality constraints and bound constraints
 *
 *  @author Lukas Pospisil
 */

#ifndef PASC_SIMPLEX_LINEQBOUNDFEASIBLESET_H
#define	PASC_SIMPLEX_LINEQBOUNDFEASIBLESET_H

#include "general/algebra/feasibleset/generalfeasibleset.h"


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
		class ExternalContent;
		friend class ExternalContent;
		ExternalContent *externalcontent;			/**< for manipulation with external-specific stuff */
		
		int T; /**< number of global disjoint simplex subsets */
		int Tlocal; /**< number of local disjoint simplex subsets */
		int K; /**< size of each simplex subset */

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
		void project(GeneralVector<VectorBase> &x);
		
//		SimplexFeasibleSet_LinEqBound<VectorBase>::ExternalContent *get_externalcontent() const;
		
};


}
} /* end of namespace */


/* ------ IMPLEMENTATION ------- */
namespace pascinference {
namespace algebra {

/* constructor */
template<class VectorBase>
SimplexFeasibleSet_LinEqBound<VectorBase>::SimplexFeasibleSet_LinEqBound(int T, int Tlocal, int K){
	LOG_FUNC_BEGIN

	this->T = T;
	this->Tlocal = Tlocal;
	this->K = K;

	//TODO

	LOG_FUNC_END
}

/* general destructor */
template<class VectorBase>
SimplexFeasibleSet_LinEqBound<VectorBase>::~SimplexFeasibleSet_LinEqBound(){
	LOG_FUNC_BEGIN
	
	
	LOG_FUNC_END	
}

/* print info about feasible set */
template<class VectorBase>
void SimplexFeasibleSet_LinEqBound<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN
	
	output <<  this->get_name() << std::endl;
	
	/* give information about presence of the data */
	output <<  " - nmb of subsets:     " << T << std::endl;
	output <<  " - size of subset:     " << K << std::endl;


	LOG_FUNC_END
}

template<class VectorBase>
std::string SimplexFeasibleSet_LinEqBound<VectorBase>::get_name() const {
	return "SimplexFeasibleSet_LinEqBound";
}

template<class VectorBase>
void SimplexFeasibleSet_LinEqBound<VectorBase>::project(GeneralVector<VectorBase> &x) {
	LOG_FUNC_BEGIN
	
	//TODO: give error, the projection is not defined for this type of feasible set

	LOG_FUNC_END
}

/* define blank external content for general VectorBase */
template<class VectorBase>
class SimplexFeasibleSet_LinEqBound<VectorBase>::ExternalContent {
};

/*SimplexFeasibleSet_LinEqBound<VectorBase>::ExternalContent *SimplexFeasibleSet_LinEqBound<VectorBase>::get_externalcontent() const {
	return externalcontent;
}
*/

}
} 

#endif
