#ifndef SIMPLEXFEASIBLESET_H
#define	SIMPLEXFEASIBLESET_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "generalfeasibleset.h"
#include "algebra.h"

namespace pascinference {

/* Simplex Feasible Set */
/* \Omega = 
 *  \lbrace x \in \mathbb{R}^{KT}: 
 *    \sum\limits_{k=0}^{K-1} x_{t+kT} = 1, \forall t = 0,\dots,T-1 
 *    x \geq 0
 *  \rbrace
 */  
template<class VectorBase>
class SimplexFeasibleSet: public GeneralFeasibleSet<VectorBase> {
	public:
		SimplexFeasibleSet();
		~SimplexFeasibleSet();

		void print(std::ostream &output) const;

		/* variables */
		int T; 
		int K; 
		
		void project(GeneralVector<VectorBase> &x);


};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
SimplexFeasibleSet<VectorBase>::SimplexFeasibleSet(){
	if(DEBUG_MODE >= 100) std::cout << "(SimplexFeasibleSet)CONSTRUCTOR" << std::endl;

	/* set initial content */
	this->T = 0;
	this->K = 0;
}

/* destructor */
template<class VectorBase>
SimplexFeasibleSet<VectorBase>::~SimplexFeasibleSet(){
	if(DEBUG_MODE >= 100) std::cout << "(SimplexFeasibleSet)DESTRUCTOR" << std::endl;
	
}

/* print info about feasible set */
template<class VectorBase>
void SimplexFeasibleSet<VectorBase>::print(std::ostream &output) const {
	output << " SimplexFeasibleSet" << std::endl;
	
	/* give information about presence of the data */
	output << "  - T:     " << T << std::endl;
	output << "  - K:     " << K << std::endl;
		
}

template<class VectorBase>
void SimplexFeasibleSet<VectorBase>::project(GeneralVector<VectorBase> &x) {
	if(DEBUG_MODE >= 100) std::cout << "(SimplexFeasibleSet)FUNCTION: project" << std::endl;
	
	x(gall) = 3.1416;
		
}



} /* end namespace */

#endif
