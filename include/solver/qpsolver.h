#ifndef QPSOLVER_H
#define	QPSOLVER_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "generalsolver.h"
#include "algebra.h"

#include "data/qpdata.h"


namespace pascinference {

/* QPSolver */ 
template<class VectorBase>
class QPSolver: public GeneralSolver {
	private:
		const QPData<VectorBase> *data; /* data for QP problem */

	public:
		QPSolver();
		~QPSolver();

		void print(std::ostream &output) const;


};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
QPSolver<VectorBase>::QPSolver(){
	if(DEBUG_MODE >= 100) std::cout << "(QPSolver)CONSTRUCTOR" << std::endl;

	/* set initial content */

}

/* destructor */
template<class VectorBase>
QPSolver<VectorBase>::~QPData(){
	if(DEBUG_MODE >= 100) std::cout << "(QPSolver)DESTRUCTOR" << std::endl;
	
}

/* print info about problem */
template<class VectorBase>
void QPSolver<VectorBase>::print(std::ostream &output) const {
	output << " QPSolver" << std::endl;
	
	/* give information about presence of the data */

		
}


} /* end namespace */

#endif
