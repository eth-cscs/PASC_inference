#ifndef PASCGENERALSOLVER_H
#define	PASCGENERALSOLVER_H

#include "common.h"
#include "algebra.h"

namespace pascinference {

class GeneralSolver {
	protected:

	public:
		GeneralSolver() {};
		~GeneralSolver() {};

		virtual void print(std::ostream &output) const {};

		friend std::ostream &operator<<(std::ostream &output, const GeneralSolver &solver); /* cannot be virtual, therefore it call virtual print() */
	
};

/* general print, call virtual print() */
std::ostream &operator<<(std::ostream &output, const GeneralSolver &solver){
	if(DEBUG_MODE >= 100) std::cout << "(GeneralSolver)OPERATOR: <<" << std::endl;
	solver.print(output);
	return output;
}



} /* end of namespace */


#endif
