#ifndef PASCGENERALPROBLEM_H
#define	PASCGENERALPROBLEM_H

#include "common.h"
#include "generaldata.h"
#include "generalsolver.h"

namespace pascinference {

class GeneralProblem {
	protected:
		bool solved; /* is the problem solved ot not? */
		Timer timer_total; /* from init to finalize */
		
		const GeneralData *data;
		const GeneralSolver *solver;
		const GeneralResult *result;
		
	public:
		GeneralProblem() {};
		~GeneralProblem() {};
	
		virtual void solve() {};
		virtual void print(std::ostream &output) const {};

		friend std::ostream &operator<<(std::ostream &output, const GeneralProblem &problem); /* cannot be virtual, therefore it call virtual print() */
};


/* general print, call virtual print() */
std::ostream &operator<<(std::ostream &output, const GeneralProblem &problem){
	if(DEBUG_MODE >= 100) coutMaster << offset <<"(GeneralProblem)OPERATOR: <<" << std::endl;
	problem.print(output);
	return output;
}


} /* end of namespace */



#endif
