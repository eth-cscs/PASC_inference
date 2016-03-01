#ifndef PASCGENERALSOLVER_H
#define	PASCGENERALSOLVER_H

#include "common.h"
#include "algebra.h"
#include "solver/list.h"
#include "generalsetting.h"
#include "generaldata.h"
#include "generalresult.h"

namespace pascinference {

/* setting class */
class GeneralSolverSetting : public GeneralSetting {
	protected:
	
	public:
		GeneralSolverSetting() {};
		~GeneralSolverSetting() {};

		virtual void print(std::ostream &output) const {

		};
	
};

/* solver class */
class GeneralSolver {
	protected:
		const GeneralData *data; /* pointer to data on which the solver operates */
		const GeneralResult *result; /* here the solver stores the results */
	public:
		GeneralSolverSetting setting; // TODO: private?

		GeneralSolver() {
			data = NULL;
			result = NULL;
		};
		GeneralSolver(const GeneralData &new_data, const GeneralResult &new_result) {
			data = &new_data;
			result = &new_result;
		};
		~GeneralSolver() {};

		virtual void print(std::ostream &output) const {};

		virtual void solve() { this->solve(SOLVER_AUTO); };
		virtual void solve(SolverType solvertype) {};

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
