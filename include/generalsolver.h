/** @file generalsolver.h
 *  @brief class for manipulation with solver
 *
 *  Defines the parent class for manipulation with solvers - setting data, algorithm for solving the problem, return solution, etc.
 *  All specific solvers implementations should be defined as inherited classes from this class.
 *
 *  @author Lukas Pospisil
 */

#ifndef PASCGENERALSOLVER_H
#define	PASCGENERALSOLVER_H

#include "common.h"
#include "algebra.h"
#include "solver/list.h"
#include "generalsetting.h"
#include "generaldata.h"

namespace pascinference {

/* setting class */
class GeneralSolverSetting : public GeneralSetting {
	protected:
		
	public:
		int debug_mode; /**< print info about the progress */
		int maxit; /**< max number of iterations */
		double eps; /**< precision */


		GeneralSolverSetting() {};
		~GeneralSolverSetting() {};

		virtual void print(std::ostream &output) const {};
};

/* solver class */
class GeneralSolver {
	protected:
		GeneralData *data; /* pointer to data on which the solver operates */
	public:
		GeneralSolverSetting setting; // TODO: private?

		GeneralSolver() {
			data = NULL;
		};
		GeneralSolver(GeneralData &new_data) {
			data = &new_data;
		};
		~GeneralSolver() {};

		virtual void print(std::ostream &output) const;
		virtual void printstatus(std::ostream &output) const;
		virtual void printtimer(std::ostream &output) const;
		virtual std::string get_name() const;

		virtual void solve() { this->solve(SOLVER_AUTO); };
		virtual void solve(SolverType solvertype) {};

		friend std::ostream &operator<<(std::ostream &output, const GeneralSolver &solver); /* cannot be virtual, therefore it call virtual print() */
	
		GeneralData *get_data() const;
};

/* general print, call virtual print() */
std::ostream &operator<<(std::ostream &output, const GeneralSolver &solver){
	if(DEBUG_MODE >= 100) coutMaster << offset <<"(GeneralSolver)OPERATOR: <<" << std::endl;
	output << solver.get_name();
	return output;
}

void GeneralSolver::print(std::ostream &output) const {
	output << offset << this->get_name() << std::endl;
}

void GeneralSolver::printstatus(std::ostream &output) const {
	output << offset << this->get_name() << std::endl;
}

void GeneralSolver::printtimer(std::ostream &output) const {
	output << offset << this->get_name() << std::endl;
}

std::string GeneralSolver::get_name() const {
	return "GeneralSolver";
}

GeneralData *GeneralSolver::get_data() const{
	return this->data;
}



} /* end of namespace */


#endif
