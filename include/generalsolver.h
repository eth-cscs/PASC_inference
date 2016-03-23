/** @file generalsolver.h
 *  @brief class for manipulation with solver
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

/** @class GeneralSolverSetting
 *  @brief manipulation with settings
 *
 *  Parent class for manipulation with settings of solvers.
 *  All specific settings of solvers implementations should be defined as inherited classes from this class.
 * 
 */ 
class GeneralSolverSetting : public GeneralSetting {
	protected:
		
	public:
		int debug_mode; /**< print info about the progress */
		int maxit; /**< max number of iterations */
		double eps; /**< precision */

		/** @brief set default setting values
		 *
		 */ 
		GeneralSolverSetting() {};

		/** @brief destroy settings
		 *
		 */ 
		~GeneralSolverSetting() {};
		
		/** @brief print setting values
		 *
		 */ 
		virtual void print(std::ostream &output) const {};
};

/** @class GeneralSolver
 *  @brief solver for solving problems
 * 
 *  Parent class for manipulation with solvers.
 *  All specific solver implementations should be defined as inherited classes from this class.
 * 
 */ 
class GeneralSolver {
	protected:
		GeneralData *data; /**< pointer to data on which the solver operates */
	public:
		GeneralSolverSetting setting; /**< settings of the solver */

		/** @brief default constructor
		 * 
		 * Set inner data to null.
		 * 
		 */ 
		GeneralSolver() {
			data = NULL;
		};
		
		/** @brief constructor from data
		 * 
		 * Set inner pointer to data on which the solver operates to given data.
		 * 
		 * @param new_data input data for solver
		 */ 
		GeneralSolver(GeneralData &new_data) {
			data = &new_data;
		};

		/** @brief general destructor
		 * 
		 */ 
		~GeneralSolver() {};

		/** @brief print basic info about solver
		 * 
		 *  Print the name of the solver and call "print" method on inner members of solver.
		 * 
		 * @param output where to print
		 */ 
		virtual void print(std::ostream &output) const;

		/** @brief print status of solver
		 * 
		 *  Print the values of inner couters of the solver (for example iteration counter).
		 * 
		 * @param output where to print
		 */ 
		virtual void printstatus(std::ostream &output) const;

		/** @brief print timers of solver
		 * 
		 *  Print the sum values of inner operation timers.
		 * 
		 * @param output where to print
		 */ 
		virtual void printtimer(std::ostream &output) const;

		/** @brief get the name of solver
		 */ 
		virtual std::string get_name() const;

		/** @brief solve the problem
		 * 
		 * @todo can be included in solve(solvertype)
		 */
		virtual void solve() { this->solve(SOLVER_AUTO); };

		/** @brief solve the problem with given solver type
		 * 
		 * @param solvertype the type of solver
		 */
		virtual void solve(SolverType solvertype) {};

		friend std::ostream &operator<<(std::ostream &output, const GeneralSolver &solver); /* cannot be virtual, therefore it call virtual print() */

		/** @brief get the pointer to inner data
		 * 
		 */ 
		GeneralData *get_data() const;
};

/* general print, call virtual print() */
std::ostream &operator<<(std::ostream &output, const GeneralSolver &solver){
	if(DEBUG_MODE >= 100) coutMaster << "(GeneralSolver)OPERATOR: <<" << std::endl;
	output << solver.get_name();
	return output;
}

void GeneralSolver::print(std::ostream &output) const {
	output <<  this->get_name() << std::endl;
}

void GeneralSolver::printstatus(std::ostream &output) const {
	output <<  this->get_name() << std::endl;
}

void GeneralSolver::printtimer(std::ostream &output) const {
	output <<  this->get_name() << std::endl;
}

std::string GeneralSolver::get_name() const {
	return "GeneralSolver";
}

GeneralData *GeneralSolver::get_data() const{
	return this->data;
}



} /* end of namespace */


#endif
