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
#include "generaldata.h"

namespace pascinference {

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
		/* settings */
		int debug_mode; /**< print info about the progress */
		int maxit; /**< max number of iterations */
		double eps; /**< precision */

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
		virtual void print(ConsoleOutput &output) const;

		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;

		/** @brief print status of solver
		 * 
		 *  Print the values of inner couters of the solver (for example iteration counter).
		 * 
		 * @param output where to print
		 */ 
		virtual void printstatus(ConsoleOutput &output) const;

		/** @brief print content of all solver data
		 * 
		 * @param output where to print
		 */ 
		virtual void printcontent(ConsoleOutput &output) const;

		/** @brief print timers of solver
		 * 
		 *  Print the sum values of inner operation timers.
		 * 
		 * @param output where to print
		 */ 
		virtual void printtimer(ConsoleOutput &output) const;

		/** @brief get the name of solver
		 */ 
		virtual std::string get_name() const;

		/** @brief solve the problem
		 * 
		 */
		virtual void solve() {};

		friend ConsoleOutput &operator<<(ConsoleOutput &output, const GeneralSolver &solver); /* cannot be virtual, therefore it call virtual print() */

		/** @brief get the pointer to inner data
		 * 
		 */ 
		GeneralData *get_data() const;

		/** @brief return number of iterations
		 * 
		 */ 
		virtual int get_it() const;
};

/* general print, call virtual print() */
ConsoleOutput &operator<<(ConsoleOutput &output, const GeneralSolver &solver){
	if(DEBUG_MODE >= 100) coutMaster << "(GeneralSolver)OPERATOR: <<" << std::endl;
	output << solver.get_name();
	return output;
}

void GeneralSolver::print(ConsoleOutput &output) const {
	output << this->get_name() << std::endl;
}

void GeneralSolver::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	output_global << this->get_name() << std::endl;
	output_global.synchronize();
}

void GeneralSolver::printstatus(ConsoleOutput &output) const {
	output << this->get_name() << ": status" << std::endl;
}

void GeneralSolver::printcontent(ConsoleOutput &output) const {
	output << this->get_name() << ": content" << std::endl;
}

void GeneralSolver::printtimer(ConsoleOutput &output) const {
	output << this->get_name() << ": timer" << std::endl;
}

std::string GeneralSolver::get_name() const {
	return "GeneralSolver";
}

GeneralData *GeneralSolver::get_data() const{
	return this->data;
}

int GeneralSolver::get_it() const{
	return 0;
}

} /* end of namespace */


#endif
