/** @file generalsolver.h
 *  @brief General solver for solving problems
 *
 *  @author Lukas Pospisil
 */

#ifndef PASCGENERALSOLVER_H
#define	PASCGENERALSOLVER_H

#include "pascinference.h"

namespace pascinference {
namespace solver {

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
		bool dump_or_not;  /**< dump the data of problem which is solved */

	public:
		/* settings */
		int debugmode; /**< print info about the progress */
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

		/** @brief print status of solver into stream
		 * 
		 *  Print the values of inner couters of the solver (for example iteration counter).
		 * 
		 * @param output stream where to print
		 */ 
		virtual void printstatus(std::ostringstream &output) const;

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

		virtual void printshort(std::ostringstream &header, std::ostringstream &values) const;
		virtual void printshort_sum(std::ostringstream &header, std::ostringstream &values) const;

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
		
		virtual int get_debugmode() const;
		virtual void set_debugmode(int debugmode);

		virtual int get_maxit() const;
		virtual void set_maxit(int maxit);

		virtual double get_eps() const;
		virtual void set_eps(double eps);

};

/* general print, call virtual print() */
ConsoleOutput &operator<<(ConsoleOutput &output, const GeneralSolver &solver){
	output << solver.get_name();
	return output;
}

void GeneralSolver::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN
	
	output << this->get_name() << std::endl;

	LOG_FUNC_END
}

void GeneralSolver::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN
	
	output_global << this->get_name() << std::endl;
	output_global.synchronize();

	LOG_FUNC_END
}

void GeneralSolver::printstatus(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN
	
	output << this->get_name() << ": status" << std::endl;

	LOG_FUNC_END
}

void GeneralSolver::printstatus(std::ostringstream &output) const {
	LOG_FUNC_BEGIN
	
	output << this->get_name() << ": status" << std::endl;

	LOG_FUNC_END
}

void GeneralSolver::printcontent(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN
	
	output << this->get_name() << ": content" << std::endl;

	LOG_FUNC_END
}

void GeneralSolver::printtimer(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN
	
	output << this->get_name() << ": timer" << std::endl;

	LOG_FUNC_END
}

void GeneralSolver::printshort(std::ostringstream &header, std::ostringstream &values) const{
	LOG_FUNC_BEGIN

	LOG_FUNC_END
}

void GeneralSolver::printshort_sum(std::ostringstream &header, std::ostringstream &values) const{
	LOG_FUNC_BEGIN

	LOG_FUNC_END
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

int GeneralSolver::get_debugmode() const {
	return this->debugmode;
}

void GeneralSolver::set_debugmode(int debugmode){
	this->debugmode = debugmode;
}

int GeneralSolver::get_maxit() const {
	return this->maxit;
}

void GeneralSolver::set_maxit(int maxit) {
	this->maxit = maxit;
}

double GeneralSolver::get_eps() const {
	return this->eps;
}

void GeneralSolver::set_eps(double eps) {
	this->eps = eps;
}



}
} /* end of namespace */


#endif
