#ifndef TSSOLVER_H
#define	TSSOLVER_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "generalsolver.h"
#include "algebra.h"

#include "data/tsdata.h"
#include "result/tsresult.h"

#define TSSOLVER_DEFAULT_MAXIT 1000;
#define TSSOLVER_DEFAULT_EPS 0.0001;

namespace pascinference {

/* settings */
class TSSolverSetting : public GeneralSolverSetting {
	protected:
		int maxit;
		double eps;

	public:
		TSSolverSetting() {
			maxit = TSSOLVER_DEFAULT_MAXIT;
			eps = TSSOLVER_DEFAULT_EPS;
			
		};
		~TSSolverSetting() {};

		virtual void print(std::ostream &output) const {
			output << " TSSolverSettings:" << std::endl;
			output << "   - maxit: " << maxit << std::endl;
			output << "   - eps: " << eps << std::endl;

		};
		
};


/* TSSolver */ 
template<class VectorBase>
class TSSolver: public GeneralSolver {
	protected:
		const TSData<VectorBase> *data; /* data on which the solver operates */
		const TSResult<VectorBase> *result; /* here solver stores results */

		GeneralSolver *gamma_solver; /* to solve inner gamma problem */
		GeneralData 	*gamma_data; /* data for gamma problem */
		GeneralResult *gamma_result; /* result for gamma problem */
		
		GeneralSolver *theta_solver; /* to solve inner theta problem */
		GeneralData 	*theta_data; /* data for theta problem */
		GeneralResult *theta_result; /* result for theta problem */


	public:
		TSSolverSetting setting;

		TSSolver();
		TSSolver(const TSData<VectorBase> &new_data, const TSResult<VectorBase> &new_result); 
		~TSSolver();

		/* with given gamma and theta solvertype */
		void solve(SolverType gamma_solvertype, SolverType theta_solvertype);
		
		/* with one solvertype */
		virtual void solve(SolverType solvertype) {
			this->solve(solvertype,solvertype); 
		};
		
		/* without given solvertype */
		virtual void solve() {
			this->solve(SOLVER_AUTO); 
		};

		virtual void print(std::ostream &output) const;
		std::string get_name() const;


};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
template<class VectorBase>
TSSolver<VectorBase>::TSSolver(){
	if(DEBUG_MODE >= 100) std::cout << "(TSSolver)CONSTRUCTOR" << std::endl;
	
	data = NULL;
	result = NULL;
	gamma_solver = NULL; /* in this time, we don't know how to solve the problem */
	theta_solver = NULL; /* in this time, we don't know how to solve the problem */

}

template<class VectorBase>
TSSolver<VectorBase>::TSSolver(const TSData<VectorBase> &new_data, const TSResult<VectorBase> &new_result){
	data = &new_data;
	result = &new_result;

	gamma_solver = NULL; /* in this time, we don't know how to solve the problem */
	theta_solver = NULL; /* in this time, we don't know how to solve the problem */


}

/* destructor */
template<class VectorBase>
TSSolver<VectorBase>::~TSSolver(){
	if(DEBUG_MODE >= 100) std::cout << "(TSSolver)DESTRUCTOR" << std::endl;

	/* destroy child solvers */
	free(gamma_solver);
	free(theta_solver);

}


/* print info about problem */
template<class VectorBase>
void TSSolver<VectorBase>::print(std::ostream &output) const {
	if(DEBUG_MODE >= 100) std::cout << "(TSSolver)FUNCTION: print" << std::endl;

	output << this->get_name() << std::endl;
	
	/* print settings */
	output << setting;
	
	/* if child solvers are specified, then print also info about it */	
	output << " Gamma Solver" << std::endl;
	if(gamma_solver){
		output << gamma_solver << std::endl;
	} else {
		output << "  - not set" << std::endl;
	}	

	output << " Theta Solver" << std::endl;
	if(theta_solver){
		output << theta_solver << std::endl;
	} else {
		output << "  - not set" << std::endl;
	}	


}

template<class VectorBase>
std::string TSSolver<VectorBase>::get_name() const {
	return "Time-Series Solver";
}

/* solve the problem */
template<class VectorBase>
void TSSolver<VectorBase>::solve(SolverType gamma_solvertype, SolverType theta_solvertype) {
	if(DEBUG_MODE >= 100) std::cout << "(TSSolver)FUNCTION: solve" << std::endl;

	/* the gamma solver wasn't specified yet */
	if(!gamma_solvertype){
		/* which specific solver we can use to solve the problem? */
		if(gamma_solvertype == SOLVER_AUTO){
			// TODO: here write sofisticated decision tree - maybe based on MODEL?
		} 

		/* prepare QP solver */
		if(gamma_solvertype == SOLVER_QP){
			/* create new instance of QP Solver to solve gamma problem */
			
			
//			gamma_solver = new QPSolver<VectorBase>(*data,*result);
		
			/* copy settings */
			//TODO
		}

	}
	/* now the child_solver should be specified and prepared */

	/* call solve function to child solver */
//	child_solver->solve();
	
}


} /* end namespace */

#endif
