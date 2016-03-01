#ifndef QPSOLVER_H
#define	QPSOLVER_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "generalsolver.h"
#include "algebra.h"

#include "data/qpdata.h"
#include "result/qpresult.h"

#define QPSOLVER_DEFAULT_MAXIT 1000;
#define QPSOLVER_DEFAULT_EPS 0.0001;

namespace pascinference {

/* settings */
class QPSolverSetting : public GeneralSolverSetting {
	protected:
		int maxit;
		double eps;
		
	public:
		QPSolverSetting() {
			maxit = QPSOLVER_DEFAULT_MAXIT;
			eps = QPSOLVER_DEFAULT_EPS;
			
		};
		~QPSolverSetting() {};

		virtual void print(std::ostream &output) const {
			output << "  QPSolverSettings:" << std::endl;
			output << "   - maxit: " << maxit << std::endl;
			output << "   - eps: " << eps << std::endl;

		};
		
};


/* QPSolver */ 
template<class VectorBase>
class QPSolver: public GeneralSolver {
	protected:
		const QPData<VectorBase> *data; /* data on which the solver operates */
		const QPResult<VectorBase> *result; /* here solver stores results */
		
	public:
		QPSolverSetting setting;

		QPSolver();
		QPSolver(const QPData<VectorBase> &new_data, const QPResult<VectorBase> &new_result); 
		~QPSolver();

		virtual void solve(SolverType solvertype);
		virtual void solve() { this->solve(SOLVER_AUTO); };

		virtual void print(std::ostream &output) const;


};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

#include "solver/cgqpsolver.h"

namespace pascinference {

/* constructor */
template<class VectorBase>
QPSolver<VectorBase>::QPSolver(){
	if(DEBUG_MODE >= 100) std::cout << "(QPSolver)CONSTRUCTOR" << std::endl;
	
	data = NULL;
	result = NULL;
}

template<class VectorBase>
QPSolver<VectorBase>::QPSolver(const QPData<VectorBase> &new_data, const QPResult<VectorBase> &new_result){
	data = &new_data;
	result = &new_result;
}

/* destructor */
template<class VectorBase>
QPSolver<VectorBase>::~QPSolver(){
	if(DEBUG_MODE >= 100) std::cout << "(QPSolver)DESTRUCTOR" << std::endl;
	
}


/* print info about problem */
template<class VectorBase>
void QPSolver<VectorBase>::print(std::ostream &output) const {
	if(DEBUG_MODE >= 100) std::cout << "(QPSolver)FUNCTION: print" << std::endl;

	output << " QPSolver" << std::endl;
	
	/* print settings */
	output << setting;
		
}

/* solve the problem */
template<class VectorBase>
void QPSolver<VectorBase>::solve(SolverType solvertype) {
	if(DEBUG_MODE >= 100) std::cout << "(QPSolver)FUNCTION: solve" << std::endl;

	GeneralSolver *solver = NULL;

	/* which specific solver we can use to solve the problem? */
	if(solvertype == SOLVER_AUTO){
		//TODO: here write sophisticated auto decision tree
	} 

	/* prepare CG solver */
	if(solvertype == SOLVER_CG){
		/* create new instance of CG Solver */
		solver = new CGQPSolver<VectorBase>(*data,*result);
		
		/* copy settings */
		
	}
	
	if(!solver){
		/* solver was not initialized - there is not a suitable solver */
		//TODO: give error
		std::cout << "something goes wrong with solvers !!!" << std::endl;
	} else {
		/* solve the problem using specific solver */
		solver->solve();

		/* destroy solver */
	}
	
}


} /* end namespace */

#endif
