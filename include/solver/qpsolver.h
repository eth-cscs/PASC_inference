#ifndef PASC_QPSOLVER_H
#define	PASC_QPSOLVER_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "generalsolver.h"
#include "algebra.h"

#include "data/qpdata.h"

#define QPSOLVER_DEFAULT_MAXIT 1000;
#define QPSOLVER_DEFAULT_EPS 0.0001;
#define QPSOLVER_DEFAULT_DEBUG_MODE 0;

namespace pascinference {

/* settings */
class QPSolverSetting : public GeneralSolverSetting {
	protected:

	public:
		QPSolverSetting() {
			this->maxit = QPSOLVER_DEFAULT_MAXIT;
			this->eps = QPSOLVER_DEFAULT_EPS;
			this->debug_mode = QPSOLVER_DEFAULT_DEBUG_MODE;
		};
		~QPSolverSetting() {};

		virtual void print(std::ostream &output) const {
			output << "  QPSolverSettings:" << std::endl;
			output << "   - maxit: " << this->maxit << std::endl;
			output << "   - eps: " << this->eps << std::endl;

		};
	

};


/* QPSolver */ 
template<class VectorBase>
class QPSolver: public GeneralSolver {
	protected:
		QPData<VectorBase> *qpdata; /* data on which the solver operates */

		QPSolver *child_solver; /* child of this solver, which actually solve the problem */
		
		double fx; /**< function value in actual iteration */
		int it; /**< actual iteration */
		int hessmult; /**< number of Hessian multiplication */
		
	public:
		QPSolverSetting setting;

		QPSolver();
		QPSolver(QPData<VectorBase> &new_qpdata); 
		~QPSolver();

		virtual void solve(SolverType solvertype);
		virtual void solve() { this->solve(SOLVER_AUTO); };

		virtual void print(std::ostream &output) const;
		virtual void printstatus(std::ostream &output) const;
		virtual void printtimer(std::ostream &output) const;
		virtual std::string get_name() const;

		virtual QPData<VectorBase> *get_data() const;
		virtual double get_fx() const;
		virtual int get_it() const;

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

#include "solver/cgqpsolver.h"
#include "solver/spgqpsolver.h"

namespace pascinference {

/* constructor */
template<class VectorBase>
QPSolver<VectorBase>::QPSolver(){
	if(DEBUG_MODE >= 100) std::cout << "(QPSolver)CONSTRUCTOR" << std::endl;
	
	qpdata = NULL;
	child_solver = NULL; /* in this time, we don't know how to solve the problem */
	
	fx = std::numeric_limits<double>::max();
}

template<class VectorBase>
QPSolver<VectorBase>::QPSolver(QPData<VectorBase> &new_qpdata){
	qpdata = &new_qpdata;

	child_solver = NULL; /* in this time, we don't know how to solve the problem */
	
	fx = std::numeric_limits<double>::max();
}

/* destructor */
template<class VectorBase>
QPSolver<VectorBase>::~QPSolver(){
	if(DEBUG_MODE >= 100) std::cout << "(QPSolver)DESTRUCTOR" << std::endl;

	/* destroy child solver */
	free(child_solver);
}


/* print info about problem */
template<class VectorBase>
void QPSolver<VectorBase>::print(std::ostream &output) const {
	if(DEBUG_MODE >= 100) std::cout << "(QPSolver)FUNCTION: print" << std::endl;

	output << this->get_name() << std::endl;
	
	/* print settings */
	output << setting;

	/* print data */
	output << *qpdata;
	
	/* if child solver is specified, then print also info about it */	
	if(child_solver){
		child_solver->print(output);
	}	
}

/* print status */
template<class VectorBase>
void QPSolver<VectorBase>::printstatus(std::ostream &output) const {
	child_solver->printstatus(output); // TODO: check existence	
}

/* print timer */
template<class VectorBase>
void QPSolver<VectorBase>::printtimer(std::ostream &output) const {
	child_solver->printtimer(output); // TODO: check existence	
}

template<class VectorBase>
std::string QPSolver<VectorBase>::get_name() const {
	return "General QP Solver";
}

/* solve the problem */
template<class VectorBase>
void QPSolver<VectorBase>::solve(SolverType solvertype) {
	if(DEBUG_MODE >= 100) std::cout << "(QPSolver)FUNCTION: solve" << std::endl;

	/* the child solver wasn't specified yet */
	if(!child_solver){
		/* which specific solver we can use to solve the problem? */
		if(solvertype == SOLVER_AUTO){
			//TODO: here write sophisticated auto decision tree
		} 

		/* prepare CG solver */
		if(solvertype == SOLVER_CG){
			/* create new instance of CG Solver */
			child_solver = new CGQPSolver<VectorBase>(*qpdata);
		
			/* copy settings */
//			child_solver->setting.maxit = setting.maxit;
//			child_solver->setting.eps = setting.eps;
//			child_solver->setting.debug_mode = DEBUG_MODE;
		}

		/* prepare SPGQP solver */
		if(solvertype == SOLVER_SPGQP){
			/* create new instance of CG Solver */
			child_solver = new SPGQPSolver<VectorBase>(*qpdata);
		
			/* copy settings */
//			child_solver->setting.maxit = setting.maxit;
//			child_solver->setting.eps = setting.eps;
		
		}
	}
	/* now the child_solver should be specified and prepared */

	/* call solve function to child solver */
	child_solver->solve();
	
}

template<class VectorBase>
double QPSolver<VectorBase>::get_fx() const {
	if(DEBUG_MODE >= 100) std::cout << "(QPSolver)FUNCTION: get_fx()" << std::endl;

	return child_solver->get_fx(); // TODO: control existence
}

template<class VectorBase>
int QPSolver<VectorBase>::get_it() const {
	if(DEBUG_MODE >= 100) std::cout << "(QPSolver)FUNCTION: get_function_value()" << std::endl;

	return child_solver->get_it(); // TODO: control existence
}

template<class VectorBase>
QPData<VectorBase> *QPSolver<VectorBase>::get_data() const {
	return qpdata;
}

} /* end namespace */

#endif