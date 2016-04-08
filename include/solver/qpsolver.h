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
		SolverType child_solvertype;
	
		QPSolverSetting() {
			this->maxit = QPSOLVER_DEFAULT_MAXIT;
			this->eps = QPSOLVER_DEFAULT_EPS;
			this->debug_mode = QPSOLVER_DEFAULT_DEBUG_MODE;
			
			this->child_solvertype = SOLVER_AUTO;
		};
		~QPSolverSetting() {};

		virtual void print(std::ostream &output) const {
			output <<  this->get_name() << std::endl;
			output <<  " - maxit: " << this->maxit << std::endl;
			output <<  " - eps:   " << this->eps << std::endl;

		};

		std::string get_name() const {
			return "General QP SolverSetting";
		};
	

};


/* QPSolver */ 
template<class VectorBase>
class QPSolver: public GeneralSolver {
	protected:
		QPData<VectorBase> *qpdata; /* data on which the solver operates */

		QPSolver *child_solver; /* child of this solver, which actually solve the problem */
		
		double fx; /**< function value in actual iteration */
		int it_sum; /**< number of all iterations */
		int it_last; /**< number of iterations during last solution */
		int hessmult_sum; /**< number of all Hessian multiplication */
		int hessmult_last; /**< number of Hessian multiplication */
		
	public:
		QPSolverSetting setting;

		QPSolver();
		QPSolver(QPData<VectorBase> &new_qpdata); 
		~QPSolver();

		virtual void solve();

		virtual void print(std::ostream &output) const;
		virtual void printcontent(std::ostream &output) const;
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
	if(setting.debug_mode >= 100) coutMaster << "(QPSolver)CONSTRUCTOR" << std::endl;
	
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
	if(setting.debug_mode >= 100) coutMaster << "(QPSolver)DESTRUCTOR" << std::endl;

	/* destroy child solver */
	free(child_solver);
}


/* print info about problem */
template<class VectorBase>
void QPSolver<VectorBase>::print(std::ostream &output) const {
	if(setting.debug_mode >= 100) coutMaster << "(QPSolver)FUNCTION: print" << std::endl;

	output <<  this->get_name() << std::endl;
	
	/* print settings */
	coutMaster.push();
	setting.print(output);
	coutMaster.pop();

	/* print data */
	if(qpdata){
		coutMaster.push();
		qpdata->print(output);
		coutMaster.pop();
	}
	
	/* if child solver is specified, then print also info about it */	
	if(child_solver){
		child_solver->print(output);
	}	
}

/* print content of solver */
template<class VectorBase>
void QPSolver<VectorBase>::printcontent(std::ostream &output) const {
	if(setting.debug_mode >= 100) coutMaster << "(CGQPSolver)FUNCTION: printcontent" << std::endl;

	output << this->get_name() << std::endl;
	
	/* print content of data */
	if(qpdata){
		output << "- data:" << std::endl;
		coutMaster.push();
		qpdata->printcontent(output);
		coutMaster.pop();
	}
		
}

/* print status */
template<class VectorBase>
void QPSolver<VectorBase>::printstatus(std::ostream &output) const {
	if(child_solver){
		child_solver->printstatus(output);
	} else {
		output <<  this->get_name() << std::endl;
	}
}

/* print timer */
template<class VectorBase>
void QPSolver<VectorBase>::printtimer(std::ostream &output) const {
	if(child_solver){
		child_solver->printtimer(output);
	} else {
		output <<  this->get_name() << std::endl;
	}
}

template<class VectorBase>
std::string QPSolver<VectorBase>::get_name() const {
	return "General QP Solver";
}

/* solve the problem */
template<class VectorBase>
void QPSolver<VectorBase>::solve() {
	if(setting.debug_mode >= 100) coutMaster << "(QPSolver)FUNCTION: solve" << std::endl;

	/* the child solver wasn't specified yet */
	if(!child_solver){
		/* which specific solver we can use to solve the problem? */
		if(setting.child_solvertype == SOLVER_AUTO){
			//TODO: here write more sophisticated auto decision tree
			if(qpdata->get_feasibleset()){
				/* with constraints */
				setting.child_solvertype = SOLVER_SPGQP;
			} else {
				/* unconstrained */
				setting.child_solvertype = SOLVER_CG;
			}
		} 

		/* prepare CG solver */
		if(setting.child_solvertype == SOLVER_CG){
			/* create new instance of CG Solver */
			child_solver = new CGQPSolver<VectorBase>(*qpdata);
		}

		/* prepare SPGQP solver */
		if(setting.child_solvertype == SOLVER_SPGQP){
			/* create new instance of CG Solver */
			child_solver = new SPGQPSolver<VectorBase>(*qpdata);
		}
	}

	/* update settings of child solver */
	child_solver->setting.debug_mode = setting.debug_mode;

	/* now the child_solver should be specified and prepared */

	/* call solve function to child solver */
	child_solver->solve();
	
}

template<class VectorBase>
double QPSolver<VectorBase>::get_fx() const {
	if(setting.debug_mode >= 100) coutMaster << "(QPSolver)FUNCTION: get_fx()" << std::endl;

	return child_solver->get_fx(); // TODO: control existence
}

template<class VectorBase>
int QPSolver<VectorBase>::get_it() const {
	if(setting.debug_mode >= 100) coutMaster << "(QPSolver)FUNCTION: get_function_value()" << std::endl;

	return child_solver->get_it(); // TODO: control existence
}

template<class VectorBase>
QPData<VectorBase> *QPSolver<VectorBase>::get_data() const {
	return qpdata;
}

} /* end namespace */

#endif
