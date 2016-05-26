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

		SolverType child_solvertype;		
	public:
		QPSolver();
		QPSolver(QPData<VectorBase> &new_qpdata); 
		~QPSolver();

		virtual void solve();

		virtual void print(ConsoleOutput &output) const;
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		virtual void printcontent(ConsoleOutput &output) const;
		virtual void printstatus(ConsoleOutput &output) const;
		virtual void printtimer(ConsoleOutput &output) const;
		virtual std::string get_name() const;

		virtual QPData<VectorBase> *get_data() const;
		virtual double get_fx() const;
		virtual int get_it() const;
		virtual int get_hessmult() const;

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
	LOG_FUNC_BEGIN
	
	qpdata = NULL;
	child_solver = NULL; /* in this time, we don't know how to solve the problem */
	
	fx = std::numeric_limits<double>::max();

	this->maxit = QPSOLVER_DEFAULT_MAXIT;
	this->eps = QPSOLVER_DEFAULT_EPS;
	this->debug_mode = QPSOLVER_DEFAULT_DEBUG_MODE;
			
	this->child_solvertype = SOLVER_AUTO;

	LOG_FUNC_END
}

template<class VectorBase>
QPSolver<VectorBase>::QPSolver(QPData<VectorBase> &new_qpdata){
	LOG_FUNC_BEGIN

	qpdata = &new_qpdata;
	child_solver = NULL; /* in this time, we don't know how to solve the problem */
	fx = std::numeric_limits<double>::max();

	this->maxit = QPSOLVER_DEFAULT_MAXIT;
	this->eps = QPSOLVER_DEFAULT_EPS;
	this->debug_mode = QPSOLVER_DEFAULT_DEBUG_MODE;
			
	this->child_solvertype = SOLVER_AUTO;

	LOG_FUNC_END
}

/* destructor */
template<class VectorBase>
QPSolver<VectorBase>::~QPSolver(){
	LOG_FUNC_BEGIN

	/* destroy child solver */
	free(child_solver);

	LOG_FUNC_END
}


/* print info about problem */
template<class VectorBase>
void QPSolver<VectorBase>::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output <<  this->get_name() << std::endl;
	
	/* print settings */
	output <<  " - maxit:      " << this->maxit << std::endl;
	output <<  " - eps:        " << this->eps << std::endl;
	output <<  " - debug_mode: " << this->debug_mode << std::endl;

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

	LOG_FUNC_END
}

template<class VectorBase>
void QPSolver<VectorBase>::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* print settings */
	output_global <<  " - maxit:      " << this->maxit << std::endl;
	output_global <<  " - eps:        " << this->eps << std::endl;
	output_global <<  " - debug_mode: " << this->debug_mode << std::endl;

	/* print data */
	if(qpdata){
		output_global.push();
		qpdata->print(output_global);
		output_global.pop();
	}
	
	/* if child solver is specified, then print also info about it */	
	if(child_solver){
		output_global.push();
		child_solver->print(output_global, output_local);
		output_global.pop();
	}	

	LOG_FUNC_END
}

/* print content of solver */
template<class VectorBase>
void QPSolver<VectorBase>::printcontent(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	output << this->get_name() << std::endl;
	
	/* print content of data */
	if(qpdata){
		output << "- data:" << std::endl;
		coutMaster.push();
		qpdata->printcontent(output);
		coutMaster.pop();
	}
		
	LOG_FUNC_END
}

/* print status */
template<class VectorBase>
void QPSolver<VectorBase>::printstatus(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	if(child_solver){
		child_solver->printstatus(output);
	} else {
		output << this->get_name() << ": status" << std::endl;
	}

	LOG_FUNC_END
}

/* print timer */
template<class VectorBase>
void QPSolver<VectorBase>::printtimer(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	if(child_solver){
		child_solver->printtimer(output);
	} else {
		output << this->get_name() << ": timer" << std::endl;
	}

	LOG_FUNC_END
}

template<class VectorBase>
std::string QPSolver<VectorBase>::get_name() const {
	return "GeneralQPSolver";
}

/* solve the problem */
template<class VectorBase>
void QPSolver<VectorBase>::solve() {
	LOG_FUNC_BEGIN

	/* the child solver wasn't specified yet */
	if(!child_solver){
		/* which specific solver we can use to solve the problem? */
		if(this->child_solvertype == SOLVER_AUTO){
			//TODO: here write more sophisticated auto decision tree
			if(qpdata->get_feasibleset()){
				/* with constraints */
				this->child_solvertype = SOLVER_SPGQP;
			} else {
				/* unconstrained */
				this->child_solvertype = SOLVER_CG;
			}
		} 

		/* prepare CG solver */
		if(this->child_solvertype == SOLVER_CG){
			/* create new instance of CG Solver */
			child_solver = new CGQPSolver<VectorBase>(*qpdata);
		}

		/* prepare SPGQP solver */
		if(this->child_solvertype == SOLVER_SPGQP){
			/* create new instance of CG Solver */
			child_solver = new SPGQPSolver<VectorBase>(*qpdata);
		}
	}

	/* update settings of child solver */
	child_solver->debug_mode = this->debug_mode;
	child_solver->eps = this->eps;

	/* call solve function to child solver */
	child_solver->solve();

	this->fx = child_solver->get_fx();
	this->it_last = child_solver->get_it();
	this->hessmult_last = child_solver->get_hessmult();
	
	LOG_IT(this->it_last)
	LOG_FX(this->fx)
		
	LOG_FUNC_END
}

template<class VectorBase>
double QPSolver<VectorBase>::get_fx() const {
	return child_solver->get_fx(); // TODO: control existence
}

template<class VectorBase>
int QPSolver<VectorBase>::get_it() const {
	return child_solver->get_it(); // TODO: control existence
}

template<class VectorBase>
int QPSolver<VectorBase>::get_hessmult() const {
	return child_solver->get_hessmult(); // TODO: control existence
}

template<class VectorBase>
QPData<VectorBase> *QPSolver<VectorBase>::get_data() const {
	return qpdata;
}

} /* end namespace */

#endif
