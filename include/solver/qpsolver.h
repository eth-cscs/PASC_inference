/** @file qpsolver.h
 *  @brief General Quadratic Programming solver
 *
 *  @author Lukas Pospisil
 */
 
#ifndef PASC_QPSOLVER_H
#define	PASC_QPSOLVER_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUGMODE;

#include "pascinference.h"
#include "data/qpdata.h"

#define QPSOLVER_DEFAULT_MAXIT 1000
#define QPSOLVER_DEFAULT_EPS 1e-9
#define QPSOLVER_DEFAULT_DEBUGMODE 0

namespace pascinference {
namespace solver {

/** \class QPSolver
 *  \brief General Quadratic Programming solver
 *
*/
template<class VectorBase>
class QPSolver: public GeneralSolver {
	protected:
		QPData<VectorBase> *qpdata; /**< data on which the solver operates */

		QPSolver *child_solver; /**< child of this solver, which actually solves the problem */
		
		double fx; /**< function value in actual iteration */
		int it_sum; /**< number of all iterations */
		int it_last; /**< number of iterations during last solution */
		int hessmult_sum; /**< number of all Hessian multiplication */
		int hessmult_last; /**< number of Hessian multiplication */

		SolverType child_solvertype;	
	public:
		/** @brief default constructor
		 * 
		 * Set inner data to null.
		 * 
		 */ 
		QPSolver();
		
		/** @brief constructor from data
		 * 
		 * Set inner pointer to data on which the solver operates to given data.
		 * 
		 * @param new_qpdata input data for solver
		 */ 		
		QPSolver(QPData<VectorBase> &new_qpdata); 

		/** @brief destructor
		 * 
		 */ 
		~QPSolver();

		virtual void solve();

		virtual void print(ConsoleOutput &output) const;
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		virtual void printcontent(ConsoleOutput &output) const;
		virtual void printstatus(ConsoleOutput &output) const;
		virtual void printstatus(std::ostringstream &output) const;
		virtual void printtimer(ConsoleOutput &output) const;
		virtual void printshort(std::ostringstream &header, std::ostringstream &values) const;
		virtual void printshort_sum(std::ostringstream &header, std::ostringstream &values) const;
		virtual std::string get_name() const;

		virtual QPData<VectorBase> *get_data() const;
		
		/** @brief return function value in last computed iteration
		* 
		*/		
		virtual double get_fx() const;

		/** @brief return number of iterations during last solve()
		* 
		*/		
		virtual int get_it() const;

		/** @brief return number of Hessian multiplication during last solve()
		* 
		*/		
		virtual int get_hessmult() const;

};


}
} /* end of namespace */

/* ------------- implementation ----------- */
//TODO: move to impls

#include "solver/cgqpsolver.h"
#include "solver/spgqpsolver.h"

namespace pascinference {
namespace solver {

/* constructor */
template<class VectorBase>
QPSolver<VectorBase>::QPSolver(){
	LOG_FUNC_BEGIN
	
	qpdata = NULL;
	child_solver = NULL; /* in this time, we don't know how to solve the problem */
	
	fx = std::numeric_limits<double>::max();

	consoleArg.set_option_value("qpsolver_maxit", &this->maxit, QPSOLVER_DEFAULT_MAXIT);
	consoleArg.set_option_value("qpsolver_eps", &this->eps, QPSOLVER_DEFAULT_EPS);
	consoleArg.set_option_value("qpsolver_debugmode", &this->debugmode, QPSOLVER_DEFAULT_DEBUGMODE);

	this->child_solvertype = SOLVER_AUTO;

	LOG_FUNC_END
}

template<class VectorBase>
QPSolver<VectorBase>::QPSolver(QPData<VectorBase> &new_qpdata){
	LOG_FUNC_BEGIN

	qpdata = &new_qpdata;
	child_solver = NULL; /* in this time, we don't know how to solve the problem */
	fx = std::numeric_limits<double>::max();

	consoleArg.set_option_value("qpsolver_maxit", &this->maxit, QPSOLVER_DEFAULT_MAXIT);
	consoleArg.set_option_value("qpsolver_eps", &this->eps, QPSOLVER_DEFAULT_EPS);
	consoleArg.set_option_value("qpsolver_debugmode", &this->debugmode, QPSOLVER_DEFAULT_DEBUGMODE);
			
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
	output <<  " - debugmode: " << this->debugmode << std::endl;

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
	output_global <<  " - maxit:       " << this->maxit << std::endl;
	output_global <<  " - eps:         " << this->eps << std::endl;
	output_global <<  " - debugmode:  " << this->debugmode << std::endl;

	output_global <<  " - child solver: " << std::endl;
	output_global.push();
	if(child_solver){
		child_solver->printstatus(output_local);
	} else {
		output_local <<  " - not set yet." << std::endl; 
	}
	output_local.synchronize();	
	output_global.pop();
	
	/* print data */
	if(qpdata){
		output_global.push();
		qpdata->print(output_global, output_local);
		output_global.pop();
	}
	
	output_global.synchronize();

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

template<class VectorBase>
void QPSolver<VectorBase>::printstatus(std::ostringstream &output) const {
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
void QPSolver<VectorBase>::printshort(std::ostringstream &header, std::ostringstream &values) const {
	LOG_FUNC_BEGIN

	if(child_solver){
		child_solver->printshort(header,values);
	}

	LOG_FUNC_END
}

template<class VectorBase>
void QPSolver<VectorBase>::printshort_sum(std::ostringstream &header, std::ostringstream &values) const {
	LOG_FUNC_BEGIN

	if(child_solver){
		child_solver->printshort_sum(header,values);
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
	child_solver->set_debugmode(this->get_debugmode());
//	child_solver->set_maxit(this->get_maxit());
//	child_solver->set_eps(this->get_eps());

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
	double fx;
	if(child_solver){
		fx = child_solver->get_fx();
	} else {
		fx = this->fx;
	}
	return fx;
}

template<class VectorBase>
int QPSolver<VectorBase>::get_it() const {
	int it;
	if(child_solver){
		it = child_solver->get_it();
	} else {
		it = this->it_last;
	}
	return it;
}

template<class VectorBase>
int QPSolver<VectorBase>::get_hessmult() const {
	int hessmult;
	if(child_solver){
		hessmult = child_solver->get_hessmult();
	} else {
		hessmult = this->hessmult_last;
	}
	return hessmult;
}

template<class VectorBase>
QPData<VectorBase> *QPSolver<VectorBase>::get_data() const {
	return qpdata;
}


}
} /* end namespace */

#endif
