#ifndef PASC_QPSOLVER_GLOBAL_H
#define	PASC_QPSOLVER_GLOBAL_H

/* for debugging, if >= 100, then print info about ach called function */
extern int DEBUG_MODE;

#include <iostream>
#include "generalsolver.h"
#include "algebra.h"

#include "solver/qpsolver.h"
#include "data/qpdata.h"

#define QPSOLVER_GLOBAL_DEFAULT_MAXIT 1000;
#define QPSOLVER_GLOBAL_DEFAULT_EPS 0.0001;
#define QPSOLVER_GLOBAL_DEFAULT_DEBUG_MODE 0;

namespace pascinference {

/* settings */
class QPSolver_GlobalSetting : public GeneralSolverSetting {
	protected:

	public:
	
		QPSolver_GlobalSetting() {
			this->maxit = QPSOLVER_GLOBAL_DEFAULT_MAXIT;
			this->eps = QPSOLVER_GLOBAL_DEFAULT_EPS;
			this->debug_mode = QPSOLVER_GLOBAL_DEFAULT_DEBUG_MODE;
		};
		~QPSolver_GlobalSetting() {};

		virtual void print(ConsoleOutput &output) const {
			output <<  this->get_name() << std::endl;
			output <<  " - maxit: " << this->maxit << std::endl;
			output <<  " - eps:   " << this->eps << std::endl;

		};

		std::string get_name() const {
			return "General Global QP SolverSetting";
		};
	

};


/* QPSolver_Global */ 
class QPSolver_Global: public QPSolver<PetscVector> {
	protected:
		const QPData<PetscVector> *qpdata; /* data on which the solver operates */
		QPData<PetscVector> *data_local; /* data of inner qp solver */		
		QPSolver<PetscVector> *solver_local; /* local solver */

		/* local data */
		GeneralVector<PetscVector> *x;
		GeneralVector<PetscVector> *x0;
		GeneralVector<PetscVector> *b;

		void GetLocalData();
		void RestoreLocalData();

	public:
		QPSolver_GlobalSetting setting;

		QPSolver_Global(const QPData<PetscVector> &new_qpdata); 
		~QPSolver_Global();

		virtual void solve();
		double get_fx() const;
		int get_it() const;
		int get_hessmult() const;

		virtual void print(ConsoleOutput &output) const;
		virtual void print(ConsoleOutput &output_global, ConsoleOutput &output_local) const;
		
		virtual void printstatus(ConsoleOutput &output) const;
		virtual void printtimer(ConsoleOutput &output) const;

		virtual void printcontent(ConsoleOutput &output) const;
		virtual std::string get_name() const;

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
QPSolver_Global::QPSolver_Global(const QPData<PetscVector> &new_qpdata){
	LOG_FUNC_BEGIN

	qpdata = &new_qpdata;

	this->fx = -1;
	this->it_last = 0; 
	this->hessmult_last = 0; 

	/* initialize local QP solver */
	this->data_local = new QPData<PetscVector>();

	int local_size;
	Vec x_local;
	Vec x0_local;
	Vec b_local;
	
	/* get local size */
	TRY( VecGetLocalSize(qpdata->get_x()->get_vector(), &local_size) );

	/* allocate local data */
	TRY( VecCreateSeq(PETSC_COMM_SELF, local_size, &x_local) );	
	TRY( VecCreateSeq(PETSC_COMM_SELF, local_size, &x0_local) );	
	TRY( VecCreateSeq(PETSC_COMM_SELF, local_size, &b_local) );	

	x = new GeneralVector<PetscVector>(x_local);
	x0 = new GeneralVector<PetscVector>(x0_local);
	b = new GeneralVector<PetscVector>(b_local);

	/* get local vectors and prepare local data */
	data_local->set_A(qpdata->get_A());
	data_local->set_b(b);
	data_local->set_x(x);
	data_local->set_x0(x0);
	data_local->set_feasibleset(qpdata->get_feasibleset());
	
	/* create new instance of local solver */
	solver_local = new QPSolver<PetscVector>(*data_local);

	LOG_FUNC_END
}

/* destructor */
QPSolver_Global::~QPSolver_Global(){
	LOG_FUNC_BEGIN

	free(this->data_local);
	free(this->solver_local);

	LOG_FUNC_END
}

void QPSolver_Global::GetLocalData(){
	LOG_FUNC_BEGIN

	TRY( VecGetLocalVector(qpdata->get_x()->get_vector(),data_local->get_x()->get_vector()) );
	TRY( VecGetLocalVector(qpdata->get_x0()->get_vector(),data_local->get_x0()->get_vector()) );
	TRY( VecGetLocalVector(qpdata->get_b()->get_vector(),data_local->get_b()->get_vector()) );

	LOG_FUNC_END
}

void QPSolver_Global::RestoreLocalData(){
	LOG_FUNC_BEGIN

	TRY( VecRestoreLocalVector(qpdata->get_x()->get_vector(),data_local->get_x()->get_vector()) );
	TRY( VecRestoreLocalVector(qpdata->get_x0()->get_vector(),data_local->get_x0()->get_vector()) );
	TRY( VecRestoreLocalVector(qpdata->get_b()->get_vector(),data_local->get_b()->get_vector()) );

	LOG_FUNC_END
}

/* print info about problem */
void QPSolver_Global::print(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

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

	LOG_FUNC_END
}

/* print info about problem */
void QPSolver_Global::print(ConsoleOutput &output_global, ConsoleOutput &output_local) const {
	LOG_FUNC_BEGIN

	output_global <<  this->get_name() << std::endl;
	
	/* print settings */
	output_global.push();
	setting.print(output_global);
	output_global.pop();

	output_global << " - local solver:" << std::endl;
	output_global.push();
	if(solver_local){
		solver_local->printstatus(output_local);
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
void QPSolver_Global::printcontent(ConsoleOutput &output) const {
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
void QPSolver_Global::printstatus(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	if(solver_local){
		solver_local->printstatus(output);
	} else {
		output <<  " - status: not set yet." << std::endl; 
	}

	LOG_FUNC_END
}


/* print timer */
void QPSolver_Global::printtimer(ConsoleOutput &output) const {
	LOG_FUNC_BEGIN

	if(solver_local){
		solver_local->printtimer(output);
	} else {
		output <<  " - timer: not set yet." << std::endl; 
	}

	LOG_FUNC_END
}

std::string QPSolver_Global::get_name() const {
	return "GeneralQPSolver_Global";
}

/* solve the problem */
void QPSolver_Global::solve() {
	LOG_FUNC_BEGIN

	/* for each process prepare QP solver and solve the problem */

	GetLocalData();
		
	/* solve local problem */
	solver_local->solve();

	this->fx = solver_local->get_fx();
	this->it_last = solver_local->get_it();
	this->hessmult_last = solver_local->get_hessmult();

	RestoreLocalData();

	LOG_IT(this->it_last)
	LOG_FX(this->fx)
	
	LOG_FUNC_END
}

double QPSolver_Global::get_fx() const {
	if(setting.debug_mode >= 100) coutMaster << "(QPSolver_Global)FUNCTION: get_fx()" << std::endl;

	return this->fx; 
}

int QPSolver_Global::get_it() const {
	if(setting.debug_mode >= 100) coutMaster << "(QPSolver_Global)FUNCTION: get_it()" << std::endl;

	return this->it_last; 
}

int QPSolver_Global::get_hessmult() const {
	if(setting.debug_mode >= 100) coutMaster << "(QPSolver_Global)FUNCTION: get_hessmult()" << std::endl;

	return this->hessmult_last; 
}


} /* end namespace */

#endif
