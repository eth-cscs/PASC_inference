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

		virtual void print(std::ostream &output) const {
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
		
	public:
		QPSolver_GlobalSetting setting;

		QPSolver_Global();
		QPSolver_Global(const QPData<PetscVector> &new_qpdata); 
		~QPSolver_Global();

		void solve();
		double get_fx() const;
		int get_it() const;
		int get_hessmult() const;

		void print(std::ostream &output) const;
		void printstatus(std::ostream &output) const;
		void printcontent(std::ostream &output) const;
		std::string get_name() const;

};

} // end of namespace

/* ------------- implementation ----------- */
//TODO: move to impls

namespace pascinference {

/* constructor */
QPSolver_Global::QPSolver_Global(){
	if(setting.debug_mode >= 100) coutMaster << "(QPSolver_Global)CONSTRUCTOR" << std::endl;
	
	qpdata = NULL;
	
	this->fx = -1;
	this->it_last = 0; 
	this->hessmult_last = 0;
	
}

QPSolver_Global::QPSolver_Global(const QPData<PetscVector> &new_qpdata){
	qpdata = &new_qpdata;

	this->fx = -1;
	this->it_last = 0; 
	this->hessmult_last = 0; 
}

/* destructor */
QPSolver_Global::~QPSolver_Global(){
	if(setting.debug_mode >= 100) coutMaster << "(QPSolver_Global)DESTRUCTOR" << std::endl;

}


/* print info about problem */
void QPSolver_Global::print(std::ostream &output) const {
	if(setting.debug_mode >= 100) coutMaster << "(QPSolver_Global)FUNCTION: print" << std::endl;

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
}

/* print content of solver */
void QPSolver_Global::printcontent(std::ostream &output) const {
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
void QPSolver_Global::printstatus(std::ostream &output) const {

	output <<  this->get_name() << std::endl;
	output <<  " - max(it):       " << this->it_last << std::endl;
	output <<  " - max(hessmult): " << this->hessmult_last << std::endl;
	output <<  " - max(fx):  " << this->fx << std::endl;	
	output <<  " - used memory:   " << MemoryCheck::get_virtual() << "%" << std::endl;
}


std::string QPSolver_Global::get_name() const {
	return "General Global QP Solver";
}

/* solve the problem */
void QPSolver_Global::solve() {
	if(setting.debug_mode >= 100) coutMaster << "(QPSolver_Global)FUNCTION: solve" << std::endl;

	/* for each process prepare QP solver and solve the problem */
	QPSolver<PetscVector> *solver_local; /* local solver */

	/* global data */
	Vec x_global = qpdata->get_x()->get_vector();
	Vec x0_global = qpdata->get_x0()->get_vector();
	Vec b_global = qpdata->get_b()->get_vector();

	/* local data */
	int local_size;
	Vec x_local;
	Vec x0_local;
	Vec b_local;

	/* allocate local data */
	TRY( VecGetLocalSize(x_global, &local_size) );
	TRY( VecCreateSeq(PETSC_COMM_SELF, local_size, &x_local) );	
	TRY( VecCreateSeq(PETSC_COMM_SELF, local_size, &x0_local) );	
	TRY( VecCreateSeq(PETSC_COMM_SELF, local_size, &b_local) );	

	/* get local vectors */
	TRY( VecGetLocalVector(x_global,x_local) );
	TRY( VecGetLocalVectorRead(x0_global,x0_local) );
	TRY( VecGetLocalVectorRead(b_global,b_local) );

	GeneralVector<PetscVector> x(x_local);
	GeneralVector<PetscVector> x0(x0_local);
	GeneralVector<PetscVector> b(b_local);

	/* get local vectors and prepare local data */
	QPData<PetscVector> data_local; /* data of inner cg solver */
	data_local.set_A(qpdata->get_A());
	data_local.set_b(&b);
	data_local.set_x(&x);
	data_local.set_x0(&x0);

	/* create new instance of local solver */
	solver_local = new QPSolver<PetscVector>(data_local);

	/* copy settings */
		
	/* solve local problem */
	solver_local->solve();

	/* restore global vectors */
	TRY( VecRestoreLocalVector(x_global,x_local) );
	TRY( VecRestoreLocalVectorRead(x0_global,x0_local) );
	TRY( VecRestoreLocalVectorRead(b_global,b_local) );

	/* destroy local storage */
//	TRY( VecDestroy(&x_local) );	
//	TRY( VecDestroy(&x0_local) );	
//	TRY( VecDestroy(&b_local) );	

	/* free local solver */
	free(solver_local); 
	
	//TODO:temp
//	TRY( VecView(x_global,PETSC_VIEWER_STDOUT_WORLD) );
	
}

double QPSolver_Global::get_fx() const {
	if(setting.debug_mode >= 100) coutMaster << "(QPSolver_Global)FUNCTION: get_fx()" << std::endl;

	return child_solver->get_fx(); // TODO: control existence
}

int QPSolver_Global::get_it() const {
	if(setting.debug_mode >= 100) coutMaster << "(QPSolver_Global)FUNCTION: get_it()" << std::endl;

	return child_solver->get_it(); // TODO: control existence
}

int QPSolver_Global::get_hessmult() const {
	if(setting.debug_mode >= 100) coutMaster << "(QPSolver_Global)FUNCTION: get_hessmult()" << std::endl;

	return child_solver->get_hessmult(); // TODO: control existence
}


//QPData<PetscVector> *QPSolver_Global::get_data() const {
//	return qpdata;
//}

} /* end namespace */

#endif
